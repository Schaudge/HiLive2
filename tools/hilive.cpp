#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/global_variables.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"

AlignmentSettings globalAlignmentSettings;

/**
 * Worker function for the alignment threads.
 * @param tasks Reference to the "to do" task queue
 * @param finished Reference to the "finished" task queue
 * @param failed Reference to the "failed" task queue
 * @param idx Pointer to the index object
 * @param surrender Control flag (threads stop if true)
 */
void worker (TaskQueue & tasks, TaskQueue & finished, TaskQueue & failed, KixRun* idx, bool & surrender ) {

    // Continue until surrender flag is set
    while ( !surrender ) {

        // Try to obtain a new task
        Task t = tasks.pop();

        // If "to do" task was found
        if ( t != NO_TASK ) {

            // Execute the task
            bool success = true;
            try {

                StreamedAlignment s (t.lane, t.tile, globalAlignmentSettings.get_root(), t.seqEl.length);
                uint64_t num_seeds;

                // Seed extension if current read is sequence fragment.
                if ( !t.seqEl.isBarcode() ) {
                	num_seeds = s.extend_alignment(t.cycle,t.seqEl.id,t.seqEl.mate,idx);
                	std::cout << "Task [" << t << "]: Found " << num_seeds << " seeds." << std::endl;

                }

                // Barcode extension if current read is barcode fragment
                else {
                	CountType mate = 1;
                	for ( ; mate <= globalAlignmentSettings.get_mates(); mate++ ) {
                		SequenceElement seqEl = globalAlignmentSettings.getSeqByMate(mate);
                		CountType current_mate_cycle = t.seqEl.id < seqEl.id ? 0 : seqEl.length;
                		s.extend_barcode(t.cycle, current_mate_cycle, t.seqEl.id, mate);
                	}
                	std::cout << "Task [" << t << "]: Extended barcode of " << --mate << " mates." << std::endl;
                }
            }
            catch (const std::exception &e) {
                std::cerr << "Failed to finish task [" << t << "]: " << e.what() << std::endl;
                success = false;
            }

            // Push the task in the correct Task Queue (Finished or Failed)
            if (success) {
                finished.push(t);
            }
            else {
                failed.push(t);
            }

        }
        else {
            // send this thread to sleep for a second
            std::this_thread::sleep_for (std::chrono::milliseconds(100));
        }
    }  
}

/**
 * Worker function for the output thread.
 * Should only be called by one single thread!
 * @param agenda Reference to the agenda that organizes the tasks
 * @param idx Pointer to the index object
 * @param surrender Control flag (threads stop if true)
 * @author Tobias Loka
 */
void output_worker( Agenda & agenda, KixRun* idx, bool & surrender ) {

	// Get the output cycles (must be sorted!)
	std::vector<CountType> output_cycles = globalAlignmentSettings.get_output_cycles();

	// Delayed surrender
	bool alignments_finished = false;

	// Continue as long as there are output cycles left
	while ( output_cycles.size() > 0 ) {

		alignments_finished = surrender;

		// Create output for the next cycle if all related tasks are finished.
		CountType next_cycle = *(output_cycles.begin());
		if ( agenda.finished(next_cycle) ) {

			try {
				if ( !alignments_to_sam(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), idx, next_cycle) )
					std::cerr << "Writing output for cycle " << std::to_string(next_cycle) << " failed." << std::endl;
				else {
					std::cout << "Wrote output for cycle " << std::to_string(next_cycle) << "." << std::endl;
				}
			}
			catch ( std::exception & e ) {
				std::cerr << "Writing output for cycle " << std::to_string(next_cycle) << " failed: " << e.what() << std::endl;
			}

			// Remove current cycle from the list of output cycles.
			output_cycles.erase(output_cycles.begin());
		}

		// If the the surrender flag set but no cycle was handled, remove the cycle from the list
		if ( alignments_finished && output_cycles.size() > 0 && next_cycle == *(output_cycles.begin()) ) {
			std::cerr << "Writing output for cycle " << std::to_string(next_cycle) << " failed: Not all tasks finished." << std::endl;
			output_cycles.erase(output_cycles.begin());
		}

		// Sleep a second
        std::this_thread::sleep_for (std::chrono::milliseconds(1000));

	}

}

/**
 * Main function that organizes the overall structure of the program.
 * @param argc Number of arguments
 * @param argv Argument array
 * @return 0 on success, other numbers on error
 */
int main(int argc, const char* argv[]) {

	// Variable for runtime measurement
    time_t t_start = time(NULL);

    // Program start output
	std::cout << std::endl << "------" << std::endl << "HiLive v"<< HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR <<
			" - Realtime Alignment of Illumina Reads" << std::endl << "------" << std::endl<< std::endl;

    // Parse command line arguments
    HiLiveArgumentParser argumentParser(argc, argv);
	int parser_returnStatus = argumentParser.parseCommandLineArguments();

	// Successful execution of "help" or "license"
	if ( parser_returnStatus == 1 ) {
		exit(EXIT_SUCCESS);
	}

	// Parsing error
	else if ( parser_returnStatus == -1 ) {
		std::cout << "Parsing of command line options failed. For help, type 'hilive --help'." << std::endl;
		exit(EXIT_FAILURE);
	}

    // Load the index
    std::cout << "Loading Index ... " << std::endl;
    KixRun* index = new KixRun();

    index->load_metadata( globalAlignmentSettings.get_index_fname() );
    index->load_fmindex( globalAlignmentSettings.get_index_fname() );

    // Report loaded k-mer properties
    std::cout << std::endl;
    std::cout << "kmer span:   " << std::to_string(globalAlignmentSettings.get_kmer_span()) << std::endl;
    std::cout << "kmer weight: " << std::to_string(globalAlignmentSettings.get_kmer_weight()) << std::endl;
    std::cout << "kmer gaps:   ";
    for ( auto gap : globalAlignmentSettings.get_kmer_gaps() ) {
    	std::cout << gap << " ";
    }
    std::cout << std::endl << std::endl;

  	// Write the alignment settings to an XML file
  	boost::property_tree::ptree xml_out = globalAlignmentSettings.to_ptree();
  	if ( ! write_xml(xml_out, get_settings_name()) )
  		exit(EXIT_FAILURE);

    // Create the overall agenda
    Agenda agenda (globalAlignmentSettings.get_cycles(), globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles());

    // Wait for the first cycle to be written
    std::cout << "Waiting for the first cycle to finish..." << std::endl;
    while ( ! agenda.cycle_available(1) ) {
        agenda.update_status();
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    // Write empty alignment file for each tile and for each sequence read
    std::cout << "Initializing Alignment files..." << std::endl;
    for (uint16_t ln : globalAlignmentSettings.get_lanes()) {
        for (uint16_t tl : globalAlignmentSettings.get_tiles()) {
            CountType mate = 1;
            for ( ; mate <= globalAlignmentSettings.get_mates(); mate++ ) {
                StreamedAlignment s (ln, tl, globalAlignmentSettings.get_root(), globalAlignmentSettings.getSeqByMate(mate).length);
                s.create_directories();
                s.init_alignment(mate);
            }
        }
    }

    std::cout << "First cycle complete. Starting alignment." << std::endl;



    // Set up the queues
    TaskQueue toDoQ;
    TaskQueue finishedQ;
    TaskQueue failedQ;

    // Create the alignment threads
    std::cout << "Creating " << globalAlignmentSettings.get_num_threads() << " threads." << std::endl;
    bool surrender = false;
    std::vector<std::thread> workers;
    for (int i = 0; i < globalAlignmentSettings.get_num_threads(); i++) {
        workers.push_back(std::thread(worker, std::ref(toDoQ), std::ref(finishedQ), std::ref(failedQ), index, std::ref(surrender)));
    }

    // Create the output thread
    std::thread output_thread = std::thread( output_worker, std::ref(agenda), index, std::ref(surrender) );

    // Process all tasks on the agenda
    while ( !agenda.finished() ) {

        // check for new BCL files and update the agenda status
        agenda.update_status();

        // fill the To Do queue with tasks from the agenda
        while(true) {
            Task t = agenda.get_task();
            if (t == NO_TASK)
                break;
            toDoQ.push(t);
            agenda.set_status(t,RUNNING);
        }

        // take a look in the finished queue and process finished tasks
        while(true) {
            Task t = finishedQ.pop();
            if (t == NO_TASK)
                break;
            agenda.set_status(t,FINISHED);
        }

        // take a look in the failed queue and process failed tasks
        while(true) {
            Task t = failedQ.pop();
            if (t == NO_TASK)
                break;
            if (agenda.get_status(t) == RUNNING) {
                // give it one more chance
                agenda.set_status(t,RETRY);
                toDoQ.push(t);
            }
            else {
                agenda.set_status(t,FAILED);
                std::cout << "Task failed! " << t << std::endl;
            }
        }

        // take a small break
        std::this_thread::sleep_for (std::chrono::milliseconds(100));
    }  

    // Halt the threads
    surrender = true;
    for (auto& w : workers) {
        w.join();
    }
    output_thread.join();

    std::cout << "All threads joined." << std::endl;
    std::cout << "Total mapping time: " << time(NULL) - t_start << " s" << std::endl << std::endl;
    delete index;

    std::cout << "Total run time: " << time(NULL) - t_start << " s" << std::endl;
    exit(EXIT_SUCCESS);
}
