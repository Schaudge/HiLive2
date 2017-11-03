#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/global_variables.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"

AlignmentSettings globalAlignmentSettings;

Task writeNextTaskToBam ( std::deque<AlnOut> & alnouts ) {

	// Search for the next task to write
	for ( auto& alnout : alnouts ) {

		// Only loop through non-finished output deques
		if ( !alnout.is_finished() ) {

			// Try to write the next task from the deque
			Task return_status = alnout.write_next();

			// Proceed with next deque if no task was written
			if ( return_status == NO_TASK ) {

				continue;

			}

			// Return the written task
			else {

				// TODO: this is not thread safe!
				if ( alnout.is_finished() ) {

					alnout.finalize();

				}

				return return_status;

			}

		}

	}

	return NO_TASK;
}

/**
 * Worker function for the alignment threads.
 * @param tasks Reference to the "to do" task queue
 * @param finished Reference to the "finished" task queue
 * @param failed Reference to the "failed" task queue
 * @param idx Pointer to the index object
 * @param surrender Control flag (threads stop if true)
 */
void worker (TaskQueue & tasks, TaskQueue & finished, TaskQueue & failed, KixRun* idx, std::deque<AlnOut> & alnouts, std::atomic<CountType> & writing_threads, bool & surrender ) {

    // Continue until surrender flag is set
    while ( !surrender ) {

    	// Start an output task if output threads and tasks available.
    	if ( ++writing_threads < globalAlignmentSettings.get_num_out_threads() ) {
    		Task written_task = writeNextTaskToBam( alnouts );
    		--writing_threads;
    		if ( written_task != NO_TASK )
    			continue;
    	}

        // Try to obtain a new task
        Task t = tasks.pop();

        // If "to do" task was found
        if ( t != NO_TASK ) {

            // Execute the task
            bool success = true;
            try {

                StreamedAlignment s (t.lane, t.tile, t.seqEl.length);
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

            	// Make previous cycle available for output. If current cycle is the last one, make current cycle available.
            	CountType seqCycle = getSeqCycle(t.cycle, t.seqEl.id);
            	CountType output_cycle = seqCycle == globalAlignmentSettings.get_cycles() ? seqCycle : seqCycle - 1;

            	if ( globalAlignmentSettings.is_output_cycle( output_cycle ) ) {
            		for ( auto& alnout : alnouts ) {
            			if ( output_cycle == alnout.get_cycle() ) {
            				alnout.set_task_available( Task(t.lane, t.tile, output_cycle));
            				break;
            			}
            		}
            	}

        		finished.push(t);

            }
            else {
                failed.push(t);
            }

        }

        // Thread is idle --> Also use it for output if the maximum number of output threads is exceeded.
        else {

    		writing_threads++;
        	writeNextTaskToBam( alnouts );
    		writing_threads--;

        }

        // send this thread to sleep for a second
        std::this_thread::sleep_for (std::chrono::milliseconds(100));

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
    index->deserialize_file(globalAlignmentSettings.get_index_fname());

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
                StreamedAlignment s (ln, tl, globalAlignmentSettings.getSeqByMate(mate).length);
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

    // Init output controller for each output cycle. TODO: check if it is possible to replace the deque by a map (cycle, alnout).
    std::deque<AlnOut> alnouts;
    for ( CountType cycle : globalAlignmentSettings.get_output_cycles() ) {
    	alnouts.emplace_back(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), cycle, index);
    }

    // Number of threads currently used for writing output.
    std::atomic<CountType> writing_threads(0);

    // Flag to stop the threads.
    bool surrender = false;

    // Create the threads
    std::cout << "Creating " << globalAlignmentSettings.get_num_threads() << " threads." << std::endl;
    std::vector<std::thread> workers;
    for (int i = 0; i < globalAlignmentSettings.get_num_threads(); i++) {
        workers.push_back(std::thread(worker, std::ref(toDoQ), std::ref(finishedQ), std::ref(failedQ), index, std::ref(alnouts), std::ref(writing_threads), std::ref(surrender)));
    }

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

    std::cout << "Finished all alignments." << std::endl;
    std::cout << "Waiting for output tasks..." << std::endl;

    for ( auto& alnout : alnouts ) {
    	while ( !alnout.is_finished() )
    		; // wait
    	alnout.finalize();
    }


    std::cout << "Finished output tasks." << std::endl;

    // Halt the threads
    surrender = true;
    for (auto& w : workers) {
        w.join();
    }

    std::cout << "All threads joined." << std::endl;
//    std::cout << "Total mapping time: " << time(NULL) - t_start << " s" << std::endl << std::endl;
    delete index;

    std::cout << "Total run time: " << time(NULL) - t_start << " s" << std::endl;
    exit(EXIT_SUCCESS);
}
