#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/global_variables.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/alnout.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"
#include "../lib/tools_static.h"

AlignmentSettings globalAlignmentSettings;
KixRun* idx;
mutex_map<std::string> fileLocks;

/**
 * If a thread is used for output, call this function to start the next available output task.
 * @param alnout The deque of output controllers for each output cycle.
 * @return The written task. NO_TASK if no task was written (e.g., when no task is available).
 */
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

				// Finalize the output for this cycle if this was the last task.
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
void worker (TaskQueue & tasks, TaskQueue & finished, TaskQueue & failed, std::deque<AlnOut> & alnouts, std::atomic<CountType> & writing_threads, bool & surrender ) {

    // Continue until surrender flag is set
    while ( !surrender ) {

    	{ // scope for block guard
    		atomic_increment_guard<CountType> block( writing_threads );

    		// Start an output task if output threads and tasks available. Allow at least 1 output thread.
    		if ( block.get_incremented_value() <= globalAlignmentSettings.get_num_out_threads() || block.get_incremented_value() == 1 ) {
    			Task written_task = writeNextTaskToBam( alnouts );
    			if ( written_task != NO_TASK ) {
    				continue;
    			}
    		}
    	}

        // Try to obtain a new task
        Task t = tasks.pop();

        // If "to do" task was found
        if ( t != NO_TASK ) {

            // Execute the task
            bool success = true;
            std::stringstream ss;

            try {

                StreamedAlignment s (t.lane, t.tile, t.seqEl.length);
                uint64_t num_seeds;

                // Seed extension if current read is sequence fragment.
                if ( !t.seqEl.isBarcode() ) {
                	num_seeds = s.extend_alignment(t.cycle,t.seqEl.id,t.seqEl.mate);
                	ss << "Task [" << t << "]: Found " << num_seeds << " seeds." << std::endl;

                }

                // Barcode extension if current read is barcode fragment
                else {
                	CountType mate = 1;
                	for ( ; mate <= globalAlignmentSettings.get_mates(); mate++ ) {
                		SequenceElement seqEl = globalAlignmentSettings.get_seq_by_mate(mate);
                		CountType current_mate_cycle = t.seqEl.id < seqEl.id ? 0 : seqEl.length;
                		s.extend_barcode(t.cycle, current_mate_cycle, t.seqEl.id, mate);
                	}
                	ss << "Task [" << t << "]: Extended barcode of " << --mate << " mates." << std::endl;
                }

            	std::cout << ss.str();

            }
            catch (const std::exception &e) {
                ss << "Failed to finish task [" << t << "]: " << e.what() << std::endl;
                std::cerr << ss.str();
                success = false;
            }

            // Push the task in the correct Task Queue (Finished or Failed)
            if (success) {

            	// Make previous cycle available for output.
            	// If current cycle is the last one of the segment and this is an output cycle
            	// or it is the very last cycle, make current cycle available.
            	std::vector<CountType> output_cycles;
            	CountType seqCycle = getSeqCycle(t.cycle, t.seqEl.id);

            	// Make previous cycle available if it is an output cycle and > 1 (to prevent double output of the last cycle of a segment)
            	if ( t.cycle > 1 && globalAlignmentSettings.is_output_cycle( seqCycle - 1 ) )
            		output_cycles.push_back(seqCycle-1);

            	// Make current cycle available if it is the very last cycle or the last segment cycle and an output cycle.
            	if ( seqCycle == globalAlignmentSettings.get_cycles() || ( globalAlignmentSettings.is_output_cycle( seqCycle ) && t.cycle == t.seqEl.length ) )
            		output_cycles.push_back(seqCycle);

            	// Actually make the output cycles available.
            	for ( auto cycle : output_cycles ) {
            		for ( auto& alnout : alnouts ) {
            			alnout.set_task_available( Task(t.lane, t.tile, cycle));
            		}
            	}

        		finished.push(t);

            }
            else {
                failed.push(t);
            }

        }

        else {
            // send this thread to sleep for a second if (and only if) no task was available
            std::this_thread::sleep_for (std::chrono::milliseconds(1000));
        }

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
	std::cout << std::endl << "__________________________________________________" << std::endl << std::endl << "HiLive v"<< HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR <<
			" - Realtime Alignment of Illumina Reads" << std::endl << "__________________________________________________" << std::endl<< std::endl;

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
    idx = new KixRun();

    idx->load_metadata( globalAlignmentSettings.get_index_fname() );
    idx->load_fmindex( globalAlignmentSettings.get_index_fname() );

  	// Write the alignment settings to an XML file
  	boost::property_tree::ptree xml_out = globalAlignmentSettings.to_ptree();
  	if ( ! write_ini(xml_out, get_config_fname()) )
  		exit(EXIT_FAILURE);

    // Create the overall agenda
    Agenda agenda (globalAlignmentSettings.get_cycles(), globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), globalAlignmentSettings.get_start_cycle());


    // Wait for the first cycle to be written
    std::cout << "Waiting for the first cycle to finish..." << std::endl;
    while ( ! agenda.cycle_available(globalAlignmentSettings.get_start_cycle()) ) {
        agenda.update_status();
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    // Write empty alignment file for each tile and for each sequence read
    std::cout << "Initializing Alignment files..." << std::endl;
    for (uint16_t ln : globalAlignmentSettings.get_lanes()) {
        for (uint16_t tl : globalAlignmentSettings.get_tiles()) {
            CountType mate = 1;
            for ( ; mate <= globalAlignmentSettings.get_mates(); mate++ ) {

            	// Don't init files if "--continue" was used to start in a later cycle.
            	if ( getMateCycle(mate, globalAlignmentSettings.get_start_cycle()) > 1 )
            			continue;


                StreamedAlignment s (ln, tl, globalAlignmentSettings.get_seq_by_mate(mate).length);
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
    	if ( cycle >= globalAlignmentSettings.get_start_cycle() )
    		alnouts.emplace_back(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), cycle);
    }

    // Number of threads currently used for writing output.
    std::atomic<CountType> writing_threads(0);

    // Flag to stop the threads.
    bool surrender = false;

    // Create the threads
    std::cout << "Creating " << globalAlignmentSettings.get_num_threads() << " threads." << std::endl;
    std::vector<std::thread> workers;
    for (int i = 0; i < globalAlignmentSettings.get_num_threads(); i++) {
        workers.push_back(std::thread(worker, std::ref(toDoQ), std::ref(finishedQ), std::ref(failedQ), std::ref(alnouts), std::ref(writing_threads), std::ref(surrender)));
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
    	while ( !alnout.is_finalized() ) {
    		; // wait
    	}
    }

    // Clear the vector will destruct all elements.
    alnouts.clear();

    std::cout << "Finished output tasks." << std::endl;

    // Halt the threads
    surrender = true;
    for (auto& w : workers) {
        w.join();
    }

    std::cout << "All threads joined." << std::endl;
//    std::cout << "Total mapping time: " << time(NULL) - t_start << " s" << std::endl << std::endl;
    delete idx;

    std::cout << "Total run time: " << time(NULL) - t_start << " s" << std::endl;
    exit(EXIT_SUCCESS);
}
