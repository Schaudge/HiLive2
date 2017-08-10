#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/global_variables.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"

std::string license =
"Copyright (c) 2015-2016, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info.\n"
"All rights reserved.\n"
"\n"
"Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:\n"
"\n"
"1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.\n"
"\n"
"2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.\n"
"\n"
"3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.\n"
"\n"
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.";

// the worker function for the threads
void worker (TaskQueue & tasks, TaskQueue & finished, TaskQueue & failed, KixRun* idx, bool & surrender ) {

    // loop that keeps on running until the surrender flag is set
    while ( !surrender ) {

        // try to obtain a new task
        Task t = tasks.pop();
        if ( t != NO_TASK ) {

            // Execute the task
            bool success = true;
            try {

                StreamedAlignment s (t.lane, t.tile, t.seqEl.length);
                uint64_t num_seeds;

                // Seed extension if current read is no barcode.
                if ( !t.seqEl.isBarcode() ) {
                	num_seeds = s.extend_alignment(t.cycle,t.seqEl.id,t.seqEl.mate,idx);
                	std::cout << "Task [" << t << "]: Found " << num_seeds << " seeds." << std::endl;

                // If current read is barcode, extend barcode sequence in all sequence read align files
                } else {
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


AlignmentSettings globalAlignmentSettings;

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

	boost::property_tree::ptree xml_out = globalAlignmentSettings.to_ptree();
	if (! write_xml(xml_out, get_xml_out_name()))
		exit(EXIT_FAILURE);

    // Load the index
    std::cout << "Loading Index ... " << std::endl;
    KixRun* index = new KixRun();
    index->deserialize_file(globalAlignmentSettings.get_index_fname());

    // Report k-mer structure that was retrieved from the index
    std::cout << "K-mer weight:        " << unsigned(globalAlignmentSettings.get_kmer_weight()) << std::endl;
    std::cout << "K-mer span:          " << unsigned(globalAlignmentSettings.get_kmer_span()) << std::endl;
    std::cout << "K-mer gap positions: ";

      if ( globalAlignmentSettings.get_kmer_gaps().size() > 0 ) {
    	  for ( auto pos : globalAlignmentSettings.get_kmer_gaps() ) {
    		  if ( pos != *(globalAlignmentSettings.get_kmer_gaps().begin()) )
    			  std::cout << ",";
    		  std::cout << (uint16_t) pos;
    	  }
    	  std::cout << std::endl;
      } else {
    	  std::cout << "-" << std::endl;
      }
      std::cout << std::endl;

    // Create the overall agenda
    Agenda agenda (globalAlignmentSettings.get_cycles(), globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles());

    // Prepare the alignment
    std::cout << "Initializing Alignment files. Waiting for the first cycle to finish." << std::endl;
    bool first_cycle_available = false;

    // Wait for the first cycle to be written. Attention - this loop will wait infinitely long if no first cycle is found
    while ( !first_cycle_available ) {

        // Check for new BCL files and update the agenda status
        agenda.update_status();

        // Check if the first cycle is available for all tiles
        first_cycle_available = true;
        for ( auto ln : globalAlignmentSettings.get_lanes() ) {
            for ( auto tl : globalAlignmentSettings.get_tiles() ) {
                if ( agenda.get_status(Task(ln,tl,globalAlignmentSettings.getSeqById(0),1)) != BCL_AVAILABLE) {
                    first_cycle_available = false;
                }
            }
        }

        // Take a small break
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    std::cout << "First cycle complete. Starting alignment." << std::endl;

    // Write empty alignment file for each tile and for each sequence read
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

    if (globalAlignmentSettings.get_temp_dir() != "" && !is_directory(globalAlignmentSettings.get_temp_dir())){
        std::cerr << "Error: Could not find temporary directory " << globalAlignmentSettings.get_temp_dir() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Set up the queues
    TaskQueue toDoQ;
    TaskQueue finishedQ;
    TaskQueue failedQ;

    // Create the threads
    std::cout << "Creating " << globalAlignmentSettings.get_num_threads() << " threads." << std::endl;
    bool surrender = false;
    std::vector<std::thread> workers;
    for (int i = 0; i < globalAlignmentSettings.get_num_threads(); i++) {
        workers.push_back(std::thread(worker, std::ref(toDoQ), std::ref(finishedQ), std::ref(failedQ), index, std::ref(surrender)));
    }

    // Process all tasks on the agenda
    while ( !agenda.finished() ) {
        // check for new BCL files and update the agenda status
        agenda.update_status();

        // fill the ToDo queue with tasks from the agenda
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

    std::cout << "All threads joined." << std::endl;
    std::cout << "Total mapping time: " << time(NULL) - t_start << " s" << std::endl;
    std::cout << "Writing output file." << std::endl;
    alignments_to_sam(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), index, globalAlignmentSettings.get_cycles());
    delete index;

    if ( globalAlignmentSettings.get_trimmedReads().size() > 0 ) {
    	std::cout << "Trimmed reads: " ;
    	for ( auto tr : globalAlignmentSettings.get_trimmedReads() ) {
    		std::cout << tr << ", ";
    	}
    }

    std::cout << std::endl;
    std::cout << "Total run time: " << time(NULL) - t_start << " s" << std::endl;
    exit(EXIT_SUCCESS);
}
