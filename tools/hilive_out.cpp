#include <boost/program_options.hpp>

#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"

namespace po = boost::program_options;

AlignmentSettings globalAlignmentSettings;
mutex_map<std::string> fileLocks;

/**
 * Main function that organizes the overall structure of the program.
 * @param argc Number of arguments
 * @param argv Argument array
 * @return 0 on success, other numbers on error
 */
int main(int argc, const char* argv[]) {

	// Program start output
	std::cout << std::endl << "------" << std::endl << "HiLive Output Tool v"<< HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR <<
			" - Output of Realtime Alignments of Illumina Reads" << std::endl << "------" << std::endl<< std::endl;

	// Parse the command line arguments
    HiLiveOutArgumentParser argumentParser(argc, argv);
	int parser_returnStatus = argumentParser.parseCommandLineArguments();

	// Successful execution of "help" or "license"
	if ( parser_returnStatus == 1 ) {
		exit(EXIT_SUCCESS);
	}

	// Parsing error
	else if ( parser_returnStatus == -1 ) {
		std::cout << "Parsing of command line options failed. For help, type 'hilive-out --help'." << std::endl;
		exit(EXIT_FAILURE);
	}

	// load the index
	std::cout << "Loading Index Header..." << std::endl;
	KixRun* index = new KixRun();
//
//	index->get_header_information(globalAlignmentSettings.get_index_fname());
//	index->store_kmer();

    index->load_metadata( globalAlignmentSettings.get_index_fname() );
    index->load_fmindex( globalAlignmentSettings.get_index_fname() );

	std::cout << "Start writing ouput." << std::endl;

	// Maximum number of output threads.
	CountType max_output_threads = std::max(CountType(1), globalAlignmentSettings.get_num_out_threads());

	std::cout << "Using " << max_output_threads << " threads." << std::endl;

	std::deque<AlnOut> alnouts;
	std::deque<std::thread> threads;

	bool all_finished = false;

	for ( CountType cycle : globalAlignmentSettings.get_output_cycles() ) {
		alnouts.emplace_back(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), cycle, index);
	}

	while ( !all_finished ) {

		all_finished = true;
		CountType num_active_threads = 0;

		for ( auto& alnout : alnouts ) {
			num_active_threads += alnout.get_task_status_num( RUNNING );
		}

		for ( auto& alnout : alnouts ) {

			// Assume that all tasks are available.
			for ( auto& lane : globalAlignmentSettings.get_lanes() ) {
				for ( auto& tile : globalAlignmentSettings.get_tiles() ) {
					alnout.set_task_available( Task(lane, tile, alnout.get_cycle()) );
				}
			}

			if ( !alnout.is_finished() ) {
				all_finished = false;

				if ( num_active_threads < max_output_threads ) {
					CountType newThreads = max_output_threads - num_active_threads;
					for ( auto i=0; i < newThreads; i++) {
							threads.emplace_back(&AlnOut::write_next, &alnout);
							++num_active_threads;
					}
				}
			} else if ( !alnout.is_finalized() ) {
				alnout.finalize();
				std::cout << "Finished output of cycle " << alnout.get_cycle() << " (" << alnout.get_task_status_num( FINISHED ) << " finished, " << alnout.get_task_status_num( FAILED ) << " failed)." << std::endl;
			}
		}
		std::this_thread::sleep_for (std::chrono::milliseconds(1000));
	}

    // Clear the vector will destruct all elements.
	alnouts.clear();

	// Ensure that all threads are finished.
	for ( auto& thread : threads )
		thread.join();

	std::cout << "Finished." << std::endl;

	delete index;

	return EXIT_SUCCESS;
}
