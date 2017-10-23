#include <boost/program_options.hpp>

#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"

namespace po = boost::program_options;

AlignmentSettings globalAlignmentSettings;

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

	index->get_header_information(globalAlignmentSettings.get_index_fname());
	index->store_kmer();

	std::cout << "Start writing ouput." << std::endl;

	for ( CountType cycle : globalAlignmentSettings.get_output_cycles() ) {
		try {
			AlnOut alnout(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), cycle, index);
			for ( auto& lane : globalAlignmentSettings.get_lanes() ) {
				for ( auto& tile : globalAlignmentSettings.get_tiles() ) {
					alnout.task_available( Task(lane, tile, cycle) );
				}
			}

			while ( !alnout.is_finished() ) {
				for ( auto i=0; i < alnout.get_task_status_num( BCL_AVAILABLE ); i++ ) {
					alnout.write_next();
				}
//				alnout.write_records();
		        std::this_thread::sleep_for (std::chrono::milliseconds(10000));
			}

			std::cout << "Num_tasks: " << std::to_string( alnout.tasks.size() ) << std::endl;
			std::cout << "Finished: " << std::to_string( alnout.get_task_status_num( FINISHED ) ) << std::endl;
			std::cout << "Failed: " << std::to_string( alnout.get_task_status_num( FAILED ) ) << std::endl;


			alnout.join();

//			if ( alignments_to_sam(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), index, cycle) )
//				std::cout << "Cycle " << std::to_string(cycle) << " ... " << "success." << std::endl;
//			else
//				std::cout << "Cycle " << std::to_string(cycle) << " ... " << "failed." << std::endl;
		} catch ( std::exception & ex ) {
			std::cout << "Cycle " << std::to_string(cycle) << " ... " << "failed: " << std::endl << ex.what() << std::endl << std::endl;
		}
	}

	std::cout << "Finished." << std::endl;

	delete index;

	return EXIT_SUCCESS;
}
