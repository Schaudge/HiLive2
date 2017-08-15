#include <boost/program_options.hpp>

#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"

namespace po = boost::program_options;

AlignmentSettings globalAlignmentSettings;

//// create a special SAM worker, that writes out a SAM file for a tile
//void sam_worker (TaskQueue & tasks, AlignmentSettings* settings, KixRun* idx) {
//
//    // loop that keeps on running until the surrender flag is set
//    while ( true ) {
//        if (tasks.size() == 0) return;
//		else {
//			// try to obtain a new task
//			Task t = tasks.pop();
//			// Execute the task
//			uint16_t test = t.cycle;
//			alignments_to_sam(t.lane,t.tile,t.root, test, idx, settings);
//        }
//    }
//}
//
//std::vector<Task> generate_sam_task(std::vector<uint16_t> lanes, std::vector<uint16_t> tiles, std::vector<uint16_t> cycles, CountType rlen, std::string root){
//	std::vector<Task> tv;
//	 // find tiles that are completely mapped
//	 for (uint16_t ln_id = 0; ln_id < lanes.size(); ++ln_id) {
//	   for (uint16_t tl_id = 0; tl_id < tiles.size(); ++tl_id) {
//		   for (uint16_t cl_id = 0; cl_id < cycles.size(); ++cl_id) {
//			 tv.push_back(Task(lanes[ln_id],tiles[tl_id],cycles[cl_id],rlen,root));
//		   }
//	   }
//	 }
//	return tv;
//}

int main(int argc, const char* argv[]) {

	// parse the command line arguments, store results in settings
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
			if ( alignments_to_sam(globalAlignmentSettings.get_lanes(), globalAlignmentSettings.get_tiles(), index, cycle) )
				std::cout << "Cycle " << std::to_string(cycle) << " ... " << "success." << std::endl;
			else
				std::cout << "Cycle " << std::to_string(cycle) << " ... " << "failed." << std::endl;
		} catch ( std::exception & ex ) {
			std::cout << "Cycle " << std::to_string(cycle) << " ... " << "failed: " << std::endl << ex.what() << std::endl << std::endl;
		}
	}

	std::cout << "Finished." << std::endl;

	delete index;
}
