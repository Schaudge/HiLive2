#include <boost/program_options.hpp>

#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/kindex.h"
#include "../lib/argument_parser.h"

AlignmentSettings globalAlignmentSettings;
mutex_map<std::string> fileLocks;


/**
 * Main function that organizes the overall structure of the program.
 * @param argc Number of arguments
 * @param argv Argument array
 * @return 0 on success, other numbers on error
 */
int main(int argc, const char* argv[]) {

	// Program start output.
	std::cout << std::endl << "________________________________________________________________________________" << std::endl << std::endl << "HiLive Index Builder v" << HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR <<
			" - Build Index for Realtime Alignment of Illumina Reads" << std::endl << "________________________________________________________________________________" << std::endl << std::endl;

	// Init argument parser
	BuildIndexArgumentParser argumentParser(argc, argv);

	// Parse command line arguments
	int parser_returnStatus = argumentParser.parseCommandLineArguments();

	// Successful execution of "help" or "license"
	if ( parser_returnStatus == 1 ) {
		exit(EXIT_SUCCESS);
	}

	// Parsing error
	else if ( parser_returnStatus == -1 ) {
		std::cout << "Parsing of command line options failed. For help, type 'hilive-build --help'." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Build index
	std::cout << "Create FM-index from file " << argumentParser.fasta_name << " ..." << std::endl;

	  KixBuild* index = new KixBuild();

	// create and save the index
	  switch( index->create_index(argumentParser.fasta_name, argumentParser.index_name, !argumentParser.do_not_convert_spaces, argumentParser.trim_ids) ) {
	  case 0: std::cerr << "Index was successfully built to file(s) " << argumentParser.index_name << std::endl; break;
	  case 1: std::cerr << "Reading input file failed (no sequences found)." << std::endl; break;
	  case 2: std::cerr << "Reading input file failed." << std::endl; break;
	  case 3: std::cerr << "Writing metadata to file(s) failed." << std::endl; break;
	  case 4: std::cerr << "Index building failed." << std::endl; break;
	  case 5: std::cerr << "Writing index to file " << argumentParser.index_name << " failed." << std::endl; break;
	  default: std::cerr << "Failed due to an unknown error." << std::endl; break;
	  }

	delete index;

	return EXIT_SUCCESS;
} 
