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
	std::cout << std::endl << "------" << std::endl << "HiLive Index Builder v" << HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR <<
			" - Build Index for Realtime Alignment of Illumina Reads" << std::endl << "------" << std::endl << std::endl;

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
	std::cout << "Create index from file " << argumentParser.fasta_name << " ..." << std::endl;
	KixBuild* index = new KixBuild();
	index->add_fasta(argumentParser.fasta_name, !argumentParser.do_not_convert_spaces, argumentParser.trim_ids);

	// Trim index
	if ( argumentParser.trim > 0) {
		uint64_t trimmed = index->trim(argumentParser.trim);
		std::cout << "Removed " << trimmed << " k-mer positions from the database." << std::endl;
	}

	// Write index to file
	std::cout << "Writing index to file " << argumentParser.index_name << std::endl;
	index->serialize_file(argumentParser.index_name);

	delete index;

	return EXIT_SUCCESS;
} 
