//#include <iostream>
#include <boost/program_options.hpp>
#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
//#include "kindex.h"
//#include "alnstream.h"
#include "parallel.h"

/**
 * Parse the command line arguments. Checks all arguments for validity and prints a summary to the command line.
 * @param license License information string.
 * @param argc The number of arguments
 * @param argv The arguments.
 * @author Jakob Schulze
 */
int parseCommandLineArguments(std::string license, int argc, char const ** argv);

/**
 * Parse the -r / -reads argument to the respective fields in globalAlignmentSettings (seqLengths, isBarcode).
 * The argument is a vector of strings, each being a number and a subsequent B (Barcode) or R (Read).
 * The order defined in the argument is also the handled order in the program.
 * @param readsArg User argument containing the read type and length information
 * @return true, if successful (no parsing error and summed length equals the given number of cycles).
 * @author Tobias Loka
 */
bool parseReadsArgument(std::vector < std::string > readsArg);

/**
 * Parse the barcode information.
 * The argument is a vector of barcodes. Multi-Barcodes (Duplex) must be separated by a "-" (e.g., AACA-AGGT = 2x4bp).
 * The reads argument must be parsed BEFORE parsing the barcodes (call parseReadsArgument() first).
 * @param barcodeArg User argument containing the barcode information
 * @return true, if successful (no parsing error and correct lengths).
 * @author Tobias Loka
 */
bool parseBarcodeArgument(std::vector < std::string > barcodeArg);
