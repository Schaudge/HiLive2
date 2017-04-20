#include <iostream>
#include <boost/program_options.hpp>
#include "headers.h"
#include "definitions.h"
#include "kindex.h"
#include "alnstream.h"
#include "parallel.h"
#include "tinyxml2.h"

/**
 * Parse the command line arguments. Checks all arguments for validity and prints a summary to the command line.
 * @param settings Object containing the program settings.
 * @param license License information string.
 * @param argc The number of arguments
 * @param argv The arguments.
 * @author Jakob Schulze
 */
int parseCommandLineArguments(AlignmentSettings & settings, std::string license, int argc, char const ** argv);

/**
 * Parse the -r / -reads argument to the respective fields in AlignmentSettings (seqLengths, isBarcode).
 * The argument is a vector of strings, each being a number and a subsequent B (Barcode) or R (Read).
 * The order defined in the argument is also the handled order in the program.
 * @param settings Object containing the program settings
 * @param readsArg User argument containing the read type and length information
 * @return true, if successful (no parsing error and summed length equals the given number of cycles).
 * @author Tobias Loka
 */
bool parseReadsArgument(AlignmentSettings & settings, std::vector < std::string > readsArg);

/**
 * Parse the barcode information.
 * The argument is a vector of barcodes. Multi-Barcodes (Duplex) must be separated by a "-" (e.g., AACA-AGGT = 2x4bp).
 * The reads argument must be parsed BEFORE parsing the barcodes (call parseReadsArgument() first).
 * @param settings Object containing the program settings
 * @param barcodeArg User argument containing the barcode information
 * @return true, if successful (no parsing error and correct lengths).
 * @author Tobias Loka
 */
bool parseBarcodeArgument(AlignmentSettings & settings, std::vector < std::string > barcodeArg);
