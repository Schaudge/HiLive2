#include <iostream>
#include <boost/program_options.hpp>
#include "headers.h"
#include "definitions.h"
#include "kindex.h"
#include "alnstream.h"
#include "parallel.h"

int parseCommandLineArguments(AlignmentSettings & settings, std::string license, int argc, char const ** argv);

/**
 * Parse the -r / -reads argument to the respective fields in AlignmentSettings (seqLengths, isBarcode).
 * @param settings Object containing the program settings
 * @param readsArg User argument containing the read type and length information
 * @return true, if successful (no parsing error and summed length equals the given number of cycles).
 * @author Tobias Loka
 */
bool parseReadsArgument(AlignmentSettings & settings, std::vector < std::string > readsArg);
