#include <iostream>
#include <boost/program_options.hpp>
#include "headers.h"
#include "definitions.h"
#include "kindex.h"
#include "alnstream.h"
#include "parallel.h"

int parseCommandLineArguments(AlignmentSettings & settings, std::string license, int argc, char const ** argv);
