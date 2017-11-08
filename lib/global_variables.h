#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

#include "../lib/headers.h"

class AlignmentSettings;

extern AlignmentSettings globalAlignmentSettings;

extern mutex_map<std::string> fileLocks;

#endif
