#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

#include "../lib/headers.h"

class AlignmentSettings;
class KixRun;

extern AlignmentSettings globalAlignmentSettings;
extern KixRun* idx;

extern mutex_map<std::string> fileLocks;

#endif
