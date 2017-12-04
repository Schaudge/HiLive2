#ifndef ILLUMINA_PARSERS_H
#define ILLUMINA_PARSERS_H

#include "headers.h"
#include "definitions.h"
#include "tools.h"

// BCL file parser
class BclParser {
  // storage for the raw binary data
  std::vector<char> data;
  // current position in data
  uint32_t position;
  // number of reads in this bcl file
  uint32_t num_reads;
 public:
  // open file and directly load all data
  uint64_t open(std::string fname);
  
  // Get the next base call
  char next();
  
  // Check if there are base calls left
  bool has_next();

  // Returns the total number of base calls in the file
  uint32_t size();
};


// filter file parser
class FilterParser {
  // storage for the raw binary data
  std::vector<char> data;
  // current position in data
  uint32_t position;
  // number of reads in this filter file
  uint32_t num_reads = 0;
 public:
  // open file and directly load all data
  uint64_t open(std::string fname);
  
  // Get the next filter flag
  bool next();
  
  // Check if there are filter flags left
  bool has_next();

  // Returns the total number of filter flags in the file
  uint32_t size();
};


// clocs file parser




#endif /* ILLUMINA_PARSERS_H */
