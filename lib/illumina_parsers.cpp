#include "illumina_parsers.h"


// Constructor takes filename and directly loads the whole file
uint64_t BclParser::open (std::string fname) {
  // read the whole file as a chunk
  data = read_binary_file(fname);
  // extract the number of reads
  memcpy(&num_reads,data.data(),4);
  // set the position pointer to the beginning of the data block
  position = 4;
  
  return data.size();
}

// Get the next base call
char BclParser::next() {
  if ( position < data.size() ) {
    position++;
    return *(data.data()+position-1);
  }
  else {
    throw std::runtime_error("Error reading BCL file: requested position is beyond EOF." );
  }
}

// Check if there are base calls left
bool BclParser::has_next() {
  return (position < data.size());
}

// Returns the total number of base calls in the file
uint32_t BclParser::size() {
  return num_reads;
}



// Constructor takes filename and directly loads the whole file
uint64_t FilterParser::open (std::string fname) {
  // read the whole file as a chunk
  data = read_binary_file(fname);
  // extract the number of reads
  memcpy(&num_reads,data.data()+8,4);
  // set the position pointer to the beginning of the data block
  position = 12;

  return data.size();
}
  
// Get the next filter flag
bool FilterParser::next() {
  if ( position < data.size() ) {
    position++;
    return (*(data.data()+position-1) > 0);
  }
  else {
    throw std::runtime_error("Error reading filter file: requested position is beyond EOF." );
  }
}

// Check if there are filter flags left
bool FilterParser::has_next() {
  return (position < data.size());
}

// Returns the total number of filter flags in the file
uint32_t FilterParser::size() {
  return num_reads;
}
