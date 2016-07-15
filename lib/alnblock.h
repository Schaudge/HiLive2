#ifndef ALNBLOCK_H
#define ALNBLOCK_H

#include "headers.h"
#include "definitions.h"
#include "kindex.h"
#include "tools.h"
#include "alnread.h"


//-------------------------------------------------------------------//
//------  The Dataset-Alignment class  ------------------------------//
//-------------------------------------------------------------------//

class DatasetAlignment {

  // dataset information
  uint16_t cycle;
  uint16_t lane;
  uint16_t tile;
  std::string root; // the BaseCalls directory
  CountType rlen;

 public:

  // getter functions for the dataset information
  uint16_t get_cycle();
  uint16_t get_lane();
  uint16_t get_tile();
  std::string get_root();
  CountType get_read_length();

  void set_cycle(uint16_t c) {cycle=c;}
  
  // get the path to the bcl file of a given cycle. cl = 0 denotes the current cycle.
  std::string get_bcl_file(uint16_t cl = 0);

  // get the path to the alignment file. The alignment file is located in
  // <root>/L00<lane>/s_<lane>_<tile>.align
  std::string get_alignment_file();

  // get the path to the filter file. The illumina filter information is located in
  // <root>/L00<lane>/s_<lane>_<tile>.filter
  std::string get_filter_file();

  // illumina validity flags for all reads
  /* std::vector<bool> flags; */

  // all read alignments
  std::vector<ReadAlignment> reads;

  // serialize the object into a char vector
  std::vector<char> serialize();

  // serialize the object into a binary file
  uint64_t serialize_file(std::string f  = "");

  // deserialize (read) data from a char vector
  uint64_t deserialize(char* d);

  // deserialize data from a binary file
  uint64_t deserialize_file(std::string f = "");

  // initialize the dataset alignment with
  // 1. BaseCalls directory
  // 2. Lane number
  // 3. Tile number
  // 4. Read length of Read 1
  void initialize(std::string rt, int ln, int tl, int rl);

  // extend alignment by one nucleotide
  void add_nucleotide(KixRun* index, AlignmentSettings* settings);

  // report alignments in SPAM format
  void report_alignments_spam(std::string fname, KixRun* idx);

}; // END class DatasetAlignment 



#endif /* ALNBLOCK_H */
