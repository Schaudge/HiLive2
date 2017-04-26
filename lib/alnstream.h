#ifndef ALNSTREAM_H
#define ALNSTREAM_H

#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "kindex.h"
#include "tools.h"
#include "alnread.h"
#include "illumina_parsers.h"
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

// Output alignment stream: write alignments to file one by one
class oAlnStream {
  // dataset information for the header
  uint16_t lane;
  uint16_t tile;
  uint16_t cycle;
  std::string root;
  CountType rlen;
  uint32_t num_reads;

  // number of reads written to file
  uint32_t num_written;

  // data buffer (don't write everything at once)
  std::vector<char> buffer;

  // size of the buffer
  uint64_t buf_size;

  // current position in the buffer
  uint64_t buf_pos;

  // output file format
  // 0: no compression
  // 1: zlib compression (level 1)
  // 11: lz4 compression (level 1)
  uint8_t format;

  // file handles
  FILE* ofile;
  gzFile ozfile;

  // write function for lz4 compression
  uint64_t lz4write(const char* buf, uint64_t size);

 public:
  // constructor initializes all member variables
  oAlnStream(uint16_t ln, uint16_t tl, uint16_t cl, std::string rt, CountType rl, uint32_t nr, uint64_t bs, uint8_t fmt);

  // open Alignment stream file and write header
  uint64_t open(std::string fname);

  // writes a read alignment to the output Alignment file. 
  // Buffering is handled internally
  uint64_t write_alignment(ReadAlignment * al);
  
  // checks if the correct number of alignments was written and closes the Alignment file
  bool close();
};



// Input alignment stream: loads read alignments from a file one by one
class iAlnStream {
  // dataset information for the header
  uint16_t lane;
  uint16_t tile;
  uint16_t cycle;
  std::string root;
  CountType rlen;
  uint32_t num_reads;

  // number of reads loaded from file
  uint32_t num_loaded;

  // data buffer (read blocks of data)
  std::vector<char> buffer;

  // size of the buffer
  uint64_t buf_size;

  // current position in the buffer
  uint64_t buf_pos;

  // output file format
  // 0: no compression
  // 1: zlib compression (level 1)
  // 11: lz4 compression (level 1)
  uint8_t format;

  // file pointer
  FILE* ifile;
  gzFile izfile;

  // read function for LZ4 compression. Reads one block of data to buffer
  uint64_t lz4read_block();

 public:
  // constructor initializes only block size and file format
  iAlnStream(uint64_t bs, uint8_t fmt);

  // open Alignment stream file and load header
  uint64_t open(std::string fname);

  // loads a read alignment from the input Alignment file. 
  // Buffering is handled internally
  ReadAlignment* get_alignment();
  
  // checks if the correct number of alignments was loaded and closes the Alignment file
  bool close();

  // get dataset information
  inline uint16_t get_lane() {return lane;};
  inline uint16_t get_tile() {return tile;};
  inline uint16_t get_cycle() {return cycle;};
  inline std::string get_root() {return root;};
  inline CountType get_rlen() {return rlen;};
  inline uint32_t get_num_reads() {return num_reads;};
  inline uint32_t get_num_loaded() {return num_loaded;};
};


//-------------------------------------------------------------------//
//------  The StreamedAlignment class  ------------------------------//
//-------------------------------------------------------------------//

class StreamedAlignment {
  
  // dataset information
  uint16_t lane;
  uint16_t tile;
  std::string root; // the BaseCalls directory
  CountType rlen;

  // fetch the next read from the input stream
  ReadAlignment get_next_read();
  
  // write an alignment to the output stream
  uint64_t write_alignment(ReadAlignment& ral);

  // get the path to the bcl file of a given cycle
  std::string get_bcl_file(uint16_t cycle, uint16_t read_number);

  // get the path to the alignment file. The alignment file is located in
  // <base>/L00<lane>/s_<lane>_<tile>.<cycle>.align
  // if base == "": base = root
  std::string get_alignment_file(uint16_t cycle, uint16_t mate, std::string base = "");

  // get the path to the filter file. The illumina filter information is located in
  // <root>/L00<lane>/s_<lane>_<tile>.filter
  std::string get_filter_file();

 public:
  StreamedAlignment(uint16_t ln, uint16_t tl, std::string rt, CountType rl): lane(ln), tile(tl), root(rt), rlen(rl) {};  

  StreamedAlignment& operator=(const StreamedAlignment& other);
  
  // create directories required to store the alignment files (only if not stored in root)
  void create_directories();

  // initialize empty alignment. Creates files for a virtual Cycle 0
  void init_alignment(uint16_t mate);
  
  // extend an existing alignment from cycle <cycle-1> to <cycle>
  uint64_t extend_alignment(uint16_t cycle, uint16_t read_no, uint16_t mate, KixRun* index);

  /**
   * Extend the barcode for all reads with the information of the current sequencing cycle.
   * @param bc_cycle The cycle of the barcode read.
   * @param read_cycle The last handled cycle for the respective mate (should always be 0 or the full length)
   * @param read_no The number of the sequence read for which the barcode will be extended (:= index in globalAlignmentSettings.seqs).
   * @param mate The read mate to extend the barcode.
   * @return
   * @author Tobias Loka
   */
  void extend_barcode(uint16_t bc_cycle, uint16_t read_cycle, uint16_t read_no, uint16_t mate);

}; /* END class StreamedAlignment */




//-------------------------------------------------------------------//
//------  Streamed SAM generation -----------------------------------//
//-------------------------------------------------------------------//

uint64_t alignments_to_sam(uint16_t ln, uint16_t tl, std::string rt, CountType rl, CountType mate, KixRun* index);


#endif /* ALNSTREAM_H */
