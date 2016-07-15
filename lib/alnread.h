#ifndef ALNREAD_H
#define ALNREAD_H

#include <chrono>

#include "headers.h"
#include "definitions.h"
#include "kindex.h"
#include "tools.h"


//-------------------------------------------------------------------//
//------  The Seed data structure  ----------------------------------//
//-------------------------------------------------------------------//


// a Seed stores the alignment of a read to a target genome
struct Seed {
  // internal sequence ID of taget genome
  GenomeIdType gid;
  // (estimated) start position of the read on the target
  PositionType start_pos;
  // number of matches
  CountType num_matches;

  // Information about matches/mismatches (similar to CIGAR). The last element is the current one
  CigarVector cigar_data;

  // get the size of the serialized object
  uint16_t serialize_size();

  // serialize the object
  std::vector<char> serialize();

  // deserialize (read) data from a char vector
  uint16_t deserialize(char* d);
};

typedef std::unique_ptr<Seed> USeed;

// compare function to sort Seed objects by position
bool seed_compare (Seed i,Seed j);
bool seed_compare_pos (const USeed & i, const USeed & j);
bool seed_compare_num_matches (const USeed & i, const USeed & j);

// std::vector of Seed pointers is much faster
typedef std::vector<USeed> SeedVec;
// a SeedVec Iterator
typedef SeedVec::iterator SeedVecIt;



//-------------------------------------------------------------------//
//------  The Read-Alignment class  ---------------------------------//
//-------------------------------------------------------------------//

class ReadAlignment {
 private:
  // read length
  CountType rlen;

  // Create new seeds from a list of kmer positions and add to current seeds
  void add_new_seeds(GenomePosListType& pos);


 public: // have everything public until the apropriate access functions are available

  // Flags for this read; 1 = read is valid (illumina flag)
  unsigned char flags;

  // the k-mer value observed in the last cycle
  HashIntoType last_kmer;

  // the last invalid cycle
  CountType last_invalid;

  // the current cycle
  CountType cycle;
  
  // a list of all found seeds
  SeedVec seeds;

  // Create a new read alignment given a certain read length
  //ReadAlignment(CountType rl): rlen(rl), last_kmer(0), last_invalid(0), cycle(0) {seeds.clear();};
  
  // assignment operator for deep copying
  //ReadAlignment& operator=(const ReadAlignment& other);
  
  // set the read_length
  void set_rlen(CountType r);
  
  // get the size of the serialized object
  uint64_t serialize_size();

  // serialize the object
  std::vector<char> serialize();

  // deserialize (read) data from a char vector
  uint64_t deserialize(char* d);

  // extend the alignment by one basecall using reference database index
  void extend_alignment(char bc, KixRun* index, AlignmentSettings* settings);

  // disable this alignment
  void disable();

  // generate the SAM flags for a seed
  uint32_t get_SAM_flags(uint32_t sd);

  // obtain start position of a seed according to SAM (leftmost) 
  PositionType get_SAM_start_pos(uint32_t sd);

  // calculate a quality score according to SAM
  uint16_t get_SAM_quality(uint32_t sd);

}; // END class ReadAlignment 






//------  Other helper functions  -----------------------------------//
//CountType num_matches(const std::vector<DiffType> &matches);
int16_t MAPQ(const SeedVec &sv);
std::string CIGAR(const Seed &seed);



#endif /* ALNREAD_H */
