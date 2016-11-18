#ifndef ALNREAD_H
#define ALNREAD_H

#include "headers.h"
#include "definitions.h"
#include "kindex.h"
#include "tools.h"
#include <seqan/basic.h>
#include <seqan/bam_io.h>


//-------------------------------------------------------------------//
//------  The Seed data structure  ----------------------------------//
//-------------------------------------------------------------------//


// a Seed stores the alignment of a read to a target genome
struct Seed {

  // internal sequence ID of taget genome
  GenomeIdType gid;

  // (estimated) start position of the read on the target
  PositionType start_pos;

  // number of matching bases
  CountType num_matches;

  // Information about matches/mismatches (similar to CIGAR). The last element is the current one
  CigarVector cigar_data;

  // return Seqans String of CigarElement
  seqan::String<seqan::CigarElement<> > returnSeqanCigarString();

  // get the size of the serialized object
  uint16_t serialize_size();

  // serialize the object
  std::vector<char> serialize();

  // deserialize (read) data from a char vector
  uint16_t deserialize(char* d);
};


typedef std::unique_ptr<Seed> USeed;
// compare function to sort Seed objects by position
bool seed_compare_pos (const USeed & i, const USeed & j);
// std::list of Seed pointers is much faster
typedef std::list<USeed> SeedVec;
// a SeedVec Iterator
typedef SeedVec::iterator SeedVecIt;



//-------------------------------------------------------------------//
//------  The Read-Alignment class  ---------------------------------//
//-------------------------------------------------------------------//


class ReadAlignment {

 private:
  // read length
  CountType rlen;

  // sequence of the read so far, saved as vector<uint8_t> so interpretation is not that trivial. Contains barcode
  CountType sequenceLen=0;
  std::vector<uint8_t> sequenceStoreVector;

  // Extend or create a placeholder seed for read with only trimmed matches
  void create_placeholder_seed();

  // convert a placeholder seed to a set of normal seeds
  void convertPlaceholder(GenomePosListType& pos);

  // Create new seeds from a list of kmer positions and add to current seeds
  void add_new_seeds(GenomePosListType& pos, std::vector<bool> & posWasUsedForExtension);

  // filter seeds based on filtering mode and q gram lemma. Also calls add_new_seeds.
  void filterAndCreateNewSeeds(GenomePosListType & pos, std::vector<bool> & posWasUsedForExtension);

  // updates cigar_data accordingly to a new matching kmer
  void addMatchingKmer(USeed & s, DiffType offset);

  // Extend an existing CIGAR string for a seed based on a new basecall. return false if last CIGAR element after extension is mismatch area (NO_MATCH), true otherwise.
  bool extendSeed(USeed & s, DiffType offset);

 public: // have everything public until the apropriate access functions are available

  // Flags for this read; 1 = read is valid (illumina flag)
  unsigned char flags = 1;

  // the last invalid cycle
  CountType last_invalid;

  // the current cycle
  CountType cycle;
  
  // a list of all found seeds
  SeedVec seeds;

  // max number of matches for this read
  CountType max_num_matches;

  // set the read_length
  void set_rlen(CountType r);
  
  // get the size of the serialized object
  uint64_t serialize_size();

  // serialize the object
  std::vector<char> serialize();

  // deserialize (read) data from a char vector
  uint64_t deserialize(char* d);

  // convert and return sequence of the read as string (without barcode)
  std::string getSequenceString();

  // convert and return sequence of the barcode
  std::string getBarcodeString();

  // append one nucleotide to sequenceStoreVector
  void appendNucleotideToSequenceStoreVector(char nuc);

  // extend the alignment by one basecall using reference database index
  void extend_alignment(char bc, KixRun* index);

  // disable this alignment
  void disable();

  // obtain start position of a seed according to SAM (leftmost) 
  PositionType get_SAM_start_pos(USeed & sd);

}; // END class ReadAlignment 



//-------------------------------------------------------------------//
//------  Other helper functions  -----------------------------------//
//-------------------------------------------------------------------//
int16_t MAPQ(const SeedVec &sv);

#endif /* ALNREAD_H */
