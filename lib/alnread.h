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
  seqan::String<seqan::CigarElement<> > returnSeqanCigarString(unsigned* nm_i);

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
  CountType total_cycles;

  // sequence of the read so far, saved as vector<uint8_t> so interpretation is not that trivial.
  CountType sequenceLen=0;
  std::vector<uint8_t> sequenceStoreVector;

  // sequence of the barcode so far, saved as vector<uint8_t> so interpretation is not that trivial
  CountType barcodeLen=0;
  std::vector<uint8_t> barcodeStoreVector;

  // Extend or create a placeholder seed for read with only trimmed matches
  void create_placeholder_seed();

  // Create new seeds from a list of kmer positions and add to current seeds
  void add_new_seeds(GenomePosListType& pos, std::vector<bool> & posWasUsedForExtension);

  /**
   * This function is the modified pigeonhole principle holding for both spaced and unspaced kmers.
   * It computes the minimum number of errors in an error region of a given CIGAR vector.
   * An error region is a region that is surrounded by MATCH elements of length >= ( kmer_span - 1 ).
   * The error region cannot contain MATCH elements of length >= ( kmer_span - 1 ).
   *
   * @param region_length Sum of all (!) elements within the error region, including involved MATCH elements.
   * @param border Number of included borders of the CIGAR vector (begin and/or end). Must be in [0,2].
   * @param Absolute number (positive) of the offset change during a region
   * @return Minimum number of errors that caused a region of the given length.
   * @author Tobias Loka, Jakob Schulze
   */
  CountType minErrors_in_region(CountType region_length, CountType border, CountType offset_change=0 );

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
  void set_total_cycles(CountType c);
  
  // get the size of the serialized object
  uint64_t serialize_size();

  // serialize the object
  std::vector<char> serialize();

  // deserialize (read) data from a char vector
  uint64_t deserialize(char* d);

  // convert and return sequence of the read as string (without barcode)
  std::string getSequenceString();

  /**
   * Convert and return sequence of the barcode. Multiple barcodes are concatenated (without delimiter).
   * @return The Barcode as string
   * @author Tobias Loka
   */
  std::string getBarcodeString();

  /**
     * Check whether the barcode of this read fulfills the criteria of at least one user-defined barcode.
     * The nucleotides are only compared pairwise, not allowing for Indels.
     * @return The index of the matching barcode in globalAlignmentSettings.multiBarcodeVector. NO_MATCH, if none.
     * Also return NO_MATCH, if demultiplexing is not activated.
     * @author 	Tobias Loka
     */
  CountType getBarcodeIndex() ;


  /**
   * Append one nucleotide to sequenceStoreVector
   * @param nucl The nucleotide. Must be 2-bit-formatted.
   * @param appendToBarcode If true, the nucleotide is appended to the barcode instead of the read sequence (default: false).
   * @return
   * @author Jakob Schulze
   */
  void appendNucleotideToSequenceStoreVector(char nuc, bool appendToBarcode=false);

  // extend the alignment by one basecall using reference database index
  void extend_alignment(char bc, KixRun* index);

  // disable this alignment
  void disable();

  // obtain start position of a seed according to SAM (leftmost) 
  PositionType get_SAM_start_pos(USeed & sd);

  /**
   * Compute the minimum number of errors for a seed by using the modified pigeonhole principle implemented in ReadAlignment::minErrors_in_region.
   *
   * @param s The seed.
   * @return The minimum number of errors for the given seed.
   * @author Tobias Loka, Jakob Schulze
   */
  CountType min_errors(USeed & s);

}; // END class ReadAlignment 



//-------------------------------------------------------------------//
//------  Other helper functions  -----------------------------------//
//-------------------------------------------------------------------//
int16_t MAPQ(const SeedVec &sv);

#endif /* ALNREAD_H */
