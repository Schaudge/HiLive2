#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "headers.h"



// bit representation of A/C/G/T.
#define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0LL : \
                         (toupper(ch)) == 'C' ? 1LL : \
                         (toupper(ch)) == 'G' ? 2LL : 3LL)

// complement bit representation of A/C/G/T.
#define twobit_comp(ch) ((toupper(ch)) == 'A' ? 3LL : \
                         (toupper(ch)) == 'C' ? 2LL : \
                         (toupper(ch)) == 'G' ? 1LL : 0LL)

// bit representation to character
#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'C' : \
                           (n) == 2 ? 'G' : 'T')




// Allowed characters in sequence
const std::string seq_chars = "ACGT";

// largest number we're going to hash into. (8 bytes/64 bits/32 nt)
// probably 32 bit/16 nt are enough here
typedef uint64_t HashIntoType;

// construct a mask to truncate a binary representation of a k-mer to length K
const HashIntoType MASK = HashIntoType(pow(4,K))-1;


// identifiers for genome sequences
typedef uint32_t GenomeIdType;
const GenomeIdType TRIMMED = std::numeric_limits<GenomeIdType>::max();

// list of genome identifiers
typedef std::vector<GenomeIdType> GenomeIdListType;

// list of strings
typedef std::vector<std::string> StringListType;


// position in a genome
typedef int32_t PositionType;

// pair of genome ID and position
struct GenomePosType {
  GenomeIdType gid;
  PositionType pos;

  GenomePosType()=default;
  GenomePosType(GenomeIdType g, PositionType p): gid(g), pos(p) {};
};
//typedef std::tuple<GenomeIdType,PositionType> GenomePosType;

// size of a pair of genome ID and position
const uint64_t GenomePos_size = sizeof(GenomeIdType) + sizeof(PositionType);

// compare function to sort GenomePosType objects by position
bool gp_compare (GenomePosType i,GenomePosType j);

// vector of ID:position pairs 
typedef std::vector<GenomePosType> GenomePosListType;

// iterator on GenomePosList
typedef GenomePosListType::iterator GenomePosListIt;


// the k-mer index array
const HashIntoType n_kmer = pow(4,K);
typedef std::array<GenomePosListType,n_kmer> KmerIndexType;


// small counters
typedef uint16_t CountType;


// difference between k-mer position in the read and matching position in the reference
typedef int16_t DiffType;

// define a mismatch as max(DiffType)
const DiffType NO_MATCH = std::numeric_limits<DiffType>::max();

// define a trimmed position as max(DiffType)-1
const DiffType TRIMMED_MATCH = std::numeric_limits<DiffType>::max()-1;


// one element in Cigar vector containing match/mismatch information about consecutive k-mers
struct CigarElement {
    CountType length;
    DiffType offset;
    CigarElement (CountType l, DiffType o): length(l), offset(o) {};
    CigarElement (): length(0), offset(NO_MATCH) {};
};

// CigarVector containing CIGAR string like information about the alignments
typedef std::vector<CigarElement> CigarVector;


// all user parameters are stored in the alignment settings
struct AlignmentSettings {
  // PARAMETER: Base Call quality cutoff, treat BC with quality < bc_cutoff as miscall
  uint8_t min_qual = 1;

  // PARAMETER: max. insert/deletion size
  DiffType window = 50;

  // PARAMETER: minimum number of errors allowed in alignment
  CountType min_errors = 6;

  // SWITCH: discard One-hit-wonders
  bool discard_ohw = true;

  // PARAMETER: first cycle to discard one-hit-wonders
  CountType start_ohw = K+5;

  // SWITCH: Best-Hit-Mode
  bool best_hit_mode = true;

  // SWITCH: Best-N-Mode
  bool best_n_mode = false;

  // PARAMETER: Best-N-Mode::N
  CountType best_n = 1;

  // SWITCH: sort positions found in index. Saves you a few seconds when turned off, but messes everything up when index is not sorted.
  bool sort_positions = false;

  // PARAMETER: temporary directory for the streamed alignment
  std::string temp_dir = "";

  // PARAMETER: directory for SAM file output
  std::string sam_dir = "";

  // SWITCH: Keep the old alignment files of previous cycles
  bool keep_aln_files = true;

  // PARAMETER: Memory block size for the input and output buffer in the streamed alignment
  uint64_t block_size = 64*1024*1024; /* 64 MB */

  // PARAMETER: Compression format for alignment files
  uint8_t compression_format = 2;

  // initialize with default values
  AlignmentSettings() : min_qual(1), 
                        window(5),
                        min_errors(2),
                        discard_ohw(true), 
                        start_ohw(K+5), 
                        best_hit_mode(true), 
                        best_n_mode(false),
                        best_n(1),
                        sort_positions(false),
                        temp_dir(""),
                        sam_dir(""),
                        keep_aln_files(true),
                        block_size(64*1024*1024),
                        compression_format(2) {};

};




#endif /* DEFINITIONS_H */
