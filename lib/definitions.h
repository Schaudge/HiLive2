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
const std::string seq_chars = "ACGTacgt";

// largest number we're going to hash into. (8 bytes/64 bits/32 nt)
// probably 32 bit/16 nt are enough here
typedef uint64_t HashIntoType;

// construct a mask to truncate a binary representation of a k-mer to length K
const HashIntoType MASK = HashIntoType(pow(4,K_HiLive))-1;

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

// vector of ID:position pairs 
typedef std::vector<GenomePosType> GenomePosListType;

// iterator on GenomePosList
typedef GenomePosListType::iterator GenomePosListIt;

// the k-mer index array
const HashIntoType n_kmer = pow(4,K_HiLive);
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
typedef std::list<CigarElement> CigarVector;




// all user parameters are stored in the alignment settings
struct AlignmentSettings {
  // HARD CODED: kmer gap structure (this is not used anywhere)
  //std::string kmer_structure = "11111110111110111";
  std::string kmer_structure = "111111111111111";

  // HARD CODED: kmer gap positions (one-based)
  //std::vector<unsigned> kmer_gaps = {8, 14};
  std::vector<unsigned> kmer_gaps;

  // HARD CODED: kmer span (kmer weight is K_HiLive)
  //unsigned kmer_span = K_HiLive+2;
  unsigned kmer_span = K_HiLive;

  // PARAMETER: Base Call quality cutoff, treat BC with quality < bc_cutoff as miscall
  CountType min_qual;

  // PARAMETER: max. insert/deletion size
  DiffType window;

  // PARAMETER: minimum number of errors allowed in alignment
  CountType min_errors;

  // SWITCH: discard One-hit-wonders
  bool discard_ohw;

  // PARAMETER: first cycle to discard one-hit-wonders
  CountType start_ohw;

  // SWITCH: Best-Hit-Mode
  bool any_best_hit_mode;

  // SWITCH: Best-Hit-Mode
  bool all_best_hit_mode;

  // SWITCH: Best-N-Mode
  bool all_best_n_scores_mode;

  // PARAMETER: Best-N-Mode::N
  CountType best_n;

  // PARAMETER: temporary directory for the streamed alignment
  std::string temp_dir;

  // SWITCH: write sam/bam output or not
  bool write_bam=false;

  // SWITCH: Keep the old alignment files of previous cycles
  bool keep_aln_files;

  // PARAMETER: Memory block size for the input and output buffer in the streamed alignment
  uint64_t block_size;

  // PARAMETER: Compression format for alignment files
  uint8_t compression_format;

  // PARAMETER: list of lanes to process
  std::vector<uint16_t> lanes;
  
  // PARAMETER: list of tiles to process
  std::vector<uint16_t> tiles;

  // PARAMETER: root directory of hilive run
  std::string root;

  // PARAMETER: path to the index file
  std::string index_fname;

  // PARAMETER: read length of all reads (including barcodes)
  CountType rlen;

  // PARAMETER: length of the sequence of all reads (excluding barcodes)
  CountType seqlen;

  // PARAMETER: lengths of the respective sequencing reads
  std::vector<CountType> seqLengths;

  // PARAMETER: true, if barcode read, false if sequence read.
  std::vector<bool> isBarcode;

  // PARAMETER: vector containing all barcodes of the reads which should be outputted
  std::vector<std::string> barcodeVector;

  // PARAMETER: directory in which to create the output directory structure 
  std::string out_dir;

  // PARAMETER: number of threads to use
  CountType num_threads;
};




#endif /* DEFINITIONS_H */
