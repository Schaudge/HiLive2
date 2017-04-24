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

/**
 * Struct containing information about the user-defined sequencing reads. This can be barcode or sequence fragments.
 * @author Tobias Loka
 */
struct SequenceElement {

	/** The id of the read. Equals the position in the argument list and in the AlignmentSettings::seqs vector (0-based). */
	CountType id;
	/** The mate number. 0 for barcodes, increasing for sequence reads in the given order (1-based). */
	CountType mate;
	/** The length of the respective read. */
	CountType length;

	/**
	 * Constructor of a SequenceElement NULL object.
	 * @author Tobias Loka
	 */
	SequenceElement () : id(0), mate(0), length(0) {};

	/**
	 * Constructor of a valid SequenceElement object.
	 * @param id The id of the read.
	 * @param m The mate number of the read (0 for barcodes, incrementing for sequence reads)
	 * @param l The length of the read
	 * @author Tobias Loka
	 */
	SequenceElement (CountType id, CountType m, CountType l): id(id), mate(m), length(l) {};

	/**
	 * Check whether the SequenceElement object is a barcode or not.
	 * @return true, if SequenceElement is a barcode. False if not.
	 * @author Tobias Loka
	 */
	bool isBarcode() { return (mate==0);}
};

/** Checks for equality of two SequenceElement objects*/
inline bool operator==(const SequenceElement l, const SequenceElement r) {return (l.length==r.length) && (l.mate==r.mate) && (l.id==r.id);}
/** Checks for inequality of two SequenceElement objects*/
inline bool operator!=(const SequenceElement l, const SequenceElement r) {return !(l==r);}
/** Defines a null-SequenceElement. */
const SequenceElement NULLSEQ = SequenceElement();

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

  // PARAMETER: total number of cycles
  CountType cycles;

  /**
   * Contains the read information of the sequencing machine (as SequenceElement objects). Includes sequence reads and barcodes.
   * Arbitrary numbers and orders of reads are supported. The summed length of all elements must equal the number of sequencing cycles.
   * @author Tobias Loka
   */
  std::vector<SequenceElement> seqs;

  // Number of mates (information taken from the seqLengths parameter), (Hint: corresponding indeces are 1-based)
  uint16_t mates;

  //PARAMETER: Stores the barcodes defined by the user. The inner vector contains the single fragments of multi-barcodes.
  std::vector<std::vector<std::string>> barcodeVector;

  // PARAMETER: number of allowed errors for the single barcodes
  std::vector<uint16_t> barcode_errors;

  // SWITCH: if true, keep all barcodes (disables barcode filtering).
  bool keep_all_barcodes;

  // PARAMETER: directory in which to create the output directory structure 
  boost::filesystem::path out_dir;

  // PARAMETER: number of threads to use
  CountType num_threads;

  /**
   * Get a SequenceElement object from the seqs vector by using the id
   * @param id The id of the SequenceElement.
   * @return The respective SequenceElement object for the given id.
   * @author Tobias Loka
   */
  SequenceElement getSeqById(CountType id) {return seqs[id];}

  /**
   * Get a SequenceElement object from the seqs vector by using the mate number
   * @param id The mate number of the SequenceElement.
   * @return The respective SequenceElement object for the given mate number. NULLSEQ if mate==0 (barcodes).
   * @author Tobias Loka
   */
  SequenceElement getSeqByMate(CountType mate) {
	  if ( mate == 0 ) return NULLSEQ;
	  for (uint16_t i = 0; i != seqs.size(); i++) {
		  if(seqs[i].mate == mate) return seqs[i];
	  }
	  return NULLSEQ;
  }
};




#endif /* DEFINITIONS_H */
