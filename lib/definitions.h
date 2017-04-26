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

// the k-mer index type
typedef std::vector<GenomePosListType> KmerIndexType;

// a lightweight type for storing the index
typedef std::vector<char*> KixRunDB;

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


#endif /* DEFINITIONS_H */
