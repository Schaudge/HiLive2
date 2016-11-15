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

#endif /* DEFINITIONS_H */
