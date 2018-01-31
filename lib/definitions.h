#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "headers.h"
#include <seqan/index.h>


/////////////////////////////////////////////
////////// Sequences / Nucleotides //////////
/////////////////////////////////////////////

/**
 * Two-bit representation of a nucleotide.
 * @param ch Nucleotide as char
 * @return 2-bit representation
 */
#define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0LL : \
                         (toupper(ch)) == 'C' ? 1LL : \
                         (toupper(ch)) == 'G' ? 2LL : 3LL)

/**
 * Complementary two-bit representation of a nucleotide.
 * @param ch Nucleotide as char
 * @return Complementary 2-bit representation
 */
#define twobit_comp(ch) ((toupper(ch)) == 'A' ? 3LL : \
                         (toupper(ch)) == 'C' ? 2LL : \
                         (toupper(ch)) == 'G' ? 1LL : 0LL)

/**
 * Nucleotide of a 2-bit representation.
 * @param n 2-bit representation of a nucleotide.
 * @return Nucleotide as char
 */
#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'C' : \
                           (n) == 2 ? 'G' : 'T')

#define revfourbit_repr(n) ((n) <= 3 ? 'N' : \
							revtwobit_repr(n & 3))

#define comp(ch) ((toupper(ch)) == 'A' ? 'T' : \
				 (toupper(ch)) == 'C' ? 'G' : \
				 (toupper(ch)) == 'G' ? 'C' : \
				 (toupper(ch)) == 'T' ? 'A' : 'N');

/**
 * Supported nucleotides.
 */
const std::string seq_chars = "ACGTacgt";


////////////////////////////////////////
////////// Genome Identifiers //////////
////////////////////////////////////////

/**
 * Type for the identifier of genomes (gid).
 */
typedef uint32_t GenomeIdType;

/**
 * Constant variable to tag a k-mer as "trimmed".
 */
const GenomeIdType TRIMMED = std::numeric_limits<GenomeIdType>::max();

/**
 * A list of Genome Ids
 */
typedef std::vector<GenomeIdType> GenomeIdListType;


//////////////////////////////////////
////////// Genome Positions //////////
//////////////////////////////////////

/**
 * Type for positions in a genome.
 */
typedef uint32_t PositionType;

// GenomeID / Position pair
typedef std::pair<GenomeIdType, PositionType> PositionPairType;

// List of Genome ID / Position pairs
typedef std::list<PositionPairType> PositionPairListType;

/**
 * A pair of genome ID and position.
 */
struct GenomePosType {

  GenomeIdType gid;
  PositionType pos;

  GenomePosType()=default;
  GenomePosType(GenomeIdType g, PositionType p): gid(g), pos(p) {};
};

/**
 * Size of a pair of genome ID and position (in bytes)
 */
const uint64_t GenomePos_size = sizeof(GenomeIdType) + sizeof(PositionType);

/**
 * A vector of GenomePosTypes.
 */
typedef std::vector<GenomePosType> GenomePosListType;

/**
 * Iterator on GenomePosListType.
 */
typedef GenomePosListType::iterator GenomePosListIt;


/////////////////////////////////
////////// K-mer index //////////
/////////////////////////////////

/**
 * Type to hash k-mers into.
 * This type also limits the k-mer weight (currently to 32).
 */
typedef uint64_t HashIntoType;

/**
 * K-mer index type.
 */
typedef std::vector<GenomePosListType> KmerIndexType;

/**
 * A lightweight type for storing the index.
 */
typedef std::vector<char*> KixRunDB;


//////////////////////////////
////////// FM index //////////
//////////////////////////////

/**
 * FM-Index Config.
 */
typedef seqan::FastFMIndexConfig<void, uint64_t,2 ,1> FMIConfig;

/**
 * FM-Index data type.
 */
typedef seqan::Index<seqan::StringSet<seqan::DnaString>, seqan::FMIndex<void, FMIConfig> > FMIndex;

/**
 * Iterator to go through the FM index structure.
 */
typedef seqan::Iterator<FMIndex, seqan::TopDown<> >::Type FMTopDownIterator;

/**
 * Data type to find a store a node in the FM index.
 */
typedef seqan::Iter<FMIndex,seqan::VSTree<seqan::TopDown<seqan::Preorder>>>::TVertexDesc FMVertexDescriptor;

inline bool operator==(const FMVertexDescriptor l, FMVertexDescriptor r) {
	return l.range == r.range;
//	return std::memcmp(&(l),&(r),sizeof(l)) == 0;
}

inline bool operator<(const FMVertexDescriptor l, FMVertexDescriptor r) {
	if ( l.range.i1 == r.range.i1)
		return l.range.i2 < r.range.i2;
	return l.range.i1 < r.range.i1;
//	return std::memcmp(&(l),&(r),sizeof(l)) < 0;
}

inline bool operator>(const FMVertexDescriptor l, FMVertexDescriptor r) {
	if ( l.range.i1 == r.range.i1)
		return l.range.i2 > r.range.i2;
	return l.range.i1 > r.range.i1;
//	return std::memcmp(&(l),&(r),sizeof(l)) > 0;
}

inline bool operator!=(const FMVertexDescriptor l, FMVertexDescriptor r) {
	return !(l==r);
//	return std::memcmp(&(l),&(r),sizeof(l)) != 0;
}


////////////////////////////////////////
////////// Integer data types //////////
////////////////////////////////////////

/**
 * Type for small counters.
 */
typedef uint16_t CountType;

/**
* Difference between k-mer position in the read and matching position in the reference.
*/
typedef int16_t DiffType;

/** Data type for alignment scores. */
typedef int16_t ScoreType;

////////////////////////////////////////
////////// Offset definitions //////////
////////////////////////////////////////

/**
 * Define a mismatch as maximum value of DiffType.
 */
const DiffType NO_MATCH = std::numeric_limits<DiffType>::max();

/**
 *  Define an insertion as maximum value of Difftype -2.
 */
const DiffType INSERTION = std::numeric_limits<DiffType>::max()-1;

/**
 * Define a deletion as maximum value of DiffType -3.
 */
const DiffType DELETION = std::numeric_limits<DiffType>::max()-2;


////////////////////////////////////
////////// CIGAR elements //////////
////////////////////////////////////

/**
 * One (internal) CIGAR element.
 */
struct CigarElement {

	/** Length of the region. */
    CountType length;

    /** Offset of the region to the original start pos (created by InDels). */
    DiffType offset;

    CigarElement (CountType l, DiffType o): length(l), offset(o) {};
    CigarElement (): length(0), offset(NO_MATCH) {};
};

/**
 * Vector of CIGAR elements, representing the alignment information for one seed.
 */
typedef std::list<CigarElement> CigarVector;


///////////////////////////////////////
////////// Sequence Elements //////////
///////////////////////////////////////

/**
 * Information about the sequences.
 * One element can be a read or a barcode.
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

/**
 * Check if two Sequence elements are equal.
 */
inline bool operator==(const SequenceElement l, const SequenceElement r) {return (l.length==r.length) && (l.mate==r.mate) && (l.id==r.id);}

/**
 * Checks if two sequence elements are not equal.
 */
inline bool operator!=(const SequenceElement l, const SequenceElement r) {return !(l==r);}

/**
 * An undefined sequence element (NULL element).
 */
const SequenceElement NULLSEQ = SequenceElement();


////////////////////////////////////////////
////////// Unmodifiable variables //////////
////////////////////////////////////////////

/**
 * Exception specialization for Unmodifiable data types.
 * @author Tobias Loka
 */
class unmodifiable_error : public std::logic_error
{
public:
	using std::logic_error::logic_error;
};

/**
 * Template to define data types that can only be set once.
 * @type T Data type of the unmodifiable object.
 * @author Tobias Loka
 */
template <typename T>
class Unmodifiable {

private:

	/** The unmodifiable object. */
	T unmodifiable_object;

	/** Flag to check if the object was already set once. */
	bool setFlag = false;

public:

	/** Constructor without setting the object (to only declare the object).*/
	Unmodifiable(){	}

	/** Constructor with setting the object (to init the object).*/
	Unmodifiable(T object) {
		unmodifiable_object = object;
	}

	/** Automatic cast to of the unmodifiable to the object type. */
	operator T() { return unmodifiable_object; }

	/**
	 * Set the unmodifiable object (will only work once!).
	 * @param object The object to be copied to this unmodifiable data type.
	 * @return true if setting was successful
	 * @author Tobias Loka
	 */
	void set(T object) {
		if ( isSet() ) {
			throw unmodifiable_error("Tried to modify unmodifiable object");
		}

		unmodifiable_object = object;
		setFlag = true;
	}

	/**
	 * Check if the object was already set.
	 * @return true if the object was already set.
	 * @author Tobias Loka
	 */
	bool isSet() {
		return setFlag;
	}

	/**
	 * Return a copy of the unmodifiable object.
	 * @param allow_unset if false, an exception is thrown when the object was not set before. Should only be true for
	 * objects that require access to certain properties before their initialization (e.g. to check a container's size
	 * without knowing if the container was already set).
	 * @return (copy/value of) the unmodifiable object
	 * @author Tobias Loka
	 */
	T get(bool allow_unset = false ) {
		if ( ! isSet() && ! allow_unset) {
			throw unmodifiable_error("Tried to access uninitialized object");
		}

		return unmodifiable_object;
	}

};

template<typename T>
bool set_unmodifiable(Unmodifiable<T> & unmodifiable, T value, std::string variable_name) {
	  try {
		  unmodifiable.set(value);
	  }
	  catch (unmodifiable_error& e) {
//		std::cerr << e.what() << " (" << variable_name << ")." << std::endl;
		  variable_name.length(); // TODO: just to remove warnings. Remove variable_name string when finished.
		 return false ;
	  }
	  return true;
}

template<typename T>
T get_unmodifiable(Unmodifiable<T> unmodifiable, std::string variable_name, bool allow_unset = false) {
	  try {
		  return unmodifiable.get(allow_unset);
	  }
	  catch (unmodifiable_error& e) {
		  std::cerr << e.what() << " (" << variable_name << ")." << std::endl;
		  return T();
	  }
}

////////////////////////////////
////////// Exceptions //////////
////////////////////////////////
class io_error : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

class file_open_error : public io_error {
public:
	using io_error::io_error;
};

class file_not_exist_error : public io_error {
public:
	using io_error::io_error;
};

class file_format_error : public io_error {
public:
	using io_error::io_error;
};

///////////////////////////////////////
////////// Other definitions //////////
///////////////////////////////////////

const CountType MAX_NUM_POSITIONS = std::numeric_limits<CountType>::max();

/**
 * A list of strings
 */
typedef std::vector<std::string> StringListType;

/**
 * The different alignment modes.
 * @author Tobias Loka
 */
enum AlignmentMode:char {
	ALL='A',
	ALLBEST='H',
	ANYBEST='B',
	BESTN='N',
	UNIQUE='U',
	UNKNOWN='Z'
};

/**
 * SAM flags as defined in the specification.
 */
enum SAMFlag:uint16_t {
	MULT_SEG			= 1,	// Template has multiple segments
	ALL_SEG_MAPPED		= 2,	// Each segment properly aligned
	SEG_UNMAPPED		= 4,	// Segment unmapped
	NEXT_SEG_UNMAPPED	= 8,	// Next segment in the template unmapped
	SEQ_RC				= 16,	// SEQ being reverse complemented
	NEXT_SEQ_RC			= 32,	// SEQ of the next segment in the template being reverse complemented
	FIRST_SEG			= 64,	// the first segment inn the template
	LAST_SEG			= 128,	// The last segment in the template
	SEC_ALIGNMENT		= 256,	// secondary alignment
	FILTER_NOT_PASSED	= 512,	// Not passing filteres, such as platform quality control
	PCR_OR_OPTICAL_DUPL	= 1024,	// PCR or optical duplicate
	SUPPL_ALIGNMENT		= 2048	// supplementary alignment
};

/**
 * Template to store a map of mutexes.
 * Ensure that a locked mutex gets always unlocked (on destruction, if necessary). If possible, use a combination of std::lock_guard and get_reference(T).
 */
template<typename K> class mutex_map {

private:
	std::map<K, std::mutex> map;
	std::mutex mut;

	typename std::map<K, std::mutex>::size_type count(K k) {
		return map.count(k);
	}

	std::mutex& try_emplace(K k) {
		{
			std::lock_guard<std::mutex> lock(mut);
			if ( !count(k) )
				map.emplace(std::piecewise_construct, std::forward_as_tuple(k), std::forward_as_tuple());
			return map.at(k);
		}
	}

public:

	void unlock(K k) {
		if ( count(k) )
			map.at(k).unlock();
	}

	void lock(K k) {
		try_emplace(k);
		map.at(k).lock();
	}

	std::mutex& at(K k){

		return try_emplace(k);
	}

};

/**
 * A data type that increments an arithmetic field for the time of it's existance.
 * This functionality can be used to block one slot of a certain capacity.
 */
template<
	typename T,
	typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type
> class block_guard {
	T& val;
	T blocked_value;
public:
	block_guard( T& value ) : val(value), blocked_value(++val){ }
	~block_guard() { --val; }
	T get_blocked_value(){ return blocked_value; }
};

/**
 * A data type that increments an atomic arithmetic field for the time of it's existance.
 * This functionality can be used to block one slot of a certain capacity.
 */
template<
	typename T,
	typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type
> class atomic_block_guard {
	std::atomic<T>& val;
	T blocked_value;
public:
	atomic_block_guard( std::atomic<T>& value ) : val(value), blocked_value(++val) { }
	~atomic_block_guard() { --val; }
	T get_blocked_value(){ return blocked_value; }
};

#endif /* DEFINITIONS_H */
