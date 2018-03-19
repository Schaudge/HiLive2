#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "headers.h"

////////////////////////////////
////////// Exceptions //////////
////////////////////////////////

/** Errors that occur during I/O operations. */
class io_error : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

/** Errors that occur when opening a file. */
class file_open_error : public io_error {
public:
	using io_error::io_error;
};

/** Errors that occur if a file does not exist. */
class file_not_exist_error : public io_error {
public:
	using io_error::io_error;
};

/** Erorrs that occur because a file has the wrong format. */
class file_format_error : public io_error {
public:
	using io_error::io_error;
};

/**
 * Exception specialization for Unmodifiable data types.
 * @author Tobias Loka
 */
class immutable_error : public std::logic_error
{
public:
	using std::logic_error::logic_error;
};

/**
 * Error that is thrown when trying to modify a value of an immutable variable that was already set.
 * @author Tobias Loka
 */
class try_to_modify_immutable_error : public immutable_error
{
public:
	using immutable_error::immutable_error;
};

/**
 * Error that is thrown when trying to modify a value of an immutable variable that was already set.
 * @author Tobias Loka
 */
class try_to_get_unset_immutable_error : public immutable_error
{
public:
	using immutable_error::immutable_error;
};


////////////////////////////////////////
////////// Integer Data Types //////////
////////////////////////////////////////

/** Type for small counters. */
typedef uint16_t CountType;

/** Difference between positions in the read and matching position in the reference. */
typedef int16_t DiffType;

/** Data type for alignment scores. */
typedef int16_t ScoreType;

/** Type for the identifier of genomes (gid). */
typedef uint32_t GenomeIdType;

/** Type for positions in a genome. */
typedef uint32_t PositionType;


///////////////////////////////////////
////////// Complex Data Types /////////
///////////////////////////////////////

/** A pair of genome ID and position. */
struct GenomePosType {

	/** Genome ID */
	GenomeIdType gid;

  	/** Position in the genome */
	PositionType pos;

	GenomePosType() : gid(0), pos(0) {};
	GenomePosType(GenomeIdType g, PositionType p): gid(g), pos(p) {};
};
inline bool operator<(const GenomePosType l, const GenomePosType r) { if(l.gid==r.gid) return l.pos<r.pos; return l.gid<r.gid;}

/**
 * Template to store a map of mutexes.
 * Attention: Ensure that a locked mutex gets always unlocked (on destruction, if necessary).
 * If possible, use a combination of std::lock_guard and get_reference(T).
 */
template<typename K> class mutex_map {

private:

	/** The actual map of mutexes. */
	std::map<K, std::mutex> map;

	/** A mutex to control access to the map of mutexes. */
	std::mutex mut;

	/**
	 * Check if a key is already used in the map.
	 * @param k The key value.
	 * @return 0, if the key is not used in the current map. >0 if the key is used.
	 */
	typename std::map<K, std::mutex>::size_type count(K k) {
		return map.count(k);
	}

	/**
	 * Put a new mutex to the map if there the key is not used yet.
	 * @param k The key value.
	 * @return Reference to the mutex for the given key. This is a newly created mutex if the key was not used before,
	 * otherwise a previously existing mutex will be returned (which my be locked!).
	 */
	std::mutex& try_emplace(K k) {
		{
			std::lock_guard<std::mutex> lock(mut);
			if ( !count(k) )
				map.emplace(std::piecewise_construct, std::forward_as_tuple(k), std::forward_as_tuple());
			return map.at(k);
		}
	}

public:

	/**
	 * Unlock the mutex for a given key.
	 * @param k The key value.
	 */
	void unlock(K k) {
		if ( count(k) )
			map.at(k).unlock();
	}

	/**
	 * Lock the mutex for a given key.
	 * @param k The key value.
	 */
	void lock(K k) {
		try_emplace(k);
		map.at(k).lock();
	}

	/**
	 * Get a reference for a mutex for a specified key.
	 * @param k The key value.
	 * @return Reference to the mutex for the specified key. If no mutex exist for this key, create a new one.
	 */
	std::mutex& at(K k){
		return try_emplace(k);
	}

};

/**
 * A data type that increments an atomic arithmetic field for the time of it's existance.
 * This functionality can be used to block one slot of a specified capacity.
 */
template<
	typename T,
	typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type
> class atomic_increment_guard {

	/** Reference to the atomic field. */
	std::atomic<T>& val;

	/** Value of the atomic field during incrementation. */
	T blocked_value;

public:

	/**
	 * Constructor. Increments the value of the given reference and stores the value after incrementing.
	 * @param value Reference to the atomic, arithmetic value.
	 */
	atomic_increment_guard( std::atomic<T>& value ) : val(value), blocked_value(++val) { }

	/** Destructor. The reference value is decremented when destructing the guard. */
	~atomic_increment_guard() { --val; }

	/**
	 Get the value of the reference after being incremented.
	 @return The value of the reference after being incremented by this guard.
	 */
	T get_incremented_value(){
		return blocked_value;
	}
};

/**
 * Template to define data types that can only be set once.
 * @type T Data type of the unmodifiable object.
 * @author Tobias Loka
 */
template <
	typename T
> class Immutable {

private:

	/** The unmodifiable object. */
	T immutable_object;

	/** Flag to check if the object was already set once. */
	bool setFlag = false;

public:

	/** Constructor without setting the object (to only declare the object).*/
	Immutable(){	}

	/** Constructor with setting the object (to init the object).*/
	Immutable(T object) {
		set(object);
	}

	/** Automatic cast of the unmodifiable to the object type. */
	operator T() { return get(); }

	/**
	 * Set the unmodifiable object (will only work once!).
	 * @param object The object to be copied to this unmodifiable data type.
	 * @return true if setting was successful
	 * @throws immutable_error if the value was already set.
	 * @author Tobias Loka
	 */
	void set(T object) {
		if ( isSet() ) {
			throw try_to_modify_immutable_error("Tried to modify an immutable instance.");
		}

		immutable_object = object;
		setFlag = true;
	}

	/**
	 * Check if the object was already set.
	 * @return true if the object was already set.
	 * @author Tobias Loka
	 */
	bool isSet() const {
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
	T get(bool allow_unset = false ) const {
		if ( ! isSet() && ! allow_unset) {
			throw try_to_get_unset_immutable_error("Tried to access the data of an unset immutable.");
		}

		return immutable_object;
	}

};


///////////////////////////////////////
////////// Const definitions //////////
///////////////////////////////////////

/** Number of characters per line for help message. */
static const unsigned default_line_width = 80;

/** Size of a pair of genome ID and position (in bytes). */
const uint64_t GenomePos_size = sizeof(GenomeIdType) + sizeof(PositionType);

/** Supported nucleotides. */
const std::string seq_chars = "ACGTacgt";

/** Supported split characters for listing multitoken arguments. */
const std::string split_chars = " ,:;&|";

///** Define a mismatch as maximum value of DiffType. */
//const DiffType NO_MATCH = std::numeric_limits<DiffType>::max();
//
///** Define an insertion as maximum value of Difftype -1. */
//const DiffType INSERTION = std::numeric_limits<DiffType>::max()-1;
//
///** Define a deletion as maximum value of DiffType -2. */
//const DiffType DELETION = std::numeric_limits<DiffType>::max()-2;

/** Define the maximal number of positions to consider. */
const CountType MAX_NUM_POSITIONS = std::numeric_limits<CountType>::max();

/** Define a mask to only consider the latter two bits of a byte. */
static const uint8_t two_bit_mask = 3;

/** Define the value for undetermined barcodes as maximum of CountType. */
static const CountType UNDETERMINED = std::numeric_limits<CountType>::max();


//////////////////////////////////
////////// Enumerations //////////
//////////////////////////////////

/** The different alignment modes. */
enum AlignmentMode:char {
	ALL='A',
	ALLBEST='H',
	ANYBEST='B',
	BESTN='N',
	UNIQUE='U',
	UNKNOWN='Z'
};

/** The output formats of HiLive. */
enum OutputFormat:char {
	SAM		= 'S',	// SAM format
	BAM		= 'B',	// BAM format
	CRAM	= 'C',	// CRAM format //TODO: not supported yet.
};

/** SAM flags as defined in the specification. */
enum SAMFlag:uint16_t {
	MULT_SEG			= 1,	// Template has multiple segments
	ALL_SEG_MAPPED		= 2,	// Each segment properly aligned
	SEG_UNMAPPED		= 4,	// Segment unmapped
	NEXT_SEG_UNMAPPED	= 8,	// Next segment in the template unmapped
	SEQ_RC				= 16,	// SEQ being reverse complemented
	NEXT_SEQ_RC			= 32,	// SEQ of the next segment in the template being reverse complemented
	FIRST_SEG			= 64,	// the first segment in the template
	LAST_SEG			= 128,	// The last segment in the template
	SEC_ALIGNMENT		= 256,	// secondary alignment
	FILTER_NOT_PASSED	= 512,	// Not passing filteres, such as platform quality control
	PCR_OR_OPTICAL_DUPL	= 1024,	// PCR or optical duplicate
	SUPPL_ALIGNMENT		= 2048	// supplementary alignment
};

/** Alignment operations. */
enum Operations:uint8_t {
	MATCH = 0,
	NO_MATCH = 1,
	INSERTION = 2,
	DELETION = 3
};


/////////////////////////////////////////////
////////// Nucleotide representations ///////
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

#define comp(ch) ((toupper(ch)) == 'A' ? 'T' : \
				 (toupper(ch)) == 'C' ? 'G' : \
				 (toupper(ch)) == 'G' ? 'C' : \
				 (toupper(ch)) == 'T' ? 'A' : 'N');


////////////////////////////////////
////////// CIGAR elements //////////
////////////////////////////////////

/**
 * One (internal) CIGAR element.
 */
struct CigarElement {

	/** Length of the region. */
    CountType length;

    /** Operation of the region (Match, mismatch, insertion or deletion) */
    Operations operation;

    CigarElement (CountType l, Operations o): length(l), operation(o) {};
    CigarElement (): length(0), operation(NO_MATCH) {};
};

/**
 * Vector of CIGAR elements, representing the alignment information for one seed.
 */
typedef std::vector<CigarElement> CigarVector;


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
	bool isBarcode() const { return (mate==0);}
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


#endif /* DEFINITIONS_H */
