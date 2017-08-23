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

/**
 * Data structure to save alignment information for a path in the FM index
 * @author Martin Lindner, Tobias Loka
 */
struct Seed {

	/**
	 * Vertex iterator for the FM index. Contains the position information of the alignment
	 */
	FMVertexDescriptor vDesc;

	/**
	 * Minimal number of errors (minimal in terms of the front softclip, number of errors after seeding is exact)
	 */
	CountType num_errors;

	/**
	 * Information about matches/mismatches (similar to CIGAR). The last element is the current one
	 */
	CigarVector cigar_data;

	/**
	 * Get the CIGAR string in SeqAn format
	 * @return CIGAR string in SeqAn format
	 */
	seqan::String<seqan::CigarElement<> > returnSeqanCigarString();

	/**
	 * Determine size of the serialized seed
	 * @return size in bytes
	 */
	uint16_t serialize_size();

	/**
	 * Serialize the seed
	 * @return serialized seed
	 */
	std::vector<char> serialize();

	/**
	 * Deserialize seed data from a char vector
	 * @param Char vector containing the seed information in binary format
	 * @return Size of loaded data in bytes
	 */
	uint16_t deserialize(char* d);

	void cout();
};

typedef std::shared_ptr<Seed> USeed;

/**
 * Shared pointer to seeds.
 */
typedef std::vector<USeed> SeedVec;

/**
 * Iterator for a list of shared pointers to seeds.
 */
typedef SeedVec::iterator SeedVecIt;

/**
 * Compare two pointers by comparing the respective target objects.
 */

template <typename T> bool PComp(const T & a, const T & b)
{
   return *a < *b;
}

/**
 * Comparator for seeds to sort them by their number of errors.
 * @author Tobias Loka
 */
inline bool seed_comparison_by_error(const USeed a, const USeed b) {

	// TODO: think about whether it makes sense to count complete softclip as mismatch or not (currently not)

	// if equal number of matches
	if (a->num_errors == b->num_errors) {

		// prefer no front softclip
		if ( a->cigar_data.front().offset != b->cigar_data.front().offset ) {
			return ( a->cigar_data.front().offset != NO_MATCH );
		}

		// prefer shorter front softclip
		if ( a->cigar_data.front().offset == NO_MATCH && b->cigar_data.front().offset == NO_MATCH) {
			if ( a->cigar_data.front().length != b->cigar_data.front().length) {
				return ( a->cigar_data.front().length < b->cigar_data.front().length );
			}
		}

		// prefer less cigar elements (if equal, don't care)
		return a->cigar_data.size() < b->cigar_data.size();
	}

	return a->num_errors < b->num_errors;
}


/**
 * Define '<'-operator for seeds.
 * Used for sorting by the range in the FM index.
 * @author Tobias Loka
 */
inline bool operator <(const Seed l, const Seed r) {

	// Lower range wins
	if ( l.vDesc < r.vDesc )
		return true;
	else if ( l.vDesc != r.vDesc )
		return false;

	// If range equal:
	else {

		auto lCig = l.cigar_data.begin();
		auto rCig = r.cigar_data.begin();

		// Different offsets -> MATCH wins.
		if ( lCig->offset != rCig->offset ) {
			return ( lCig->offset != NO_MATCH );
		}

		// Equal offsets
		else {

			// both match or equal length -> less errors win
			if ( lCig->offset != NO_MATCH || lCig->length == rCig->length )
				return l.num_errors < r.num_errors;

			// both NO_MATCH and different length -> lower length wins
			else
				return lCig->length < rCig->length;
		}

	}
	return false;
}


//-------------------------------------------------------------------//
//------  The Read-Alignment class  ---------------------------------//
//-------------------------------------------------------------------//

/**
 * Class for the organization of the complete alignment of a single read.
 */
class ReadAlignment {

 private:

	/** Number of cycles for the current alignment */
	CountType total_cycles;

	/** Length of the read sequence */
	CountType sequenceLen=0;

	/** Sequence of the read */
	std::vector<uint8_t> sequenceStoreVector;

	/** Length of the barcode */
	CountType barcodeLen=0;

	/** Barcode Sequence of the read */
	std::vector<uint8_t> barcodeStoreVector;

	/**
	 * Extend a single seed. Creates new seeds for all possible extensions that occur in the reference genome (Match, InDel, SNP)
	 * @param base The next nucleotide
	 * @param s The seed to be extended
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void extendSeed(char base, USeed s, CountType allowedErrors, KixRun* index, SeedVec & newSeeds);

	/**
	 * Extend a seed by alignment matches (Match or SNP)
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed to be extended
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void getMatchSeeds(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds);

	/**
	 * Extend a seed by an insertion
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed to be extended
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
     * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void getInsertionSeeds(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds);


	/**
	 * Extend a seed by deletions
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed to be extended
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void getDeletionSeeds(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds);

	/**
	 * Add deletions to the alignment up to the permitted number of errors
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed of the previous iteration
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void recursive_goDown(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds);

	/** Create new seeds
	 * @param index The FM index
	 * @param settings The alignment settings
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void createSeeds(KixRun* index, SeedVec & newSeeds);



 public:

	/** Flag of the illumina read (1=valid) */
	unsigned char flags = 1;

	/** The current alignment cycle */
	CountType cycle;
  
	/** List of all seeds for the respective read */
	SeedVec seeds;
  
	/**
	 * Get the size of the serialized alignment
	 * @return Size of the serialized alignment in bytes
	 */
	uint64_t serialize_size();

	/**
	 * Serialize the alignment
	 * @return Serialized alignment.
	 */
	std::vector<char> serialize();

	/**
	 * Deserialize alignment data from a char vector
	 * @param d The alignment data
	 * @return Number of handled bytes
	 */
	uint64_t deserialize(char* d);

	/**
	 * Get read sequence as string (without barcode).
	 * @return The read sequence
	 */
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
	void appendNucleotideToSequenceStoreVector(char bc, bool appendToBarcode=false);

	/**
	 * Extend alignment by 1 base.
	 * @param bc The next base (A,C,G or T)
	 * @param index The FM index
	 * @param settings The alignment settings
	 * @author Tobias Loka
	 */
	void extend_alignment(char bc, KixRun* index, bool testPrint=false);

	/**
	 * Disable the alignment.
	 */
	void disable();

	/**
	 * Determine the maximal number of errors that are permitted for a seed
	 * @param s The seed
	 * @return Maximal number of permitted errors
	 * @author Tobias Loka
	 */
	CountType getMaxNumErrors(USeed s);

	/**
	 * Obtain the start position of the alignment with SAM specifications (most left position)
	 * @param index The FM index
	 * @param p The position as stored in the index
	 * @param sd The seed containing the alignment information (CIGAR)
	 * @return The position in SAM specification
	 * @author Tobias Loka
	 */
	PositionType get_SAM_start_pos(KixRun* index, PositionPairType p, USeed & sd);


	/**
	 * Get all positions of a seed
	 * @param index The FM index
	 * @param sd The respective seed
	 * @param position_list Reference to the list of positions (all resulting positions are added to this list)
	 * @author Tobias Loka
	 */
	void getPositions(KixRun* index, USeed sd, PositionPairListType & position_list);


	/**
	 * Get all seeds for the respective read sorted by errors
	 * @return Vector containing sorted seeds
	 * @author Tobias Loka
	 */
	void getSeeds_errorsorted(SeedVec & seeds_sorted);


	/**
	 * Constructor for an object of type ReadAlignment
	 * @param tot_cyc Number of total alignment cycles
	 * @param cyc Current alignment cycle
	 * @author Tobias Loka
	 */
	ReadAlignment(CountType tot_cyc, CountType cyc){total_cycles=tot_cyc; cycle=cyc;};

}; // END class ReadAlignment 



//-------------------------------------------------------------------//
//------  Other helper functions  -----------------------------------//
//-------------------------------------------------------------------//
int16_t MAPQ(const SeedVec &sv);

#endif /* ALNREAD_H */
