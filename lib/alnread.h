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

	/** Vertex iterator for the FM index. Contains the position information of the alignment. */
	FMVertexDescriptor vDesc;

	/** Maximal AS:i score that can be reached based on the current alignment. */
	ScoreType max_as;

	/** Information about matches/mismatches (similar to CIGAR). The last element is the current one. */
	CigarVector cigar_data;

	/** Vector to store nucleotides of the index sequence when different to the read. */
	std::vector<uint8_t> mdz_nucleotides;

	/** Length of the mdz nucleotide list. */
	uint8_t mdz_length;

	/**
	 * Get the CIGAR string in SeqAn format
	 * @return CIGAR string in SeqAn format
	 * TODO: give extended CIGAR information as a parameter instead of inside the function.
	 */
	seqan::String<seqan::CigarElement<> > returnSeqanCigarString() const;

	/**
	 * Get the alignment score of a seed (AS:i).
	 * @return Alignment score.
	 */
	ScoreType get_as() const;

	/**
	 * Get the number of errors of a seed (NM:i).
	 * @return Number of errors.
	 */
	CountType get_nm() const;

	/**
	 * Get length of the front softclip.
	 * @return Length of the front softclip.
	 */
	CountType get_softclip_length() const;

	/**
	 * Add a nucleotide to the MDZ vector.
	 * @param nucl Nucleotide to be added to the MDZ vector (as char [ACGT])
	 */
	void add_mdz_nucleotide(char nucl);

	/**
	 * Get the String for the MD:Z tag.
	 * @return The MDZ tag as string.
	 */
	std::string getMDZString() const;

	/**
	 * Get all positions for this seed.
	 * @return A container containing the positions.
	 */
	std::vector<GenomePosType>  getPositions( ) const {
		return getPositions( 0, MAX_NUM_POSITIONS );
	}

	/**
	 * Get the positions for this seed. Returns all positions starting at a given index.
	 * @param firstPosition Index of the first position to obtain.
	 * @return A container containing the positions.
	 */
	std::vector<GenomePosType>  getPositions( CountType firstPosition ) const {
		return getPositions( firstPosition, MAX_NUM_POSITIONS );
	}

	/**
	 * Get the positions for this seed. The span is half-open: [firstPosition,lastPosition).
	 * @param firstPosition Index of the first position to obtain positions from.
	 * @param lastPosition Index of the last position to obtain positions from.
	 *        If the given index is higher than the last index, all positions until the last index are returned.
	 * @return A container containing the positions.
	 */
	std::vector<GenomePosType>  getPositions( CountType firstPosition, CountType lastPosition ) const;

	/**
	 * Get the number of positions where this seed aligns.
	 * @return The number of positions for the seed.
	 */
	CountType getNumPositions() const {
		return vDesc.range.i2 - vDesc.range.i1;
	}

	/**
	 * Determine size of the serialized seed
	 * @return size in bytes
	 */
	uint16_t serialize_size() const;

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

	/**
	 * [Debug] Print the seed to the command line.
	 */
	void cout() const;
};

/** Data type of USeed. */
typedef std::shared_ptr<Seed> USeed;

/** Shared pointer to seeds. */
typedef std::vector<USeed> SeedVec;

/** Iterator for a list of shared pointers to seeds. */
typedef SeedVec::iterator SeedVecIt;

/**
 * Compare two pointers by comparing the respective target objects.
 * @param a "left" pointer
 * @param b "right" pointer
 * @return true, if the object a points at is "smaller" (in terms of the comparison) than the object b points at.
 */
template <typename T> bool PComp(const T & a, const T & b)
{
	return *a < *b;
}

/**
 * Comparator for seeds to sort them by their alignment score.
 * @param a "left" USeed
 * @param b "right" USeed
 * @return true if the score of a is greater than the score of b. For equal scores,
 * the seed with no or shorter softclip (second criterion) or with less CIGAR
 * elements (third criterion) is ranked better.
 * @author Tobias Loka
 */
inline bool seed_comparison_by_as(const USeed a, const USeed b) {

	// if equal number of matches
	if (a->get_as() == b->get_as()) {

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

	return a->get_as() > b->get_as();
}

/**
 * Define '<'-operator for seeds.
 * Used for sorting by the range in the FM index.
 * @param l "left" Seed
 * @param r "right" Seed
 * @return true if the position of l in the index is smaller than the position of r.
 * @author Tobias Loka
 */
inline bool operator <(const Seed l, const Seed r) {

	if ( l.vDesc != r.vDesc ) {
		return l.vDesc < r.vDesc;
	}
	else {
		if ( l.get_softclip_length() != r.get_softclip_length() ) {
			return l.get_softclip_length() < r.get_softclip_length();
		}
		else {
			return l.max_as > r.max_as;
		}
	}

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
	void extendSeed(char base, USeed s, SeedVec & newSeeds);

	/**
	 * Extend a seed by alignment matches (Match or SNP)
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed to be extended
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void getMatchSeeds(CountType base_repr, USeed origin, SeedVec & newSeeds);

	/**
	 * Extend a seed by an insertion
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed to be extended
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void getInsertionSeeds(CountType base_repr, USeed origin, SeedVec & newSeeds);

	/**
	 * Extend a seed by deletions
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed to be extended
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void getDeletionSeeds(CountType base_repr, USeed origin, SeedVec & newSeeds);

	/**
	 * Add deletions to the alignment up to the permitted number of errors
	 * @param base_repr Binary representation of the current nucleotide
	 * @param origin The seed of the previous iteration
	 * @param allowedErrors Number of permitted errors for this seed
	 * @param index The FM index
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void recursive_goDown(CountType base_repr, USeed origin, SeedVec & newSeeds);

	/** Create new seeds
	 * @param index The FM index
	 * @param settings The alignment settings
	 * @param newSeeds Reference to the list of seeds (all resulting seeds are added to this list)
	 * @author Tobias Loka
	 */
	void createSeeds(SeedVec & newSeeds);


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
	uint64_t serialize_size() const;

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
	std::string getSequenceString() const;

	/**
	 * Convert and return sequence of the barcode. Multiple barcodes are concatenated (without delimiter).
	 * @return The Barcode as string
	 * @author Tobias Loka
	 */
	std::string getBarcodeString() const;

	/**
	 * Get base call quality of the read as string (without barcode).
	 * @return The read qualities in PHRED33 syntax.
	 */
	std::string getQualityString() const;

	/**
	 * Check whether the barcode of this read fulfills the criteria of at least one user-defined barcode.
	 * The nucleotides are only compared pairwise, not allowing for Indels.
	 * @return The index of the matching barcode in globalAlignmentSettings.multiBarcodeVector. NO_MATCH, if none.
	 * Also return NO_MATCH, if demultiplexing is not activated.
	 * @author 	Tobias Loka
	 */
	CountType getBarcodeIndex() const;

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
	void extend_alignment(char bc);

	/**
	 * Disable the alignment.
	 */
	void disable();

	/** Check if this read is disabled. */
	bool is_disabled() const;

	/**
	 * Obtain the start position of the alignment with SAM specifications (most left position)
	 * @param index The FM index
	 * @param p The position as stored in the index
	 * @param sd The seed containing the alignment information (CIGAR)
	 * @return The position in SAM specification
	 * @author Tobias Loka
	 */
	PositionType get_SAM_start_pos(GenomePosType p, USeed & sd) const;

	/**
	 * Sort the seeds by their alignment score.
	 * @author Tobias Loka
	 */
	void sort_seeds_by_as();

	/**
	 * Constructor for an object of type ReadAlignment
	 * @param tot_cyc Number of total alignment cycles
	 * @param cyc Current alignment cycle
	 * @author Tobias Loka
	 */
	ReadAlignment(CountType tot_cyc, CountType cyc){total_cycles=tot_cyc; cycle=cyc;};

	/**
	 * Get a vector of MAPQ values, one value for each seed.
	 * Size of the vector equals the number of seeds (not positions!).
	 * @return Vector containing the MAPQ values for each seed.
	 * @author Tobias Loka
	 */
	std::vector<uint8_t> getMAPQs() const;

	void addReadInfoToRecord(seqan::BamAlignmentRecord & record) const;

}; // END class ReadAlignment 



//-------------------------------------------------------------------//
//------  Other helper functions  -----------------------------------//
//-------------------------------------------------------------------//


#endif /* ALNREAD_H */
