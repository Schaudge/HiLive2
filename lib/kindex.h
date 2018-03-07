#ifndef KINDEX_H
#define KINDEX_H
 
#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "tools.h"
#include <seqan/sequence.h>
#include <seqan/index.h>


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
}

inline bool operator<(const FMVertexDescriptor l, FMVertexDescriptor r) {
	if ( l.range.i1 == r.range.i1)
		return l.range.i2 < r.range.i2;
	return l.range.i1 < r.range.i1;
}

inline bool operator>(const FMVertexDescriptor l, FMVertexDescriptor r) {
	if ( l.range.i1 == r.range.i1)
		return l.range.i2 > r.range.i2;
	return l.range.i1 > r.range.i1;
}

inline bool operator!=(const FMVertexDescriptor l, FMVertexDescriptor r) {
	return !(l==r);
}


//-------------------------------------------------------------------//
//------  The k-mer index builder: KixBuild -------------------------//
//-------------------------------------------------------------------//

class KixBuild {

 public:

	/**
	 * Function to build the fm-index from a input (multi-)fasta file.
	 * @param fname Name of the input fasta file
	 * @param iname Name of the output index file(s)
	 * @param convert_spaces If true, inner spaces are converted to underscores
	 * @param trim_ids If true, the id is trimmed after the first whitespace
	 * @return 0 if successful, error code otherwise
	 * @author Tobias Loka
	 */
	int create_index( const std::string &fname, const std::string &iname, bool convert_spaces, bool trim_ids);

	/** Names of the sequences in the index. */
	std::vector<std::string> seq_names;

	/** Lengths of the sequences in the index. */
	std::vector<uint32_t> seq_lengths;

	/** Sequences in the index. */
	seqan::StringSet<seqan::DnaString> seqs;

 private:

	/**
	 * Read a fasta file and load the sequences to the HiLive data structure.
	 * @param fname Name of the input fasta file
	 * @param convert_spaces If true, inner spaces are converted to underscores
	 * @param trim_ids If true, the id is trimmed after the first whitespace
	 * @return Number of loaded sequences
	 * @author Martin Lindner, Tobias Loka
	 */
	int add_fasta(const std::string &fname, bool convert_spaces, bool trim_ids);

	/**
	 * Save the SeqAn FM-index to file.
	 * @param idx The SeqAn FM-index structure
	 * @param iname The name of the output index file(s)
	 * @return 0 if successful, error code otherwise
	 * @author Tobias Loka
	 */
	int save_fmindex(FMIndex & idx, const std::string &iname);

	/**
	 * Save the metadata to file.
	 * @param iname The name of the output index file(s)
	 * @return 0 if successful, error code otherwise
	 * @author Tobias Loka
	 */
	int save_metadata(const std::string &iname);

	/**
	 * Save the Sequence names to file.
	 * @param iname The name of the output index file(s)
	 * @return 0 if successful, error code otherwise
	 * @author Tobias Loka
	 */
	int save_seqnames(const std::string &iname);

	/**
	 * Save the Sequence lengths to file.
	 * @param iname The name of the output index file(s)
	 * @return 0 if successful, error code otherwise
	 * @author Tobias Loka
	 */
	int save_seqlengths(const std::string &iname);

};  // END class KixBuild


//-------------------------------------------------------------------//
//------  The k-mer runtime index: KixRun ---------------------------//
//-------------------------------------------------------------------//

class KixRun {

 public:

	/**
	 * Load the complete fm-index from file (incl. meta data)
	 * @param index_name: path to the index file(s)
	 * @return 0 on success, other value on failure
	 */
	int load_fmindex(std::string index_name);

	/**
	 * Load the meta data of the fm index
	 * @param iname Name of the index file(s)
	 * @return 0 on success, other value on failure
	 */
	int load_metadata( const std::string &iname );

	/**
	 * Load the sequence lengths of the fm index
	 * @param iname Name of the index file(s)
	 * @return 0 on success, other value on failure
	 */
	int load_seqlengths( const std::string &iname );

	/**
	 * Load the sequence names of the fm index
	 * @param iname Name of the index file(s)
	 * @return 0 on success, other value on failure
	 */
	int load_seqnames( const std::string &iname );

	/**
	 * Get the sequence name by using the fm index sequence id.
	 * !! The gid does NOT equal the index of the seq_names array !!
	 * @param gid Sequence id as given in the fm index
	 * @return Name of the sequence
	 */
	std::string getSequenceName(uint32_t gid) {return seq_names[gid/2]; };

	/**
	 * Get the sequence length by using the fm index sequence id.
	 * !! The gid does NOT equal the index of the seq_names array !!
	 * @param gid Sequence id as given in the fm index
	 * @return Length of the sequence
	 */
	uint32_t getSequenceLength(CountType gid) {return seq_lengths[gid/2]; };

	/**
	 * Get the sequence names list
	 * @return The complete list of sequence names
	 */
	std::vector<std::string> & getSeqNames(){ return seq_names; };

	/**
	 * Get the sequence lengths vector
	 * @return The complete vector of sequence lengths
	 */
	std::vector<uint32_t> & getSeqLengths(){ return seq_lengths; };

	/**
	 * Get the number of sequences
	 * @return The number of sequences in the index (real number of sequences, not the number of sequences in the fm index)
	 */
	CountType getNumSequences() { return seq_lengths.size(); };

	/**
	 * Check if an alignment is on the reverse strand.
	 * @param gid Sequence id as given in the fm index
	 * @return true, if alignment is on the reverse strand
	 */
	bool isReverse(uint32_t gid){ return gid%2==1; };

	/** The FM index itself. */
	FMIndex idx;

 private:

	/** Names of the sequences in the index. */
	std::vector<std::string> seq_names;

	/** Lengths of the sequences in the index. */
	std::vector<uint32_t> seq_lengths;

};  // END class KixRun


#endif /* KINDEX_H */
