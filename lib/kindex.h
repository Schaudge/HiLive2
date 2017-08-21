#ifndef KINDEX_H
#define KINDEX_H
 
#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "tools.h"
#include <seqan/sequence.h>
#include <seqan/index.h>



//-------------------------------------------------------------------//
//------  The k-mer index builder: KixBuild -------------------------//
//-------------------------------------------------------------------//

class KixBuild {

 public:
  
  // constructor resizing db (see below) to match the number of possible k-mers
  KixBuild();
  
  /**
   * Create an FM index.
   * @param fname Name of the input FASTA file
   * @param iname Name of the output index file
   * @param convert_spaces If true, spaces in the sequence headers are converted
   * @param trim_ids If true, sequence headers are trimed after the first white space
   * @return 0 if successful, error code otherwise.
   * @author Tobias Loka
   */
  int create_index( const std::string &fname, const std::string &iname, bool convert_spaces, bool trim_ids);

	/** Names of the sequences in the index. */
	StringListType seq_names;

	/** Lengths of the sequences in the index. */
	std::vector<uint32_t> seq_lengths;

	/** Sequences in the index. */
  seqan::StringSet<seqan::DnaString> seqs;

private:

  /**
   * Add a FASTA sequence to the index.
   * @param fname Name of the input FASTA file.
   * @param convert_spaces If true, spaces in the sequence headers are converted
   * @param trim_ids If true, sequence headers are trimed after the first white space
   * @return Number of loaded sequences
   * @author Tobias Loka
   */
  int add_fasta(const std::string &fname, bool convert_spaces, bool trim_ids);

  /**
   * Write the FM index to the hard drive.
   * @param idx The internal FM index object
   * @param iname Name of the output index file
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

	/** The FM index itself. */
	FMIndex idx;

	/**
	 * Load an FM-index from the hard drive.
	 * @param index_name Name of the input index file.
	 * @return 0 on success, other value on error
	 * @author Tobias Loka
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
	StringListType getSeqNames(){ return seq_names; };

	/**
	 * Get the sequence lengths vector
	 * @return The complete vector of sequence lengths
	 */
	std::vector<uint32_t> getSeqLengths(){ return seq_lengths; };

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

 private:

	/** Names of the sequences in the index. */
	StringListType seq_names;

	/** Lengths of the sequences in the index. */
	std::vector<uint32_t> seq_lengths;

};  // END class KixRun


#endif /* KINDEX_H */
