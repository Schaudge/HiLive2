/**
 * This class provides functions that are dependent of the alignmentSettings!
 * For functions that are not dependent of any other HiLive class, please use the tools_static class!
 * Please do NOT add further includes to this file since this will lead to unwanted dependencies!
 */

#ifndef TOOLS_H
#define TOOLS_H

/* DONT ADD ANY INCLUDES */
#include "tools_static.h"
#include "alignmentSettings.h"
#include "global_variables.h"
/* DONT ADD ANY INCLUDES */


//////////////////////////////////////
////////// Build file names //////////
//////////////////////////////////////

/**
 * Get the name of a bcl file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @param cl The sequencing cycle.
 * @return Path to the bcl file.
 */
inline std::string get_bcl_fname(uint16_t ln, uint16_t tl, uint16_t cl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L" << to_N_digits(ln,3) << "/C" << cl << ".1/s_"<< ln <<"_" << tl << ".bcl";
  return path_stream.str();
}

/**
 * Get the name of an alignment file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @param cl The cycle for the respective mate.
 * @param mt The mate number.
 * @return Path to the alignment file.
 */
inline std::string get_align_fname(uint16_t ln, uint16_t tl, uint16_t cl, uint16_t mt){
  std::ostringstream path_stream;
  std::string base = globalAlignmentSettings.get_temp_dir() != "" ? globalAlignmentSettings.get_temp_dir() : globalAlignmentSettings.get_root();
  path_stream << base << "/L" << to_N_digits(ln,3) << "/s_"<< ln << "_" << tl << "." << mt << "."<< cl << ".align";
  return path_stream.str();
}
/**
 * Get the name of a filter file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @return Path to the filter file.
 */
inline std::string get_filter_fname(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L" << to_N_digits(ln,3) << "/s_"<< ln << "_" << tl << ".filter";
  return path_stream.str();
}
/**
 * Get the name of a clocs file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @return Path to the clocs file.
 */
inline std::string get_clocs_fname(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "../L" << to_N_digits(ln,3) << "/s_"<< ln << "_" << tl << ".clocs";
  return path_stream.str();
}
/**
 * Get the name of the settings file.
 * @return Path to the settings file.
 */
inline std::string get_config_fname() {
	std::ostringstream path_stream;
	path_stream << globalAlignmentSettings.get_temp_dir() << "/hilive_config.ini";
	return path_stream.str();
}
/** Get the current sequencing cycle using the current alignment cycle and read number.
 * @param cycle The read cycle.
 * @param seq_id The sequence id (:= id of the respective element in globalAlignmentSettings::seqs)
 * @return The sequencing cycle.
 * @author Tobias Loka
 */
inline uint16_t getSeqCycle(uint16_t cycle, uint16_t seq_id) {
	uint16_t seq_cycle = cycle;
	for ( int i = 0; i < seq_id; i++ )
		seq_cycle += globalAlignmentSettings.get_seq_by_id(i).length;
	return seq_cycle;
}
/**
 * Get the cycle of a mate for a given sequencing cycle.
 * When the mate is completely finished in the given cycle, return its total sequence length.
 * @param mate_number Mate of interest.
 * @param seq_cycle The sequencing cycle.
 * @return Cycle of the mate in the given sequencing cycle.
 * @author Tobias Loka
 */
inline uint16_t getMateCycle( uint16_t mate_number, uint16_t seq_cycle ) {

	// Invalid mate
	if ( mate_number == 0 || mate_number > globalAlignmentSettings.get_mates() )
		return 0;

	// Iterate through all sequence elements (including barcodes)
	for ( CountType id = 0; id < globalAlignmentSettings.get_seqs().size(); id++ ) {

		// Current sequence element
		SequenceElement seq = globalAlignmentSettings.get_seq_by_id(id);

		// Seq is mate of interest
		if ( seq.mate == mate_number )
			return ( seq.length > seq_cycle ? seq_cycle : seq.length );

		// Not enough cycles left to reach mate of interest
		else if ( seq.length >= seq_cycle )
			return 0;

		// Reduce number of cycles by the Seq length
		else
			seq_cycle -= seq.length;

	}

	// Should not be reached
	return 0;
}
////////////////////////////////////
////////// SAM/BAM output //////////
////////////////////////////////////

/**
 * Get the header for a SAM/BAM output file.
 * @return The BAM header.
 * @author Tobias Loka
 */
seqan::BamHeader getBamHeader();

/**
 * Name of a temporary SAM/BAM file (for the time it is written).
 * @param barcode Barcode of the output file (or "undetermined" for undetermined reads)
 * @param cycle The output cycle.
 * @return Name of the temporary output file for writing.
 * @author Tobias Loka
 */
inline std::string getBamTempFileName(std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_output_format() == OutputFormat::BAM ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_out_dir() << "/hilive_out_" << "cycle" << std::to_string(cycle) << "_" << barcode << ".temp" << file_suffix;
	return fname.str();
}

/**
 * Final name of a SAM/BAM file.
 * @param barcode Barcode of the output file (or "undetermined" for undetermined reads)
 * @param cycle The output cycle.
 * @return Name of the final output file.
 * @author Tobias Loka
 */
inline std::string getBamFileName(std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_output_format() == OutputFormat::BAM ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_out_dir() << "/hilive_out_" << "cycle" << std::to_string(cycle) << "_" << barcode << file_suffix;
	return fname.str();
}
/** Reverse a MD:Z tag for reverse alignments. */
std::string reverse_mdz(std::string mdz);


/////////////////////////////
////////// Scoring //////////
/////////////////////////////

inline uint16_t getMinSingleErrorPenalty() {

	// Mismatch: +1 mismatch, -1 match
	uint16_t mismatch_penalty = globalAlignmentSettings.get_mismatch_penalty() + globalAlignmentSettings.get_match_score();

	// Deletion: +1 deletion, maximum number of matches can still be reached
	uint16_t deletion_penalty = globalAlignmentSettings.get_deletion_opening_penalty() + globalAlignmentSettings.get_deletion_extension_penalty();

	// Insertion: +1 insertion, -1 match
	uint16_t insertion_penalty = globalAlignmentSettings.get_insertion_opening_penalty() + globalAlignmentSettings.get_insertion_extension_penalty() + globalAlignmentSettings.get_match_score();

	return std::min(mismatch_penalty, std::min(insertion_penalty, deletion_penalty));
}

inline uint16_t getMaxSingleErrorPenalty() {

	// Mismatch: +1 mismatch, -1 match
	uint16_t mismatch_penalty = globalAlignmentSettings.get_mismatch_penalty() + globalAlignmentSettings.get_match_score();

	// Deletion: +1 deletion, maximum number of matches can still be reached
	uint16_t deletion_penalty = globalAlignmentSettings.get_deletion_opening_penalty() + globalAlignmentSettings.get_deletion_extension_penalty();

	// Insertion: +1 insertion, -1 match
	uint16_t insertion_penalty = globalAlignmentSettings.get_insertion_opening_penalty() + globalAlignmentSettings.get_insertion_extension_penalty() + globalAlignmentSettings.get_match_score();

	return std::max(mismatch_penalty, std::max(insertion_penalty, deletion_penalty));
}

inline ScoreType getMaxPossibleScore( CountType cycles ) {
	return cycles * globalAlignmentSettings.get_match_score();
}

inline CountType getMinSoftclipPenalty( CountType softclip_length ) {
	return ceil( float(softclip_length) / globalAlignmentSettings.get_anchor_length() ) * getMinSingleErrorPenalty();
}

inline ScoreType getMinCycleScore( CountType cycle, CountType read_length ) {

	if ( cycle < globalAlignmentSettings.get_anchor_length() )
		return globalAlignmentSettings.get_min_as();

	ScoreType maxScore = getMaxPossibleScore(read_length);
	ScoreType minCycleScore = maxScore - ( ceil((cycle - globalAlignmentSettings.get_anchor_length()) / float(globalAlignmentSettings.get_error_rate())) * getMinSingleErrorPenalty() );
	return std::max(minCycleScore, globalAlignmentSettings.get_min_as());
}

/////////////////////////////////
////////// Other stuff //////////
/////////////////////////////////

/**
 * Copy a file while locking them in the global fileLocks.
 */
int atomic_rename( const char *oldname, const char *newname );

/**
 * Check if a cycle is a seeding cycle.
 */
inline bool isSeedingCycle(CountType cycle) {

	// Don't seed cycles smaller than the anchor length
	if ( cycle < globalAlignmentSettings.get_anchor_length() )
		return false;

	// Don't seed cycles larger than the maximal softclip length
	if ( cycle > globalAlignmentSettings.get_max_softclip_length())
		return false;

	// Create seeds when reaching the anchor length for the first time
	if ( cycle == globalAlignmentSettings.get_anchor_length() )
		return true;

	// Create seeds every seeding_interval cycle after the first anchor
	if ( ( cycle - globalAlignmentSettings.get_anchor_length() ) % globalAlignmentSettings.get_seeding_interval() == 0 )
		return true;

	return false;

}

#endif /* TOOLS_H */
