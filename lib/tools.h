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


///////////////////////////////////
////////// K-mer Hashing //////////
///////////////////////////////////

/**
 * Calculate the sequence from a hash value.
 * @param myHash The input hash value.
 * @param hashLen Length (weight) of the hashed sequence.
 * @return The unhashed sequence.
 */
std::string unhash(HashIntoType myHash, unsigned hashLen);


////////////////////////////////////////////
////////// File name construction //////////
////////////////////////////////////////////

/**
 * Get the name of a bcl file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @param cl The sequencing cycle.
 * @return Path to the bcl file.
 */
std::string bcl_name(uint16_t ln, uint16_t tl, uint16_t cl);

/**
 * Get the name of an alignment file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @param cl The cycle for the respective mate.
 * @param mt The mate number.
 * @return Path to the alignment file.
 */
std::string alignment_name(uint16_t ln, uint16_t tl, uint16_t cl, uint16_t mt);

/**
 * Get the name of a filter file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @return Path to the filter file.
 */
std::string filter_name(uint16_t ln, uint16_t tl);

/**
 * Get the name of a clocs file.
 * @param ln The lane number.
 * @param tl The tile number.
 * @return Path to the clocs file.
 */
std::string position_name(uint16_t ln, uint16_t tl);

/**
 * Get the name of the settings file.
 * @return Path to the settings file.
 */
std::string get_settings_name();

/**
 * Get the name of the output log file.
 * @return Path to the output log file.
 */
std::string get_out_log_name();

/** Get the current sequencing cycle using the current alignment cycle and read number.
 * @param cycle The read cycle.
 * @param seq_id The sequence id (:= id of the respective element in globalAlignmentSettings::seqs)
 * @return The sequencing cycle.
 * @author Tobias Loka
 */
uint16_t getSeqCycle(uint16_t cycle, uint16_t seq_id=1);

/**
 * Get the cycle of a mate for a given sequencing cycle.
 * When the mate is completely finished in the given cycle, return its total sequence length.
 * @param mate_number Mate of interest.
 * @param seq_cycle The sequencing cycle.
 * @return Cycle of the mate in the given sequencing cycle.
 * @author Tobias Loka
 */
uint16_t getMateCycle( uint16_t mate_number, uint16_t seq_cycle );


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
std::string getBamTempFileName(std::string barcode, CountType cycle);

/**
 * Final name of a SAM/BAM file.
 * @param barcode Barcode of the output file (or "undetermined" for undetermined reads)
 * @param cycle The output cycle.
 * @return Name of the final output file.
 * @author Tobias Loka
 */
std::string getBamFileName(std::string barcode, CountType cycle);

/** Reverse a MD:Z tag for reverse alignments. */
std::string reverse_mdz(std::string mdz);


/////////////////////////////
////////// Scoring //////////
/////////////////////////////

uint16_t getMinSingleErrorPenalty();

uint16_t getMaxSingleErrorPenalty();

uint16_t getMaxPossibleScore( CountType cycles );

uint16_t getMinSoftclipPenalty( CountType softclip_length );

ScoreType getMinCycleScore( CountType cycle, CountType read_length );

#endif /* TOOLS_H */
