/**
 * This class provides functions that are dependent of the alignmentSettings!
 * For functions that are not dependent of any other HiLive class, please use the tools_static class!
 * Please do NOT add further includes to this file since this will lead to unwanted dependencies!
 */

#ifndef TOOLS_H
#define TOOLS_H

/* DONT ADD ANY INCLUDES */
#include "headers.h"
#include "definitions.h"
#include "tools_static.h"
#include "alignmentSettings.h"
#include "global_variables.h"
#include "kindex.h" // TODO: remove dependency. Functions that are dependent on kindex should be placed there.
/* DONT ADD ANY INCLUDES */


//------  Hashing helper functions  ---------------------------------//
HashIntoType hash(const char * kmer, HashIntoType& _h, HashIntoType& _r);
std::string::const_iterator hash_fw(std::string::const_iterator it, std::string::const_iterator end, HashIntoType& _h);
//HashIntoType rc(HashIntoType fw); 
/* returns the sequence of a k-mer */
std::string unhash(HashIntoType myHash, unsigned hashLen=globalAlignmentSettings.get_kmer_weight());

// file name construction functions
std::string bcl_name(uint16_t ln, uint16_t tl, uint16_t cl);
std::string alignment_name(uint16_t ln, uint16_t tl, uint16_t cl, uint16_t mt);
std::string filter_name(uint16_t ln, uint16_t tl);
std::string position_name(uint16_t ln, uint16_t tl);
std::string sam_tile_name(uint16_t ln, uint16_t tl, bool write_bam);
std::string get_xml_out_name();

/** Get the current sequencing cycle using the current alignment cycle and read number.
 * @param cycle The read cycle.
 * @param read_number The read number (:= index of globalAlignmentSettings.seqs)
 * @return The sequencing cycle.
 * @author Tobias Loka
 */
uint16_t getSeqCycle(uint16_t cycle, uint16_t read_number=1);

#endif /* TOOLS_H */
