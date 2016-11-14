#ifndef TOOLS_H
#define TOOLS_H

#include "headers.h"
#include "definitions.h"
#include "kindex.h"


// compare function to sort GenomePosType objects by position
bool gp_compare (GenomePosType i,GenomePosType j);

// calculates the total size of a file in bytes
std::ifstream::pos_type get_filesize(const std::string &fname);

// checks if a directory with the given name exists
bool is_directory(const std::string &path);

// checks if a file exists
bool file_exists(const std::string &fname);

// reads a binary file from hdd and stores its raw content in a char vector
std::vector<char> read_binary_file(const std::string &fname);

// extract the number of reads from a BCL file
uint32_t num_reads_from_bcl(std::string bcl);

// writes a char vector into a binary file
uint64_t write_binary_file(const std::string &fname, const std::vector<char> & data);

//------  Hashing helper functions  ---------------------------------//
HashIntoType hash(const char * kmer, HashIntoType& _h, HashIntoType& _r, AlignmentSettings & settings);
std::string::const_iterator hash_fw(std::string::const_iterator it, std::string::const_iterator end, HashIntoType& _h, AlignmentSettings & settings);
//HashIntoType rc(HashIntoType fw); 
/* returns the sequence of a k-mer */
std::string unhash(HashIntoType myHash, unsigned hashLen=K_HiLive);

// file name construction functions
std::string bcl_name(std::string rt, uint16_t ln, uint16_t tl, uint16_t cl);
std::string alignment_name(uint16_t ln, uint16_t tl, uint16_t cl, uint16_t mt, std::string base);
std::string filter_name(std::string rt, uint16_t ln, uint16_t tl);
std::string position_name(std::string rt, uint16_t ln, uint16_t tl);
std::string sam_tile_name(std::string rt, uint16_t ln, uint16_t tl, uint16_t mate, bool write_bam);
std::string sam_lane_name(std::string rt, uint16_t ln, bool write_bam);

/** Get the current sequencing cycle using the current alignment cycle and read number.
 * @param cycle The read cycle.
 * @param settings Object containing the program settings.
 * @param read_number The read number (:= index of settings.seqs)
 * @return The sequencing cycle.
 * @author Tobias Loka
 */
uint16_t getSeqCycle(uint16_t cycle, AlignmentSettings* settings, uint16_t read_number=1);

/**
 * Split a string by a delimiter.
 * @param s String to split.
 * @param delim Delimiter for splitting
 * @param elems Vector to push the fragments
 * @return
 * @author Tobias Loka
 */
void split(const std::string &s, char delim, std::vector<std::string> &elems);


#endif /* TOOLS_H */
