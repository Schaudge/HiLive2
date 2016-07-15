#ifndef TOOLS_H
#define TOOLS_H

#include "headers.h"
#include "definitions.h"
#include "kindex.h"

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
HashIntoType hash(const char * kmer, HashIntoType& _h, HashIntoType& _r);
HashIntoType hash_fw(const char * kmer, HashIntoType& _h);
HashIntoType rc(HashIntoType fw); 

// Update an existing kmer by left-shifting all nucleotides and appending new nucleotide
void update_kmer(HashIntoType &kmer, HashIntoType nuc);

// file name construction functions
std::string bcl_name(std::string rt, uint16_t ln, uint16_t tl, uint16_t cl);
std::string alignment_name(std::string rt, uint16_t ln, uint16_t tl, uint16_t cl);
std::string filter_name(std::string rt, uint16_t ln, uint16_t tl);
std::string position_name(std::string rt, uint16_t ln, uint16_t tl);
std::string sam_tile_name(std::string rt, uint16_t ln, uint16_t tl);
std::string sam_lane_name(std::string rt, uint16_t ln);





#endif /* TOOLS_H */
