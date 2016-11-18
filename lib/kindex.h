#ifndef KINDEX_H
#define KINDEX_H
 
#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "tools.h"



//-------------------------------------------------------------------//
//------  The k-mer index builder: KixBuild -------------------------//
//-------------------------------------------------------------------//

class KixBuild {
  // add a single k-mer to the database
  // Note: the index uses 1-based positions (to allow for negative positions)
  int add_kmer(HashIntoType kmer, GenomeIdType id, PositionType pos);

 public:
  
  // add k-mers of all sequences in FASTA file
  int add_fasta(const std::string &fname, GenomeIdListType &ids, bool convert_spaces, bool trim_ids);
  int add_fasta(const std::string &fname, bool convert_spaces, bool trim_ids);
  
  // add all k-mers in a string sequence to the database
  GenomeIdType start_sequence(const std::string &s, std::string& tailingKmer, PositionType& sequencePosition);
  GenomeIdType continue_sequence(const std::string &s, std::string& tailingKmer, PositionType& sequencePosition);

  // trim the database: remove kmers with more than max_count occurrences
  uint64_t trim(uint64_t max_count);

  // serialize the KixBuild
  std::vector<char> serialize();
  
  // serialize and store the KixBuild to a file
  uint64_t serialize_file(std::string f);

  // deserialize KixBuild
  uint64_t deserialize(char* d);
  
  // load and deserialize KixBuild from file
  uint64_t deserialize_file(std::string f);


  GenomeIdType num_seq; // total number of sequences in the database
  KmerIndexType db; // the database structure itself
  StringListType seq_names; // names of the sequences in the database
  std::vector<uint32_t> seq_lengths; // lengths of the sequences in the database

};  // END class KixBuild



//-------------------------------------------------------------------//
//------  The k-mer runtime index: KixRun ---------------------------//
//-------------------------------------------------------------------//


typedef std::array<char*,n_kmer> KixRunDB;

class KixRun {
 public:
  // pointer to the matching positions for a k-mer
  char* kmer(HashIntoType kmer);
  
  // retrieve all fwd and rc occurrences of kmer in the index
  GenomePosListType retrieve_positions(std::string kmerSpan);

  // deserialize Kix
  uint64_t deserialize(char* d);
  
  // load and deserialize Kix from file
  uint64_t deserialize_file(std::string f);

  // Database content
  GenomeIdType num_seq; // total number of sequences in the database
  StringListType seq_names; // names of the sequences in the database
  std::vector<uint32_t> seq_lengths; // lengths of the sequences in the database
  KixRunDB db;  // the lightweight database structure itself, pointing to sdata
  std::vector<char> sdata; // actual chunk of data
};  // END class KixRun


#endif /* KINDEX_H */
