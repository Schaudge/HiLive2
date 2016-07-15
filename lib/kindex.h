#ifndef KINDEX_H
#define KINDEX_H
 
#include "headers.h"
#include "definitions.h"
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
  int add_fasta(const std::string &fname, GenomeIdListType &ids);
  int add_fasta(const std::string &fname);
  
  // add all k-mers in a string sequence to the database
  GenomeIdType add_sequence(const std::string &s);

  // trim the database: remove kmers with more than max_count occurrences
  uint64_t trim(uint64_t max_count);

  // set the name of a sequence
  int set_name(GenomeIdType id, const std::string &name);

  // get the name of a sequence
  std::string get_name(GenomeIdType id);

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

};  // END class KIindex



//-------------------------------------------------------------------//
//------  The k-mer runtime index: KixRun ---------------------------//
//-------------------------------------------------------------------//


typedef std::array<char*,n_kmer> KixRunDB;

class KixRun {
 public:
  // pointer to the matching positions for a k-mer
  char* kmer(HashIntoType kmer);

  // get the name of a sequence
  std::string get_name(GenomeIdType id);
  
  // retrieve all fwd and rc occurrences of kmer in the index
  GenomePosListType retrieve_positions(HashIntoType kmer);

  // deserialize Kix
  uint64_t deserialize(char* d);
  
  // load and deserialize Kix from file
  uint64_t deserialize_file(std::string f);

  // generate a SAM compliant header string
  std::string get_SAM_header();

  // Database content
  GenomeIdType num_seq; // total number of sequences in the database
  StringListType seq_names; // names of the sequences in the database
  std::vector<uint32_t> seq_lengths; // lengths of the sequences in the database
  KixRunDB db;  // the lightweight database structure itself, pointing to sdata
  std::vector<char> sdata; // actual chunk of data
};  // END class KixRun










#endif /* KINDEX_H */
