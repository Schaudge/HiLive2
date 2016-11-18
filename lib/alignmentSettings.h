#ifndef ALIGNMENTSETTINGS_H
#define ALIGNMENTSETTINGS_H

// all user parameters are stored in the alignment settings
class AlignmentSettings {
 public:
  // HARD CODED: kmer gap structure (this is not used anywhere)
  //std::string kmer_structure = "11111110111110111";
  std::string kmer_structure = "111111111111111";

  // HARD CODED: kmer gap positions (one-based)
  //std::vector<unsigned> kmer_gaps = {8, 14};
  std::vector<unsigned> kmer_gaps;

  // HARD CODED: kmer span (kmer weight is K_HiLive)
  //unsigned kmer_span = K_HiLive+2;
  unsigned kmer_span = K_HiLive;

  // PARAMETER: Base Call quality cutoff, treat BC with quality < bc_cutoff as miscall
  CountType min_qual;

  // PARAMETER: max. insert/deletion size
  DiffType window;

  // PARAMETER: minimum number of errors allowed in alignment
  CountType min_errors;

  // SWITCH: discard One-hit-wonders
  bool discard_ohw;

  // PARAMETER: first cycle to discard one-hit-wonders
  CountType start_ohw;

  // SWITCH: Any-Best-Hit-Mode
  bool any_best_hit_mode;

  // SWITCH: Any-Best-Hit-Mode
  bool all_best_hit_mode;

  // SWITCH: All-Best-N-Scores-Mode
  bool all_best_n_scores_mode;

  // PARAMETER: All-Best-N-Scores-Mode::N
  CountType best_n;

  // PARAMETER: temporary directory for the streamed alignment
  std::string temp_dir;

  // SWITCH: write sam/bam output or not
  bool write_bam=false;

  // SWITCH: Keep the old alignment files of previous cycles
  bool keep_aln_files;

  // PARAMETER: Memory block size for the input and output buffer in the streamed alignment
  uint64_t block_size;

  // PARAMETER: Compression format for alignment files
  uint8_t compression_format;

  // PARAMETER: list of lanes to process
  std::vector<uint16_t> lanes;
  
  // PARAMETER: list of tiles to process
  std::vector<uint16_t> tiles;

  // PARAMETER: root directory of hilive run
  std::string root;

  // PARAMETER: path to the index file
  std::string index_fname;

  // PARAMETER: read length of all reads (including barcodes)
  CountType rlen;

  // PARAMETER: length of the sequence of all reads (excluding barcodes)
  CountType seqlen;

  // PARAMETER: vector containing all barcodes of the reads which should be outputted
  std::vector<std::string> barcodeVector;

  // PARAMETER: directory in which to create the output directory structure 
  std::string out_dir;

  // PARAMETER: number of threads to use
  CountType num_threads;
};

#endif
