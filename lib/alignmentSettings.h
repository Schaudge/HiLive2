#ifndef ALIGNMENTSETTINGS_H
#define ALIGNMENTSETTINGS_H

// all user parameters are stored in the alignment settings
class AlignmentSettings {
 private:
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
  bool      min_qual_setFlag=false;

  // PARAMETER: max. insert/deletion size
  DiffType window;
  bool     window_setFlag=false;

  // PARAMETER: minimum number of errors allowed in alignment
  CountType min_errors;
  bool      min_errors_setFlag=false;

  // SWITCH: discard One-hit-wonders
  bool discard_ohw;
  bool discard_ohw_setFlag=false;

  // PARAMETER: first cycle to discard one-hit-wonders
  CountType start_ohw;
  bool      start_ohw_setFlag=false;

  // SWITCH: Any-Best-Hit-Mode
  bool any_best_hit_mode;
  bool any_best_hit_mode_setFlag=false;

  // SWITCH: Any-Best-Hit-Mode
  bool all_best_hit_mode;
  bool all_best_hit_mode_setFlag=false;

  // SWITCH: All-Best-N-Scores-Mode
  bool all_best_n_scores_mode;
  bool all_best_n_scores_mode_setFlag=false;

  // PARAMETER: All-Best-N-Scores-Mode::N
  CountType best_n;
  bool      best_n_setFlag=false;

  // PARAMETER: temporary directory for the streamed alignment
  std::string temp_dir;
  bool        temp_dir_setFlag=false;

  // SWITCH: write sam/bam output or not
  bool write_bam=false;
  bool write_bam_setFlag=false;

  // SWITCH: Keep the old alignment files of previous cycles
  bool keep_aln_files;
  bool keep_aln_files_setFlag=false;

  // PARAMETER: Memory block size for the input and output buffer in the streamed alignment
  uint64_t block_size;
  bool     block_size_setFlag=false;

  // PARAMETER: Compression format for alignment files
  uint8_t compression_format;
  bool    compression_format_setFlag=false;

  // PARAMETER: list of lanes to process
  std::vector<uint16_t> lanes;
  bool                  lanes_setFlag=false;
  
  // PARAMETER: list of tiles to process
  std::vector<uint16_t> tiles;
  bool                  tiles_setFlag=false;

  // PARAMETER: root directory of hilive run
  std::string root;
  bool        root_setFlag=false;

  // PARAMETER: path to the index file
  std::string index_fname;
  bool        index_fname_setFlag=false;

  // PARAMETER: read length of all reads (including barcodes)
  CountType rlen;
  bool      rlen_setFlag=false;

  // PARAMETER: length of the sequence of all reads (excluding barcodes)
  CountType seqlen;
  bool      seqlen_setFlag=false;

  // PARAMETER: vector containing all barcodes of the reads which should be outputted
  std::vector<std::string> barcodeVector;
  bool                     barcodeVector_setFlag=false;

  // PARAMETER: directory in which to create the output directory structure 
  std::string out_dir;
  bool        out_dir_setFlag=false;

  // PARAMETER: number of threads to use
  CountType num_threads;
  bool      num_threads_setFlag=false;

 



 // getter and setter, all build up the same way, except the first three, hardcoded ones
 public:
  std::string get_kmer_structure() {
      return(this->kmer_structure);
  }


  std::vector<unsigned> get_kmer_gaps() {
      return(this->kmer_gaps);
  }


  unsigned get_kmer_span() {
      return(this->kmer_span);
  }


  void set_min_qual(CountType value) {
      if (!min_qual_setFlag) {
          min_qual_setFlag = true;
          this->min_qual = value;
      }
      else
          std::cerr << "Warning: min_qual can only be set once." << std::endl;
  }
  CountType get_min_qual() {
      return(this->min_qual);
  }


  void set_window(DiffType value) {
      if (!window_setFlag) {
          window_setFlag = true;
          this->window = value;
      }
      else
          std::cerr << "Warning: window can only be set once." << std::endl;
  }
  DiffType get_window() {
      return(this->window);
  }


  void set_min_errors(CountType value) {
      if (!min_errors_setFlag) {
          min_errors_setFlag = true;
          this->min_errors = value;
      }
      else
          std::cerr << "Warning: min_errors can only be set once." << std::endl;
  }
  CountType get_min_errors() {
      return(this->min_errors);
  }


  void set_discard_ohw(bool value) {
      if (!discard_ohw_setFlag) {
          discard_ohw_setFlag = true;
          this->discard_ohw = value;
      }
      else
          std::cerr << "Warning: discard_ohw can only be set once." << std::endl;
  }
  bool get_discard_ohw() {
      return(this->discard_ohw);
  }


  void set_start_ohw(CountType value) {
      if (!start_ohw_setFlag) {
          start_ohw_setFlag = true;
          this->start_ohw = value;
      }
      else
          std::cerr << "Warning: start_ohw can only be set once." << std::endl;
  }
  CountType get_start_ohw() {
      return(this->start_ohw);
  }


  void set_any_best_hit_mode(bool value) {
      if (!any_best_hit_mode_setFlag) {
          any_best_hit_mode_setFlag = true;
          this->any_best_hit_mode = value;
      }
      else
          std::cerr << "Warning: any_best_hit_mode can only be set once." << std::endl;
  }
  bool get_any_best_hit_mode() {
      return(this->any_best_hit_mode);
  }


  void set_all_best_hit_mode(bool value) {
      if (!all_best_hit_mode_setFlag) {
          all_best_hit_mode_setFlag = true;
          this->all_best_hit_mode = value;
      }
      else
          std::cerr << "Warning: all_best_hit_mode can only be set once." << std::endl;
  }
  bool get_all_best_hit_mode() {
      return(this->all_best_hit_mode);
  }


  void set_all_best_n_scores_mode(bool value) {
      if (!all_best_n_scores_mode_setFlag) {
          all_best_n_scores_mode_setFlag = true;
          this->all_best_n_scores_mode = value;
      }
      else
          std::cerr << "Warning: all_best_n_scores_mode can only be set once." << std::endl;
  }
  bool get_all_best_n_scores_mode() {
      return(this->all_best_n_scores_mode);
  }


  void set_best_n(CountType value) {
      if (!best_n_setFlag) {
          best_n_setFlag = true;
          this->best_n = value;
      }
      else
          std::cerr << "Warning: best_n can only be set once." << std::endl;
  }
  CountType get_best_n() {
      return(this->best_n);
  }


  void set_temp_dir(std::string value) {
      if (!temp_dir_setFlag) {
          temp_dir_setFlag = true;
          this->temp_dir = value;
      }
      else
          std::cerr << "Warning: temp_dir can only be set once." << std::endl;
  }
  std::string get_temp_dir() {
      return(this->temp_dir);
  }


  void set_write_bam=false(bool value) {
      if (!write_bam=false_setFlag) {
          write_bam=false_setFlag = true;
          this->write_bam=false = value;
      }
      else
          std::cerr << "Warning: write_bam=false can only be set once." << std::endl;
  }
  bool get_write_bam=false() {
      return(this->write_bam=false);
  }


  void set_keep_aln_files(bool value) {
      if (!keep_aln_files_setFlag) {
          keep_aln_files_setFlag = true;
          this->keep_aln_files = value;
      }
      else
          std::cerr << "Warning: keep_aln_files can only be set once." << std::endl;
  }
  bool get_keep_aln_files() {
      return(this->keep_aln_files);
  }


  void set_block_size(uint64_t value) {
      if (!block_size_setFlag) {
          block_size_setFlag = true;
          this->block_size = value;
      }
      else
          std::cerr << "Warning: block_size can only be set once." << std::endl;
  }
  uint64_t get_block_size() {
      return(this->block_size);
  }


  void set_compression_format(uint8_t value) {
      if (!compression_format_setFlag) {
          compression_format_setFlag = true;
          this->compression_format = value;
      }
      else
          std::cerr << "Warning: compression_format can only be set once." << std::endl;
  }
  uint8_t get_compression_format() {
      return(this->compression_format);
  }


  void set_lanes(std::vector<uint16_t> value) {
      if (!lanes_setFlag) {
          lanes_setFlag = true;
          this->lanes = value;
      }
      else
          std::cerr << "Warning: lanes can only be set once." << std::endl;
  }
  std::vector<uint16_t> get_lanes() {
      return(this->lanes);
  }


  void set_tiles(std::vector<uint16_t> value) {
      if (!tiles_setFlag) {
          tiles_setFlag = true;
          this->tiles = value;
      }
      else
          std::cerr << "Warning: tiles can only be set once." << std::endl;
  }
  std::vector<uint16_t> get_tiles() {
      return(this->tiles);
  }


  void set_root(std::string value) {
      if (!root_setFlag) {
          root_setFlag = true;
          this->root = value;
      }
      else
          std::cerr << "Warning: root can only be set once." << std::endl;
  }
  std::string get_root() {
      return(this->root);
  }


  void set_index_fname(std::string value) {
      if (!index_fname_setFlag) {
          index_fname_setFlag = true;
          this->index_fname = value;
      }
      else
          std::cerr << "Warning: index_fname can only be set once." << std::endl;
  }
  std::string get_index_fname() {
      return(this->index_fname);
  }


  void set_rlen(CountType value) {
      if (!rlen_setFlag) {
          rlen_setFlag = true;
          this->rlen = value;
      }
      else
          std::cerr << "Warning: rlen can only be set once." << std::endl;
  }
  CountType get_rlen() {
      return(this->rlen);
  }


  void set_seqlen(CountType value) {
      if (!seqlen_setFlag) {
          seqlen_setFlag = true;
          this->seqlen = value;
      }
      else
          std::cerr << "Warning: seqlen can only be set once." << std::endl;
  }
  CountType get_seqlen() {
      return(this->seqlen);
  }


  void set_barcodeVector(std::vector<std::string> value) {
      if (!barcodeVector_setFlag) {
          barcodeVector_setFlag = true;
          this->barcodeVector = value;
      }
      else
          std::cerr << "Warning: barcodeVector can only be set once." << std::endl;
  }
  std::vector<std::string> get_barcodeVector() {
      return(this->barcodeVector);
  }


  void set_out_dir(std::string value) {
      if (!out_dir_setFlag) {
          out_dir_setFlag = true;
          this->out_dir = value;
      }
      else
          std::cerr << "Warning: out_dir can only be set once." << std::endl;
  }
  std::string get_out_dir() {
      return(this->out_dir);
  }


  void set_num_threads(CountType value) {
      if (!num_threads_setFlag) {
          num_threads_setFlag = true;
          this->num_threads = value;
      }
      else
          std::cerr << "Warning: num_threads can only be set once." << std::endl;
  }
  CountType get_num_threads() {
      return(this->num_threads);
  }
};

#endif
