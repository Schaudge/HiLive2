#ifndef ALIGNMENTSETTINGS_H
#define ALIGNMENTSETTINGS_H

#include "definitions.h"

// all user parameters are stored in the alignment settings
class AlignmentSettings {
 private:
  // HARD CODED: kmer gap structure (this is not used anywhere, instead the kmer_gaps are parsed)
  std::string kmer_structure = "111111101111100111";
//  std::string kmer_structure = "111111111111111";

  // HARD CODED: kmer gap positions (one-based)
  std::vector<unsigned> kmer_gaps = {8, 14, 15};
  std::vector<unsigned> rev_kmer_gaps = {4, 5, 11};
//  std::vector<unsigned> kmer_gaps;
//  std::vector<unsigned> rev_kmer_gaps;

  // PARAMETER: kmer span (automatically computed from kmer_weight and kmer_gaps)
  uint8_t kmer_span;
  bool     kmer_span_setFlag=false;

  // PARAMETER: Weight of the k-mers
  uint8_t kmer_weight;
  bool    kmer_weight_setFlag=false;

  // VARIABLE: maximum number of consecutive gaps in the gap pattern (will be computed at runtime)
  CountType max_consecutive_gaps;
  bool max_consecutive_gaps_setFlag=false;

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
  CountType cycles;
  bool      cycles_setFlag=false;

  // PARAMETER: read length of all reads (including barcodes)
  std::string runInfo_fname;
  bool        runInfo_fname_setFlag=false;

  //PARAMETER: Stores the barcodes defined by the user. The inner vector contains the single fragments of multi-barcodes.
  std::vector<std::vector<std::string>> barcodeVector;
  bool                                  barcodeVector_setFlag=false;

  // PARAMETER: directory in which to create the output directory structure 
  boost::filesystem::path out_dir;
  bool                    out_dir_setFlag=false;

  // PARAMETER: number of threads to use
  CountType num_threads;
  bool      num_threads_setFlag=false;

  // VARIABLE: list of trimmed reads for output generation
  std::vector<CountType> trimmedReads;
  bool                   trimmedReads_setFlag=false;

  // SWITCH: activate extended CIGAR annotation
  bool extended_cigar;
  bool extended_cigar_setFlag=false;

  /**
   * Contains the read information of the sequencing machine (as SequenceElement objects). Includes sequence reads and barcodes.
   * Arbitrary numbers and orders of reads are supported. The summed length of all elements must equal the number of sequencing cycles.
   * @author Tobias Loka
   */
  std::vector<SequenceElement> seqs;
  bool                         seqs_setFlag=false;

  // Number of mates (information taken from the seqLengths parameter), (Hint: corresponding indeces are 1-based)
  uint16_t mates;
  bool     mates_setFlag=false;

  // PARAMETER: number of allowed errors for the single barcodes
  std::vector<uint16_t> barcode_errors;
  bool                  barcode_errors_setFlag=false;

  // SWITCH: if true, keep all barcodes (disables barcode filtering).
  bool keep_all_barcodes;
  bool keep_all_barcodes_setFlag=false;



 public:
  /**
   * Get a SequenceElement object from the seqs vector by using the id
   * @param id The id of the SequenceElement.
   * @return The respective SequenceElement object for the given id.
   * @author Tobias Loka
   */
  SequenceElement getSeqById(CountType id) {return seqs[id];}

  /**
   * Get a SequenceElement object from the seqs vector by using the mate number
   * @param id The mate number of the SequenceElement.
   * @return The respective SequenceElement object for the given mate number. NULLSEQ if mate==0 (barcodes).
   * @author Tobias Loka
   */
  SequenceElement getSeqByMate(CountType mate) {
	  if ( mate == 0 ) return NULLSEQ;
	  for (uint16_t i = 0; i != seqs.size(); i++) {
		  if(seqs[i].mate == mate) return seqs[i];
	  }
	  return NULLSEQ;
  }
 



 // getter and setter, all build up the same way, except the first four
  std::string get_kmer_structure() {
      return(this->kmer_structure);
  }


  std::vector<unsigned> get_kmer_gaps() {
      return(this->kmer_gaps);
  }


  std::vector<unsigned> get_rev_kmer_gaps() {
      return(this->rev_kmer_gaps);
  }


  void set_kmer_span(uint8_t value) {
      if (!kmer_span_setFlag) {
          kmer_span_setFlag = true;
          this->kmer_span = value;
      }
      else
          std::cerr << "Warning: kmer_span can only be set once." << std::endl;
  }
  uint8_t get_kmer_span() {
      return(this->kmer_span);
  }


  void set_kmer_weight(uint8_t value) {
      if (!kmer_weight_setFlag) {
          kmer_weight_setFlag = true;
          this->kmer_weight = value;
          this->set_kmer_span(this->kmer_weight + this->kmer_gaps.size());
      }
      else
          std::cerr << "Warning: kmer_weight can only be set once." << std::endl;
  }
  uint8_t get_kmer_weight() {
      return(this->kmer_weight);
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


  void set_write_bam(bool value) {
      if (!write_bam_setFlag) {
          write_bam_setFlag = true;
          this->write_bam = value;
      }
      else
          std::cerr << "Warning: write_bam can only be set once." << std::endl;
  }
  bool get_write_bam() {
      return(this->write_bam);
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


  void set_cycles(CountType value) {
      if (!cycles_setFlag) {
          cycles_setFlag = true;
          this->cycles = value;
      }
      else
          std::cerr << "Warning: cycles can only be set once." << std::endl;
  }
  CountType get_cycles() {
      return(this->cycles);
  }


  void set_runInfo_fname(std::string value) {
      if (!runInfo_fname_setFlag) {
          runInfo_fname_setFlag = true;
          this->runInfo_fname = value;
      }
      else
          std::cerr << "Warning: runInfo_fname can only be set once." << std::endl;
  }
  std::string get_runInfo_fname() {
      return(this->runInfo_fname);
  }


  void set_barcodeVector(std::vector<std::vector<std::string> > value) {
      if (!barcodeVector_setFlag) {
          barcodeVector_setFlag = true;
          this->barcodeVector = value;
      }
      else
          std::cerr << "Warning: barcodeVector can only be set once." << std::endl;
  }
  std::vector<std::vector<std::string> > get_barcodeVector() {
      return(this->barcodeVector);
  }


  void set_out_dir(boost::filesystem::path value) {
      if (!out_dir_setFlag) {
          out_dir_setFlag = true;
          this->out_dir = value;
      }
      else
          std::cerr << "Warning: out_dir can only be set once." << std::endl;
  }
  boost::filesystem::path get_out_dir() {
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


  void set_seqs(std::vector<SequenceElement> value) {
      if (!seqs_setFlag) {
          seqs_setFlag = true;
          this->seqs = value;
      }
      else
          std::cerr << "Warning: seqs can only be set once." << std::endl;
  }
  std::vector<SequenceElement> get_seqs() {
      return(this->seqs);
  }


  void set_mates(uint16_t value) {
      if (!mates_setFlag) {
          mates_setFlag = true;
          this->mates = value;
      }
      else
          std::cerr << "Warning: mates can only be set once." << std::endl;
  }
  uint16_t get_mates() {
      return(this->mates);
  }


  void set_barcode_errors(std::vector<uint16_t> value) {
      if (!barcode_errors_setFlag) {
          barcode_errors_setFlag = true;
          this->barcode_errors = value;
      }
      else
          std::cerr << "Warning: barcode_errors can only be set once." << std::endl;
  }
  std::vector<uint16_t> get_barcode_errors() {
      return(this->barcode_errors);
  }


  void set_keep_all_barcodes(bool value) {
      if (!keep_all_barcodes_setFlag) {
          keep_all_barcodes_setFlag = true;
          this->keep_all_barcodes = value;
      }
      else
          std::cerr << "Warning: keep_all_barcodes can only be set once." << std::endl;
  }
  bool get_keep_all_barcodes() {
      return(this->keep_all_barcodes);
  }


  void set_trimmedReads(std::vector<CountType> value) {
      if (!trimmedReads_setFlag) {
          trimmedReads_setFlag = true;
          this->trimmedReads = value;
      }
      else
          std::cerr << "Warning: trimmedReads can only be set once." << std::endl;
  }
  std::vector<CountType> get_trimmedReads() {
      return(this->trimmedReads);
  }
  // TODO this does not belong in alignmentSettings
  void add_trimmedRead(CountType value) {
      this->trimmedReads.push_back(value);
  }


  void set_extended_cigar(bool value) {
      if (!extended_cigar_setFlag) {
          extended_cigar_setFlag = true;
          this->extended_cigar = value;
      }
      else
          std::cerr << "Warning: extended_cigar can only be set once." << std::endl;
  }
  bool get_extended_cigar() {
      return(this->extended_cigar);
  }


  void set_max_consecutive_gaps(CountType value) {
      if (!max_consecutive_gaps_setFlag) {
          max_consecutive_gaps_setFlag = true;
          this->max_consecutive_gaps = value;
      }
      else
          std::cerr << "Warning: max_consecutive_gaps can only be set once." << std::endl;
  }
  CountType get_max_consecutive_gaps() {
      return(this->max_consecutive_gaps);
  }
};

#endif
