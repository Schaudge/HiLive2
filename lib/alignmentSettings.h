#ifndef ALIGNMENTSETTINGS_H
#define ALIGNMENTSETTINGS_H

#include "definitions.h"

// Data structure to store the alignment settings
class AlignmentSettings {

private:

  // kmer gap positions
  Unmodifiable<std::vector<unsigned>> kmer_gaps;

  // reverse gap positions
  Unmodifiable<std::vector<unsigned>> rev_kmer_gaps;

  bool kmer_setFlag=false;

  // PARAMETER: kmer span (automatically computed from kmer_weight and kmer_gaps)
  Unmodifiable<uint8_t> kmer_span;

  // PARAMETER: Weight of the k-mers
  Unmodifiable<uint8_t> kmer_weight;

  // VARIABLE: maximum number of consecutive gaps in the gap pattern (will be computed at runtime)
  Unmodifiable<CountType> max_consecutive_gaps;

  // PARAMETER: Base Call quality cutoff, treat BC with quality < bc_cutoff as miscall
  Unmodifiable<CountType> min_qual;

  // PARAMETER: max. insert/deletion size
  Unmodifiable<DiffType> window;

  // PARAMETER: minimum number of errors allowed in alignment
  Unmodifiable<CountType> min_errors;

  // SWITCH: discard One-hit-wonders
  Unmodifiable<bool> discard_ohw;

  // PARAMETER: first cycle to discard one-hit-wonders
  Unmodifiable<CountType> start_ohw;

  // SWITCH: Any-Best-Hit-Mode
  Unmodifiable<bool> any_best_hit_mode;

  // SWITCH: Any-Best-Hit-Mode
  Unmodifiable<bool> all_best_hit_mode;

  // SWITCH: All-Best-N-Scores-Mode
  Unmodifiable<bool> all_best_n_scores_mode;

  // PARAMETER: All-Best-N-Scores-Mode::N
  Unmodifiable<CountType> best_n;

  // PARAMETER: temporary directory for the streamed alignment
  Unmodifiable<std::string> temp_dir;

  // SWITCH: write sam/bam output or not
  Unmodifiable<bool> write_bam=false;

  // SWITCH: Keep the old alignment files of previous cycles
  Unmodifiable<bool> keep_aln_files;

  // PARAMETER: Memory block size for the input and output buffer in the streamed alignment
  Unmodifiable<uint64_t> block_size;

  // PARAMETER: Compression format for alignment files
  Unmodifiable<uint8_t> compression_format;

  // PARAMETER: list of lanes to process
  Unmodifiable<std::vector<uint16_t>> lanes;
  
  // PARAMETER: list of tiles to process
  Unmodifiable<std::vector<uint16_t>> tiles;

  // PARAMETER: root directory of hilive run
  Unmodifiable<std::string> root;

  // PARAMETER: path to the index file
  Unmodifiable<std::string> index_fname;

  // PARAMETER: read length of all reads (including barcodes)
  Unmodifiable<CountType> cycles;

  // PARAMETER: read length of all reads (including barcodes)
  Unmodifiable<std::string> runInfo_fname;

  //PARAMETER: Stores the barcodes defined by the user. The inner vector contains the single fragments of multi-barcodes.
  Unmodifiable<std::vector<std::vector<std::string>>> barcodeVector;

  // PARAMETER: directory in which to create the output directory structure 
  Unmodifiable<boost::filesystem::path> out_dir;

  // PARAMETER: number of threads to use
  Unmodifiable<CountType> num_threads;

  // VARIABLE: list of trimmed reads for output generation
  std::vector<CountType> trimmedReads;

  // SWITCH: activate extended CIGAR annotation
  Unmodifiable<bool> extended_cigar;

  /**
   * Contains the read information of the sequencing machine (as SequenceElement objects). Includes sequence reads and barcodes.
   * Arbitrary numbers and orders of reads are supported. The summed length of all elements must equal the number of sequencing cycles.
   * @author Tobias Loka
   */
  Unmodifiable<std::vector<SequenceElement>> seqs;

  // Number of mates (information taken from the seqLengths parameter), (Hint: corresponding indeces are 1-based)
  Unmodifiable<uint16_t> mates;

  // PARAMETER: number of allowed errors for the single barcodes
  Unmodifiable<std::vector<uint16_t>> barcode_errors;

  // SWITCH: if true, keep all barcodes (disables barcode filtering).
  Unmodifiable<bool> keep_all_barcodes;

  template<typename T>
  void set_unmodifiable(Unmodifiable<T> & unmodifiable, T value, std::string variable_name) {
	  try {
		  unmodifiable.set(value);
	  }
	  catch (unmodifiable_error& e) {
		std::cerr << e.what() << " (" << variable_name << ")." << std::endl;
	  }
  }

  template<typename T>
  T get_unmodifiable(Unmodifiable<T> unmodifiable, std::string variable_name, bool allow_unset = false) {
	  try {
		  return unmodifiable.get(allow_unset);
	  }
	  catch (unmodifiable_error& e) {
		  std::cerr << e.what() << " (" << variable_name << ")." << std::endl;
		  return T();
	  }
  }

 public:
  /**
   * Get a SequenceElement object from the seqs vector by using the id
   * @param id The id of the SequenceElement.
   * @return The respective SequenceElement object for the given id.
   * @author Tobias Loka
   */
  SequenceElement getSeqById(CountType id) {return seqs.get()[id];}

  /**
   * Get a SequenceElement object from the seqs vector by using the mate number
   * @param id The mate number of the SequenceElement.
   * @return The respective SequenceElement object for the given mate number. NULLSEQ if mate==0 (barcodes).
   * @author Tobias Loka
   */
  SequenceElement getSeqByMate(CountType mate) {
	  if ( mate == 0 ) return NULLSEQ;
	  auto the_seq = seqs.get();
	  for (uint16_t i = 0; i != the_seq.size(); i++) {
		  if(the_seq[i].mate == mate) return the_seq[i];
	  }
	  return NULLSEQ;
  }

  std::vector<unsigned> get_kmer_gaps() {
      return get_unmodifiable(kmer_gaps, "kmer_gaps", true);
  }


  std::vector<unsigned> get_rev_kmer_gaps() {
      return get_unmodifiable(rev_kmer_gaps, "rev_kmer_gaps", true);
  }

  bool set_kmer( uint8_t kmer_weight, std::vector<unsigned> gaps ) {

	  if ( gaps.size() > 0 ) {
		  // Prepare user-defined list of gap positions (sort and erase duplicates)
		  std::sort(gaps.begin(), gaps.end());
		  gaps.erase( std::unique(gaps.begin(), gaps.end()), gaps.end());

		  // Weight and gap positions not consistent
		  if ( kmer_weight + gaps.size() <= *(std::max_element(gaps.begin(), gaps.end())) || *(std::min_element(gaps.begin(), gaps.end())) <= 1 ) {
			  std::cerr << "Warning: k-mer weight and gap pattern not consistent. Ensure that the first gap positions is >1 and" <<
					  "the maximal gap positions is lower than the total length of the k-mer pattern." << std::endl;
			  return false;
		  }
	  }

	  // Set k-mer variables
	  set_unmodifiable(this->kmer_weight, kmer_weight, "kmer_weight");
	  set_unmodifiable(this->kmer_gaps, gaps, "kmer_gaps");
	  set_unmodifiable(this->kmer_span, uint8_t(kmer_weight + gaps.size()), "kmer_span");

	  std::vector<unsigned> rev_kmer_gaps;
	  for ( auto gap:gaps ) {
		  rev_kmer_gaps.push_back(this->kmer_span - gap + 1);
	  }
	  std::reverse(rev_kmer_gaps.begin(), rev_kmer_gaps.end());
	  set_unmodifiable(this->rev_kmer_gaps, rev_kmer_gaps, "rev_kmer_gaps");

	  // Compute maximal consecutive gaps in gap pattern
	  CountType current_consecutive_gaps = 0;
	  CountType last_gap = 0;
	  CountType temp_max_consecutive_gaps = 0;

	  for ( unsigned el : this->get_kmer_gaps() ) {

		  // init first gap
		  if ( last_gap == 0 ) {
			  current_consecutive_gaps = 1;
			  last_gap = el;
			  continue;
		  }

		  // handle consecutive gaps
		  else if ( el == unsigned( last_gap + 1 ) ){
			  current_consecutive_gaps += 1;
			  last_gap = el;
		  }

		  // handle end of gap region
		  else {
			  temp_max_consecutive_gaps = std::max ( temp_max_consecutive_gaps, current_consecutive_gaps );
			  current_consecutive_gaps = 1;
			  last_gap = el;
		  }

	  }
	  set_unmodifiable(this->max_consecutive_gaps, std::max ( temp_max_consecutive_gaps, current_consecutive_gaps ), "max_consecutive_gaps");

	  return true;
  }

  uint8_t get_kmer_span() {
      return get_unmodifiable(kmer_span, "kmer_span");
  }

  uint8_t get_kmer_weight() {
      return get_unmodifiable(kmer_weight, "kmer_weight");
  }

  void set_min_qual(CountType value) {
    	  set_unmodifiable(min_qual, value, "min_qual");
  }

  CountType get_min_qual() {
      return get_unmodifiable(min_qual, "min_qual");
  }

  void set_window(DiffType value) {
      set_unmodifiable(window, value, "window");
  }

  DiffType get_window() {
      return get_unmodifiable(window, "window");
  }


  void set_min_errors(CountType value) {
	  set_unmodifiable(min_errors, value, "min_errors");
  }

  CountType get_min_errors() {
      return get_unmodifiable(min_errors, "min_errors");
  }


  void set_discard_ohw(bool value) {
	  set_unmodifiable(discard_ohw, value, "discard_ohw");
  }

  bool get_discard_ohw() {
      return get_unmodifiable(discard_ohw, "discard_ohw");
  }


  void set_start_ohw(CountType value) {
	  set_unmodifiable(start_ohw, value, "start_ohw");
  }

  CountType get_start_ohw() {
      return get_unmodifiable(start_ohw, "start_ohw");
  }

  void set_any_best_hit_mode(bool value) {
	  set_unmodifiable(any_best_hit_mode, value, "any_best_hit_mode");
  }

  bool get_any_best_hit_mode() {
      return get_unmodifiable(any_best_hit_mode, "any_best_hit_mode");
  }

  void set_all_best_hit_mode(bool value) {
	  set_unmodifiable(all_best_hit_mode, value, "all_best_hit_mode");
  }

  bool get_all_best_hit_mode() {
      return get_unmodifiable(all_best_hit_mode, "all_best_hit_mode");
  }


  void set_all_best_n_scores_mode(bool value) {
	  set_unmodifiable(all_best_n_scores_mode, value, "all_best_n_scores_mode");
  }

  bool get_all_best_n_scores_mode() {
      return get_unmodifiable(all_best_n_scores_mode, "all_best_n_scores_mode");
  }


  void set_best_n(CountType value) {
	  set_unmodifiable(best_n, value, "best_n");
  }

  CountType get_best_n() {
      return get_unmodifiable(best_n, "best_n");
  }

  void set_temp_dir(std::string value) {
      set_unmodifiable(temp_dir, value, "temp_dir");
  }

  std::string get_temp_dir() {
      return get_unmodifiable(temp_dir, "temp_dir");
  }


  void set_write_bam(bool value) {
	  set_unmodifiable(write_bam, value, "write_bam");
  }

  bool get_write_bam() {
      return get_unmodifiable(write_bam, "write_bam");
  }

  void set_keep_aln_files(bool value) {
	  set_unmodifiable(keep_aln_files, value, "keep_aln_files");
  }
  bool get_keep_aln_files() {
      return get_unmodifiable(keep_aln_files, "keep_aln_files");
  }

  void set_block_size(uint64_t value) {
	  set_unmodifiable(block_size, value, "block_size");
  }

  uint64_t get_block_size() {
      return get_unmodifiable(block_size, "block_size");
  }

  void set_compression_format(uint8_t value) {
	  set_unmodifiable(compression_format, value, "compression_format");
  }

  uint8_t get_compression_format() {
      return get_unmodifiable(compression_format, "compression_format");
  }


  void set_lanes(std::vector<uint16_t> value) {
	  set_unmodifiable(lanes, value, "lanes");
  }

  std::vector<uint16_t> get_lanes() {
      return get_unmodifiable(lanes, "lanes", true);
  }


  void set_tiles(std::vector<uint16_t> value) {
	  set_unmodifiable(tiles, value, "tiles");
  }

  std::vector<uint16_t> get_tiles() {
      return get_unmodifiable(tiles, "tiles", true);
  }


  void set_root(std::string value) {
	  set_unmodifiable(root, value, "root");
  }

  std::string get_root() {
      return get_unmodifiable(root, "root");
  }

  void set_index_fname(std::string value) {
	  set_unmodifiable(index_fname, value, "index_fname");
  }

  std::string get_index_fname() {
      return get_unmodifiable(index_fname, "index_fname");
  }

  void set_cycles(CountType value) {
	  set_unmodifiable(cycles, value, "cycles");
  }

  CountType get_cycles() {
      return get_unmodifiable(cycles, "cycles");
  }

  void set_runInfo_fname(std::string value) {
	  set_unmodifiable(runInfo_fname, value, "runInfo_fname");
  }

  std::string get_runInfo_fname() {
      return get_unmodifiable(runInfo_fname, "runInfo_fname");
  }

  void set_barcodeVector(std::vector<std::vector<std::string> > value) {
	  set_unmodifiable(barcodeVector, value, "barcodeVector");
  }

  std::vector<std::vector<std::string> > get_barcodeVector() {
      return get_unmodifiable(barcodeVector, "barcodeVector", true);
  }

  void set_out_dir(boost::filesystem::path value) {
	  set_unmodifiable(out_dir, value, "out_dir");
  }

  boost::filesystem::path get_out_dir() {
      return get_unmodifiable(out_dir, "out_dir");
  }


  void set_num_threads(CountType value) {
	  set_unmodifiable(num_threads, value, "num_threads");
  }

  CountType get_num_threads() {
      return get_unmodifiable(num_threads, "num_threads");
  }

  void set_seqs(std::vector<SequenceElement> value) {
	  set_unmodifiable(seqs, value, "seqs");
  }

  std::vector<SequenceElement> get_seqs() {
      return get_unmodifiable(seqs, "seqs", true);
  }

  void set_mates(uint16_t value) {
	  set_unmodifiable(mates, value, "mates");
  }

  uint16_t get_mates() {
      return get_unmodifiable(mates, "mates", true);
  }


  void set_barcode_errors(std::vector<uint16_t> value) {
	  set_unmodifiable(barcode_errors, value, "barcode_errors");
  }

  std::vector<uint16_t> get_barcode_errors() {
      return get_unmodifiable(barcode_errors, "barcode_errors", true);
  }

  void set_keep_all_barcodes(bool value) {
	  set_unmodifiable(keep_all_barcodes, value, "keep_all_barcodes");
  }

  bool get_keep_all_barcodes() {
	  if ( get_barcodeVector().size() == 0 )
		  return true;
      return get_unmodifiable(keep_all_barcodes, "keep_all_barcodes");
  }

  void set_trimmedReads(std::vector<CountType> value) {
	  trimmedReads = value;
  }

  std::vector<CountType> get_trimmedReads() {
      return trimmedReads;
  }

  // TODO this does not belong in alignmentSettings
  void add_trimmedRead(CountType value) {
      trimmedReads.push_back(value);
  }

  void set_extended_cigar(bool value) {
	  set_unmodifiable(extended_cigar, value, "extended_cigar");
  }

  bool get_extended_cigar() {
      return get_unmodifiable(extended_cigar, "extended_cigar");
  }

  CountType get_max_consecutive_gaps() {
      return get_unmodifiable(max_consecutive_gaps, "max_consecutive_gaps");
  }
};

#endif
