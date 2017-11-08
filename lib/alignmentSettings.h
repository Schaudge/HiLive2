#ifndef ALIGNMENTSETTINGS_H
#define ALIGNMENTSETTINGS_H

#include "headers.h"
#include "definitions.h"
#include "tools_static.h"

// Data structure to store the alignment settings
class AlignmentSettings {

private:

  // kmer gap positions
  Unmodifiable<std::vector<unsigned>> kmer_gaps;

  // reverse gap positions
  Unmodifiable<std::vector<unsigned>> rev_kmer_gaps;

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

  // PARAMETER: All-Best-N-Scores-Mode::N
  Unmodifiable<CountType> best_n;

  // PARAMETER: temporary directory for the streamed alignment
  Unmodifiable<std::string> temp_dir;

  // SWITCH: write sam/bam output or not
  Unmodifiable<bool> write_bam;

  // PARAMETER: Cycles for intermediate SAM/BAM output
  Unmodifiable<std::vector<uint16_t>> output_cycles;

  // SWITCH: Keep the old alignment files of previous cycles
  Unmodifiable<std::vector<uint16_t>> keep_aln_files;

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

  // PARAMETER: the first cycle to handle. Should be 1 by default.
  Unmodifiable<CountType> start_cycle;

  // PARAMETER: read length of all reads (including barcodes)
  Unmodifiable<CountType> cycles;

  //PARAMETER: Stores the barcodes defined by the user. The inner vector contains the single fragments of multi-barcodes.
  Unmodifiable<std::vector<std::vector<std::string>>> barcodeVector;

  // PARAMETER: directory in which to create the output directory structure 
  Unmodifiable<std::string> out_dir;

  // PARAMETER: number of threads to use
  Unmodifiable<CountType> num_threads;

  // PARAMETER: max. amount of threads used for output
  Unmodifiable<CountType> num_out_threads;

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

  Unmodifiable<AlignmentMode> mode;

  Unmodifiable<float> min_as_ratio;

  template<typename T>
  bool set_unmodifiable(Unmodifiable<T> & unmodifiable, T value, std::string variable_name) {
	  try {
		  unmodifiable.set(value);
	  }
	  catch (unmodifiable_error& e) {
//		std::cerr << e.what() << " (" << variable_name << ")." << std::endl;
		  variable_name.length(); // TODO: just to remove warnings. Remove variable_name string when finished.
		 return false ;
	  }
	  return true;
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

  std::vector<std::string> xmlParse_barcodeVector() {
	  std::vector<std::string> bc_strings;
	  for ( CountType i = 0; i < get_barcodeVector().size(); i++ ) {
		  bc_strings.push_back( get_barcodeString(i) );
	  }
	  return bc_strings;
  }

  std::vector<std::string> xmlParse_seqs() {
	  std::vector<std::string> seq_vector;
	  for ( auto el : get_seqs() ) {
		  std::string seq_string;
		  seq_string += std::to_string(el.length);
		  seq_string += el.mate == 0 ? "B" : "R";
		  seq_vector.push_back(seq_string);
	  }
	  return seq_vector;
  }

  void set_barcodeVector(std::vector<std::vector<std::string> > value) {
	  set_unmodifiable(barcodeVector, value, "barcodeVector");
  }

  void set_seqs(std::vector<SequenceElement> value) {
  	  set_unmodifiable(seqs, value, "seqs");
    }

  void set_block_size(uint64_t value) {
	  set_unmodifiable(block_size, value, "block_size");
  }

  void set_mates(uint16_t value) {
  	  set_unmodifiable(mates, value, "mates");
    }

  AlignmentMode get_mode() {
	  return get_unmodifiable(mode, "mode");
  }

  void set_mode(AlignmentMode value, CountType bestn = 0) {
	  set_unmodifiable(mode, value, "mode");
	  set_unmodifiable(best_n, bestn, "best_n");
  }

public:

  /**
   * Create a property tree that is filled with all (relevant) settings.
   * @return Property tree containing all settings.
   * @author Tobias Loka
   */
  boost::property_tree::ptree to_ptree() {

	  boost::property_tree::ptree xml_out;

	  // General settings
	  xml_out.add_child("settings.lanes", getXMLnode_vector ( get_lanes() ));
	  xml_out.add_child("settings.tiles", getXMLnode_vector ( get_tiles() ));
	  xml_out.add_child("settings.min_errors", getXMLnode (get_min_errors() ));
	  xml_out.add_child("settings.cycles", getXMLnode ( get_cycles() ));
	  xml_out.add_child("settings.sequences", getXMLnode_vector ( xmlParse_seqs() ));

	  // Barcode settings
	  xml_out.add_child("settings.barcodes.sequences", getXMLnode_vector( xmlParse_barcodeVector() ));
	  xml_out.add_child("settings.barcodes.errors", getXMLnode_vector ( get_barcode_errors() ));
	  xml_out.add_child("settings.barcodes.keep_all", getXMLnode ( get_keep_all_barcodes() ));

	  // Alignment mode
	  std::string mode = std::string(1, char(get_mode()));
	  if ( get_mode() == AlignmentMode::BESTN )
		  mode += std::to_string(get_best_n());
	  xml_out.add_child("settings.mode", getXMLnode ( mode ));

	  // Paths
	  xml_out.add_child("settings.paths.temp_dir", getXMLnode ( get_temp_dir() ));
	  xml_out.add_child("settings.paths.out_dir", getXMLnode ( get_out_dir() ));
	  xml_out.add_child("settings.paths.root", getXMLnode ( get_root() ));
	  xml_out.add_child("settings.paths.index", getXMLnode ( get_index_fname() ));

	  // Output settings
	  xml_out.add_child("settings.out.bam", getXMLnode ( get_write_bam() ));
	  xml_out.add_child("settings.out.cycles", getXMLnode_vector ( get_output_cycles() ));
	  xml_out.add_child("settings.out.extended_cigar", getXMLnode ( get_extended_cigar() ));
	  xml_out.add_child("settings.out.min_as_ratio", getXMLnode ( get_min_as_ratio()) );

	  // Technical settings
	  xml_out.add_child("settings.technical.num_threads", getXMLnode ( get_num_threads() ));
	  xml_out.add_child("settings.technical.num_out_threads", getXMLnode ( get_num_out_threads() ));
	  xml_out.add_child("settings.technical.keep_aln_files", getXMLnode_vector ( get_keep_aln_files() ));
	  xml_out.add_child("settings.technical.block_size", getXMLnode ( get_block_size() ));
	  xml_out.add_child("settings.technical.compression_format", getXMLnode ( get_compression_format() ));

	  // Alignment algorithm settings
	  xml_out.add_child("settings.align.min_qual", getXMLnode (get_min_qual() ));
	  xml_out.add_child("settings.align.window", getXMLnode (get_window() ));
	  xml_out.add_child("settings.align.discard_ohw", getXMLnode ( get_discard_ohw() ));
	  xml_out.add_child("settings.align.start_ohw", getXMLnode ( get_start_ohw() ));

	  return xml_out;
  }

 void set_barcodes ( std::vector< std::string > barcodeArg ) {

	 // Get the barcode length(s) from the seqs vector
  	std::vector<uint16_t> barcode_lengths;
  	for ( uint16_t seq_num = 0; seq_num < get_seqs().size(); seq_num++ ) {
  		if ( getSeqById(seq_num).isBarcode() )
  	  		barcode_lengths.push_back( getSeqById(seq_num).length );
  	}

  	// Fill 2D-vector for internal storage
    std::vector<std::vector<std::string> > barcodeVector;
  	for ( auto barcode = barcodeArg.begin(); barcode != barcodeArg.end(); ++barcode) {

  		// Check if all characters in the current barcode are valid
  		std::string valid_chars = seq_chars + "-";
  		for(CountType i = 0; i != (*barcode).length(); i++){
  			char c = (*barcode)[i];
  			if ( valid_chars.find(c) == std::string::npos )
  				throw std::runtime_error("Invalid character '" + std::string(1,c) + "' in barcode sequence " + *barcode + ".");
  		}

  		// Split barcode string into fragments
  		std::vector<std::string> fragments;
  		split(*barcode, '-', fragments);

  		// Check correct length of all fragments
  		if ( barcode_lengths.size() != fragments.size())
			throw std::runtime_error("Wrong fragment length in barcode " + *barcode);

  		// Check correct number of fragments
  		for ( uint16_t num = 0; num != fragments.size(); num++ ) {
  			if ( fragments[num].length() != barcode_lengths[num] ) {
  				throw std::runtime_error("Wrong number of fragments for barcode " + *barcode);
  			}
  		}

  		// Push barcode to the final 2D vector
  		barcodeVector.push_back(fragments);
  	}

  	set_barcodeVector(barcodeVector);

  }

 std::vector<std::vector<std::string> > get_barcodeVector() {
     return get_unmodifiable(barcodeVector, "barcodeVector", true);
 }

 std::string format_barcode(std::string unformatted_barcode) {
	  CountType pos = 0;
	  for ( auto el : get_seqs() ) {
		  if ( el.mate == 0 ) {
			  pos+=el.length;
			  if ( unformatted_barcode.length() >= pos )
				  unformatted_barcode.insert(pos++, "-");
			  else
				  break;
		  }
	  }

	  return unformatted_barcode.substr(0,pos-1);

 }

 std::string get_barcodeString(CountType index) {

	  // invalid index
	  if ( index >= get_barcodeVector().size() ) {
		  return "";
	  }

	  else {

		  std::vector<std::string> bc_vec = get_barcodeVector()[index];

		  std::stringstream ss;
		  for ( auto fragment : bc_vec ) {
			  ss << fragment;
		  }
		  std::string barcode_string = ss.str();
		  return format_barcode(barcode_string);
	  }
 }

 void set_read_structure ( std::vector<std::string> read_argument ) {

	 // Init variables
	 CountType lenSum = 0;
	 CountType length = 0;
	 std::string length_string = "";
	 char type;
	 unsigned mates = 0;
	 std::vector<SequenceElement> temp;

	 // Iterate through input vector
	 for ( auto read = read_argument.begin(); read != read_argument.end(); ++read ) {

		 // Split string into fragment length and fragment type (B or R)
		 length_string = (*read).substr(0,(*read).length()-1);
		 type = (*(*read).rbegin());

		 if ( length_string.find_first_not_of("0123456789")!=std::string::npos ) {
			 throw std::runtime_error("Invalid length for read fragment " + *read  + ". Please only use unsigned integer values.");
		 }

		 try{
			 length = CountType(std::atol(length_string.c_str()));
		 } catch( std::bad_cast & ex ){
			 std::cerr << "Error while casting length " << length_string << " to type uint16_t." << std::endl;
			 throw ex;
		 }

		 if ( type!='B' && type!='R' ) {
			 std::stringstream ss;
			 ss << "\'" << type << "\' is no valid read type. Please use \'R\' for sequencing reads or \'B\' for barcode reads.";
			 throw std::runtime_error(ss.str());

		 }

		 temp.push_back(SequenceElement(temp.size(), (type == 'R') ? ++mates : 0, length));
		 lenSum += length;

	 }
	 set_seqs(temp);
	 set_mates(mates);

	 if ( lenSum!=get_cycles() ) {
		 throw std::runtime_error("Sum of defined reads does not equal the given number of cycles.");
	 }

 }

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

  void disable_ohw(bool value) {
	  set_unmodifiable(discard_ohw, !value, "discard_ohw");
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

  bool get_any_best_hit_mode() {
      return (get_unmodifiable(mode, "mode")==AlignmentMode::ANYBEST);
  }

  bool get_all_hit_mode() {
      return (get_mode()==AlignmentMode::ALL);
  }

  bool get_all_best_hit_mode() {
      return (get_mode()==AlignmentMode::ALLBEST);
  }

  bool get_all_best_n_scores_mode() {
      return (get_mode()==AlignmentMode::BESTN);
  }

  bool get_unique_hit_mode() {
      return (get_mode()==AlignmentMode::UNIQUE);
  }

  CountType get_best_n() {
      return get_unmodifiable(best_n, "best_n", true);
  }

  void set_temp_dir(std::string value) {
      set_unmodifiable(temp_dir, value, "temp_dir");
  }

  std::string get_temp_dir() {
	  std::string dir = get_unmodifiable(temp_dir, "temp_dir");
	  dir = dir!="" ? dir : get_unmodifiable(root, "root");
	  return dir;
  }

  void set_write_bam(bool value) {
	  set_unmodifiable(write_bam, value, "write_bam");
  }

  bool get_write_bam() {
      return get_unmodifiable(write_bam, "write_bam");
  }

  void set_output_cycles(std::vector<uint16_t> cycles) {
	  std::vector<uint16_t> the_cycles;
	  for ( auto it = cycles.begin(); it != cycles.end(); ++it ) {
		  if ( *it > get_cycles() )
			  the_cycles.push_back(get_cycles());
		  else
			  the_cycles.push_back(*it);
	  }
	  std::sort(the_cycles.begin(), the_cycles.end());
	  the_cycles.erase( std::unique(the_cycles.begin(), the_cycles.end()), the_cycles.end());
	  set_unmodifiable(output_cycles, the_cycles, "output_cycles");
  }

  std::vector<uint16_t> get_output_cycles() {
	  return get_unmodifiable(output_cycles, "output_cycles", true);
  }

  bool is_output_cycle(CountType cycle) {
	  auto out_cycles = get_output_cycles();
	  if ( std::find(out_cycles.begin(), out_cycles.end(), cycle) == out_cycles.end() )
		  return false;
	  return true;
  }

  void set_keep_aln_files(std::vector<uint16_t> value) {
	  set_unmodifiable(keep_aln_files, value, "keep_aln_files");
  }

  std::vector<uint16_t> get_keep_aln_files() {
      return get_unmodifiable(keep_aln_files, "keep_aln_files");
  }

  bool is_keep_aln_files_cycle(CountType cycle) {
	  auto aln_files_cycles = get_keep_aln_files();
	  if ( std::find(aln_files_cycles.begin(), aln_files_cycles.end(), cycle) == aln_files_cycles.end() )
		  return false;
	  return true;
  }

  void set_block_size(std::string value) {

	  uint64_t size;
	  char type = 'B';

	  // Split value to size and type
	  if ( value.find_first_of("BKM") != std::string::npos ) {
		  type = *value.rbegin();
		  value = value.substr(0,value.length()-1);
	  }

	  if ( value.find_first_not_of("0123456789")!=std::string::npos ) {
		  throw std::runtime_error("Invalid block size " + value + ". Please only use unsigned integer values.");
	  }

	  try{
		  size = uint64_t(std::atol(value.c_str()));
	  } catch( std::bad_cast & ex ){
		  std::cerr << "Error while casting length " << value << " to type uint16_t." << std::endl;
		  throw ex;
	  }

	  if ( type == 'B' )
		  set_block_size(size);
	  else if ( type == 'K' )
		  set_block_size(size*1024);
	  else if ( type == 'M' )
		  set_block_size(size*1024*1024);
	  else
		  throw std::runtime_error("Invalid block size type. Only 'B' (Bytes), 'K' (Kilobytes) or 'M' (Megabytes) are permitted.");

  }

  uint64_t get_block_size() {
      return get_unmodifiable(block_size, "block_size");
  }

  void set_compression_format(uint16_t value) {
	  if ( value > 2 )
		  value = 2;
	  uint8_t one_byte_value = value;
	  set_unmodifiable(compression_format, one_byte_value, "compression_format");
  }

  uint8_t get_compression_format() {
      return get_unmodifiable(compression_format, "compression_format");
  }

  void set_lanes(std::vector<uint16_t> value) {
	  std::sort( value.begin(), value.end() );
	  value.erase( std::unique( value.begin(), value.end() ), value.end() );
	  set_unmodifiable(lanes, value, "lanes");
  }

  std::vector<uint16_t> get_lanes() {
      return get_unmodifiable(lanes, "lanes", true);
  }

  void set_mode(std::string value) {

	  // All hit mode
	  if ( value == "ALL" || value == "A" ) {
		  set_mode(AlignmentMode::ALL);
	  }

	  // Unique mode
	  if ( value == "UNIQUE" || value == "U" ) {
		  set_mode(AlignmentMode::UNIQUE);
	  }

	  // Best N scores mode
	  else if ( value.substr(0,5) == "BESTN" || value.substr(0,1) == "N" ) {

		  std::string bestn = value.substr(0,5) == "BESTN" ? value.substr(5) : value.substr(1);

		  if ( bestn.find_first_not_of("0123456789")!=std::string::npos ) {
			  throw std::runtime_error("Invalid alignment mode: " + value + ".");
		  }
		  try{
			  set_mode(AlignmentMode::BESTN, CountType(std::atol(bestn.c_str())));
		  } catch( std::bad_cast & ex ){
			  std::cerr << "Error while casting length " << bestn << " to type uint16_t." << std::endl;
			  throw ex;
		  }
	  }

	  // All best mode
	  else if ( value == "ALLBEST" || value == "H" ) {
		  set_mode(AlignmentMode::ALLBEST);
	  }

	  // All hit mode
	  else if ( value == "ANYBEST" || value == "B" ) {
		  set_mode(AlignmentMode::ANYBEST);
	  }

	  // Unknown mode
	  else {
		  throw std::runtime_error("Invalid alignment mode: " + value + ".");
	  }
  }

  void set_tiles(std::vector<uint16_t> value) {
	  std::sort( value.begin(), value.end() );
	  value.erase( std::unique( value.begin(), value.end() ), value.end() );
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

  void set_start_cycle(CountType value) {
	  set_unmodifiable(start_cycle, value, "start_cycle");
  }

  CountType get_start_cycle() {
      CountType ret = get_unmodifiable(start_cycle, "start_cycle", true);

      // Value if not set.
      if ( ret == 0 )
    	  return 1;

      return ret;
  }

  void set_out_dir(std::string value) {
	  set_unmodifiable(out_dir, value, "out_dir");
  }

  std::string get_out_dir() {
      return get_unmodifiable(out_dir, "out_dir");
  }

  void set_num_threads(CountType value) {
	  set_unmodifiable(num_threads, value, "num_threads");
  }

  CountType get_num_threads() {
      return get_unmodifiable(num_threads, "num_threads");
  }

  void set_num_out_threads(CountType value) {
	  set_unmodifiable(num_out_threads, value, "num_out_threads");
  }

  CountType get_num_out_threads() {
      return get_unmodifiable(num_out_threads, "num_out_threads");
  }

  std::vector<SequenceElement> get_seqs() {
      return get_unmodifiable(seqs, "seqs", true);
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

  void set_extended_cigar(bool value) {
	  set_unmodifiable(extended_cigar, value, "extended_cigar");
  }

  bool get_extended_cigar() {
      return get_unmodifiable(extended_cigar, "extended_cigar");
  }

  CountType get_max_consecutive_gaps() {
      return get_unmodifiable(max_consecutive_gaps, "max_consecutive_gaps");
  }

  float get_min_as_ratio() {
        return get_unmodifiable(min_as_ratio, "min_as_ratio");
    }

  void set_min_as_ratio(float value) {
	  if ( value > 1.0f )
		  value = 1.0f;
	  if ( value < 0.0f )
		  value = 0.0f;
	  set_unmodifiable(min_as_ratio, value, "min_as_ratio");
  }

};

#endif
