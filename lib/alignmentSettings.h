#ifndef ALIGNMENTSETTINGS_H
#define ALIGNMENTSETTINGS_H

#include "headers.h"
#include "definitions.h"
#include "tools_static.h"

// Data structure to store the alignment settings
class AlignmentSettings {

private:

  // PARAMETER: Base Call quality cutoff, treat BC with quality < bc_cutoff as miscall
  Unmodifiable<CountType> min_qual;

  // PARAMETER: All-Best-N-Scores-Mode::N
  Unmodifiable<CountType> best_n;

  // PARAMETER: temporary directory for the streamed alignment
  Unmodifiable<std::string> temp_dir;

  // SWITCH: write sam/bam output or not
  Unmodifiable<bool> write_bam;

  // PARAMETER: Cycles for intermediate SAM/BAM output
  Unmodifiable<std::vector<uint16_t>> output_cycles;

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

  //PARAMETER: Stores the barcodes defined by the user. The inner vector contains the single fragments of multi-barcodes.
  Unmodifiable<std::vector<std::vector<std::string>>> barcodeVector;

  // PARAMETER: directory in which to create the output directory structure 
  Unmodifiable<std::string> out_dir;

  // PARAMETER: number of threads to use
  Unmodifiable<CountType> num_threads;

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

  // Definitions for FM index.
  Unmodifiable<CountType> anchorLength;
  Unmodifiable<CountType> errorRate;

  // Scoring scheme
  Unmodifiable<ScoreType> min_as;
  Unmodifiable<CountType> match_score;
  Unmodifiable<CountType> mismatch_penalty;
  Unmodifiable<CountType> insertion_opening_penalty;
  Unmodifiable<CountType> deletion_opening_penalty;
  Unmodifiable<CountType> insertion_extension_penalty;
  Unmodifiable<CountType> deletion_extension_penalty;
  Unmodifiable<float> softclip_opening_penalty;
  Unmodifiable<float> softclip_extension_penalty;
  Unmodifiable<CountType> max_gap_length;

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
	  xml_out.add_child("settings.out.min_as", getXMLnode ( get_min_as()) );

	  // Technical settings
	  xml_out.add_child("settings.technical.num_threads", getXMLnode ( get_num_threads() ));
	  xml_out.add_child("settings.technical.keep_aln_files", getXMLnode ( get_keep_aln_files() ));
	  xml_out.add_child("settings.technical.block_size", getXMLnode ( get_block_size() ));
	  xml_out.add_child("settings.technical.compression_format", getXMLnode ( get_compression_format() ));

	  // Alignment algorithm settings
	  xml_out.add_child("settings.align.min_qual", getXMLnode (get_min_qual() ));
	  xml_out.add_child("settings.align.anchor", getXMLnode ( get_anchor_length() ));

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

  void set_min_qual(CountType value) {
    	  set_unmodifiable(min_qual, value, "min_qual");
  }

  CountType get_min_qual() {
      return get_unmodifiable(min_qual, "min_qual");
  }

  void set_min_as(ScoreType value) {
	  set_unmodifiable(min_as, value, "min_as");
  }

  ScoreType get_min_as() {
	  return get_unmodifiable(min_as, "min_as");
  }

  void set_max_gap_length(CountType value) {
	  set_unmodifiable(max_gap_length, value, "max_gap_length");
  }

  CountType get_max_gap_length() {
	  return get_unmodifiable(max_gap_length, "max_gap_length");
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

  CountType get_best_n() {
      return get_unmodifiable(best_n, "best_n", true);
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

  void set_keep_aln_files(bool value) {
	  set_unmodifiable(keep_aln_files, value, "keep_aln_files");
  }

  bool get_keep_aln_files() {
      return get_unmodifiable(keep_aln_files, "keep_aln_files");
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

  CountType get_anchor_length() {
	  return get_unmodifiable(anchorLength, "anchorLength");
  }

  void set_anchor_length( CountType value ) {
	  set_unmodifiable(anchorLength, value, "anchorLength");
	  set_unmodifiable(errorRate, CountType(value/2), "errorRate");
  }

  CountType get_error_rate() {
	  return get_unmodifiable(errorRate, "errorRate");
  }

  CountType get_match_score() {
	  return get_unmodifiable(match_score, "match_score");
  }

  void set_match_score(CountType value) {
	  set_unmodifiable(match_score, value, "match_score");
  }

  CountType get_mismatch_penalty() {
	  return get_unmodifiable(mismatch_penalty, "mismatch_penalty");
  }

  void set_mismatch_penalty(CountType value) {
	  set_unmodifiable(mismatch_penalty, value, "mismatch_penalty");
  }

  CountType get_insertion_opening_penalty() {
	  return get_unmodifiable(insertion_opening_penalty, "insertion_opening_penalty");
  }

  void set_insertion_opening_penalty(CountType value) {
	  set_unmodifiable(insertion_opening_penalty, value, "insertion_opening_penalty");
  }

  CountType get_deletion_opening_penalty() {
	  return get_unmodifiable(deletion_opening_penalty, "deletion_opening_penalty");
  }

  void set_deletion_opening_penalty(CountType value) {
	  set_unmodifiable(deletion_opening_penalty, value, "deletion_opening_penalty");
  }

  CountType get_insertion_extension_penalty() {
	  return get_unmodifiable(insertion_extension_penalty, "insertion_extension_penalty");
  }

  void set_insertion_extension_penalty(CountType value) {
	  set_unmodifiable(insertion_extension_penalty, value, "insertion_extension_penalty");
  }

  CountType get_deletion_extension_penalty() {
	  return get_unmodifiable(deletion_extension_penalty, "deletion_extension_penalty");
  }

  void set_deletion_extension_penalty(CountType value) {
	  set_unmodifiable(deletion_extension_penalty, value, "deletion_extension_penalty");
  }

  float get_softclip_extension_penalty() {
	  return get_unmodifiable(softclip_extension_penalty, "softclip_extension_penalty");
  }

  void set_softclip_extension_penalty(float value) {
	  set_unmodifiable(softclip_extension_penalty, value, "softclip_extension_penalty");
  }

  float get_softclip_opening_penalty() {
	  return get_unmodifiable(softclip_opening_penalty, "softclip_opening_penalty");
  }

  void set_softclip_opening_penalty(float value) {
	  set_unmodifiable(softclip_opening_penalty, value, "softclip_opening_penalty");
  }

};

#endif
