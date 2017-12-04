#ifndef ALIGNMENTSETTINGS_H
#define ALIGNMENTSETTINGS_H

#include "headers.h"
#include "definitions.h"
#include "tools_static.h"

/**
 * Data structure to store the alignment settings.
 * @author Martin Lindner
 */
class AlignmentSettings {

private:

///////////////////////////////////////
////////// Technical options //////////
///////////////////////////////////////

    /** Memory block size for the input and output buffer in the streamed alignment. */
    Unmodifiable<uint64_t> block_size;

    /** Compression format for alignment files. */
    Unmodifiable<uint8_t> compression_format;

    /** Number of threads to use. */
    Unmodifiable<CountType> num_threads;

    // SWITCH: Keep the old alignment files of previous cycles
    Unmodifiable<std::vector<uint16_t>> keep_aln_files;


////////////////////////////////////
////////// Output options //////////
////////////////////////////////////

	/** Value N for best_n output mode. */
    Unmodifiable<CountType> best_n;

    /** Switch: Write BAM instead of SAM output. */
    Unmodifiable<bool> write_bam;

    /** Cycles for (intermediate) SAM/BAM output. */
    Unmodifiable<std::vector<uint16_t>> output_cycles;

    /** Switch: Activate extended CIGAR annotation. */
    Unmodifiable<bool> extended_cigar;

    /** Output mode. */
    Unmodifiable<AlignmentMode> mode;

    /** Switch: store all sequences. */
    Unmodifiable<bool> keep_all_sequences;


////////////////////////////////////
////////// Scoring scheme //////////
////////////////////////////////////

    /** Minimal alignment score to keep an alignment. */
    Unmodifiable<ScoreType> min_as;

    /** Score for an alignment match. */
    Unmodifiable<CountType> match_score;

    /** Penalty for an alignment mismatch. */
    Unmodifiable<CountType> mismatch_penalty;

    /** Penalty for opening an insertion region. */
    Unmodifiable<CountType> insertion_opening_penalty;

    /** Penalty for opening a deletion region. */
    Unmodifiable<CountType> deletion_opening_penalty;

    /** Penalty for extending an insertion region. */
    Unmodifiable<CountType> insertion_extension_penalty;

    /** Penalty for extending a deletion region. */
    Unmodifiable<CountType> deletion_extension_penalty;

    /** Penalty for opening a softclip region (only for output, not for filtering). */
    Unmodifiable<float> softclip_opening_penalty;

    /** Penalty for extending a softclip region (only for output, not for filtering). */
    Unmodifiable<float> softclip_extension_penalty;

    /** Maximal consecutive number of insertions or deletions. */
    Unmodifiable<CountType> max_gap_length;

    // PARAMETER: max. amount of threads used for output
    Unmodifiable<CountType> num_out_threads;

	// PARAMETER: the first cycle to handle. Should be 1 by default.
	Unmodifiable<CountType> start_cycle;

    /** Maximal relative length of softclip region. */
    Unmodifiable<float> max_softclip_ratio;


///////////////////////////////////////
////////// Alignment options //////////
///////////////////////////////////////

    /** Length of the alignment anchor. */
    Unmodifiable<CountType> anchorLength;

    /** Interval to increase the number of tolerated errors. */
    Unmodifiable<CountType> errorRate;

    /** Interval to create new seeds. */
    Unmodifiable<CountType> seeding_interval;

	/** Minimal quality for a base call to be valid. */
	Unmodifiable<CountType> min_qual;
	Unmodifiable<bool> force_resort;

    /** List of lanes to process. */
    Unmodifiable<std::vector<uint16_t>> lanes;

    /** List of tiles to process. */
    Unmodifiable<std::vector<uint16_t>> tiles;

    /** Total length of the sequencing procedure. */
    Unmodifiable<CountType> cycles;

    /** List of all sequence fragments. */
    Unmodifiable<std::vector<SequenceElement>> seqs;

    /** Total number of mates. */
    Unmodifiable<uint16_t> mates;


////////////////////////////////////////////
////////// Demultiplexing options //////////
////////////////////////////////////////////

    /** List of barcode sequences (split in an inner vector for multibarcodes). */
    Unmodifiable<std::vector<std::vector<std::string>>> barcodeVector;

    /** Number of tolerated mismatches for barcodes. */
    Unmodifiable<std::vector<uint16_t>> barcode_errors;

    /** Switch: Keep all barcodes (disables barcode filtering). */
    Unmodifiable<bool> keep_all_barcodes;


//////////////////////////////////////////
////////// File/Directory paths //////////
//////////////////////////////////////////

    /** Temporary directory for the streamed alignment. */
    Unmodifiable<std::string> temp_dir;

    /** Root directory of hilive run. */
    Unmodifiable<std::string> root;

    /** Index file(s) prefix. */
    Unmodifiable<std::string> index_fname;

    /** Output directory. */
    Unmodifiable<std::string> out_dir;


///////////////////////////////////////
////////// Private functions //////////
///////////////////////////////////////

    /**
     * Parse the barcode vector to a string vector that can be used for XML output.
     * @return Vector with one string for each barcode. Subsequences are delimited with a "-" character.
     * @author Tobias Loka
     */
    std::vector<std::string> xmlParse_barcodeVector() {
    	std::vector<std::string> bc_strings;
    	for ( CountType i = 0; i < get_barcodeVector().size(); i++ ) {
    		bc_strings.push_back( get_barcodeString(i) );
    	}
    	return bc_strings;
    }

    /**
     * Parse the seq vector to a string vector that can be used for XML output.
     * @return Vector with one string for each element in format [0-9]*[BR]
     * @author Tobias Loka
     */
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

    /**
     * Set the barcode vector. Can only be set once (unmodifiable object).
     * @param value The barcode split in a 2D-vector, where first dimension are the different barcodes and second dimension the fragments of multi-barcodes.
     * @author Tobias Loka
     */
    void set_barcodeVector(std::vector<std::vector<std::string> > value) {
    	set_unmodifiable(barcodeVector, value, "barcodeVector");
    }

    /**
     * Set the seq vector. Can only be set once (unmodifiable object).
     * @param value Vector of SequenceElement objects that contain all information about the read structure.
     * @author Tobias Loka
     */
    void set_seqs(std::vector<SequenceElement> value) {
    	set_unmodifiable(seqs, value, "seqs");
    }

    /**
     * Set the block size. Can only be set once (unmodifiable object).
     * @param value The block size as unsigned integer.
     * @author Tobias Loka
     */
    void set_block_size(uint64_t value) {
    	set_unmodifiable(block_size, value, "block_size");
    }

    /**
     * Set the number of mates. Can only be set once (unmodifiable object).
     * @param value The number of mates as unsigned integer.
     * @author Tobias Loka
     */
    void set_mates(uint16_t value) {
    	set_unmodifiable(mates, value, "mates");
    }

    /**
     * Get the output mode.
     * @return The output mode as AlignmentMode object.
     * @author Tobias Loka
     */
    AlignmentMode get_mode() {
    	return get_unmodifiable(mode, "mode");
    }

    /**
     * Set the output mode. Can only be set once (unmodifiable object).
     * @param value The output mode as AlignmentMode object.
     * @param bestn (optional) Defines the N for the best_n mode (0 if not defined).
     * @author Tobias Loka
     */
    void set_mode(AlignmentMode value, CountType bestn = 0) {
    	set_unmodifiable(mode, value, "mode");
    	set_unmodifiable(best_n, bestn, "best_n");
    }


public:

///////////////////////////////////////////
////////// General functionality //////////
///////////////////////////////////////////

    /**
     * Create a property tree representing the AlignmentSettings object that is filled with all (relevant) settings.
     * @return Property tree containing all relevant settings.
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
    	xml_out.add_child("settings.out.keep_all_sequences", getXMLnode ( get_keep_all_sequences() ));

    	// Technical settings
    	xml_out.add_child("settings.technical.num_threads", getXMLnode ( get_num_threads() ));
    	xml_out.add_child("settings.technical.num_out_threads", getXMLnode ( get_num_out_threads() ));
    	xml_out.add_child("settings.technical.keep_aln_files", getXMLnode_vector ( get_keep_aln_files() ));
    	xml_out.add_child("settings.technical.block_size", getXMLnode ( get_block_size() ));
    	xml_out.add_child("settings.technical.compression_format", getXMLnode ( get_compression_format() ));

    	// Scoring scheme
    	xml_out.add_child("settings.scores.match_score", getXMLnode ( get_match_score() ));
    	xml_out.add_child("settings.scores.mismatch_penalty", getXMLnode ( get_mismatch_penalty() ));
    	xml_out.add_child("settings.scores.insertion_opening_penalty", getXMLnode ( get_insertion_opening_penalty() ));
    	xml_out.add_child("settings.scores.insertion_extension_penalty", getXMLnode ( get_insertion_extension_penalty() ));
    	xml_out.add_child("settings.scores.deletion_opening_penalty", getXMLnode ( get_deletion_opening_penalty() ));
    	xml_out.add_child("settings.scores.deletion_extension_penalty", getXMLnode ( get_deletion_extension_penalty() ));
    	xml_out.add_child("settings.scores.max_gap_length", getXMLnode ( get_max_gap_length() ));
    	xml_out.add_child("settings.scores.softclip_opening_penalty", getXMLnode ( get_softclip_opening_penalty() ));
    	xml_out.add_child("settings.scores.softclip_extension_penalty", getXMLnode ( get_softclip_extension_penalty() ));
    	xml_out.add_child("settings.scores.max_softclip_ratio", getXMLnode ( get_max_softclip_ratio() ));

    	// Alignment algorithm settings
    	xml_out.add_child("settings.align.min_qual", getXMLnode (get_min_qual() ));
    	xml_out.add_child("settings.align.anchor", getXMLnode ( get_anchor_length() ));
    	xml_out.add_child("settings.align.error_interval", getXMLnode ( get_error_rate() ));
    	xml_out.add_child("settings.align.seeding_interval", getXMLnode ( get_seeding_interval() ));

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

     		// Check correct number of fragments
     		if ( barcode_lengths.size() != fragments.size()) {
   			throw std::runtime_error("Wrong number of fragments for barcode " + *barcode + " (should have " + std::to_string(barcode_lengths.size()) + " fragments).");
     		}

     		// Check correct length of all fragments
     		for ( uint16_t num = 0; num != fragments.size(); num++ ) {
     			if ( fragments[num].length() != barcode_lengths[num] ) {
     				throw std::runtime_error("Wrong fragment length in barcode " + *barcode);
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
   	  else if ( value == "UNIQUE" || value == "U" ) {
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

    bool get_keep_all_sequences() {
    	return get_unmodifiable(keep_all_sequences, "keep-all-sequences");
    }

    void set_keep_all_sequences(bool value) {
    	set_unmodifiable(keep_all_sequences, value, "kepp-all-sequences");
    }



////////// Scoring scheme  //////////

    /**
     * Get the minimum alignment score for an alignment to be considered.
     * @return Minimum score for an alignment.
     * @author Tobias Loka
     */
    ScoreType get_min_as() {
    	return get_unmodifiable(min_as, "min_as");
    }

    /**
     * Set the minimum alignment score for an alignment to be considered.
     * @param value Minimum score for an alignment.
     * @author Tobias Loka
     */
    void set_min_as(ScoreType value) {
    	set_unmodifiable(min_as, value, "min_as");
    }

    /**
     * Get the maximumlength of consecutive gaps (insertions or deletions).
     * @return Maximum gap length.
     * @author Tobias Loka
     */
    CountType get_max_gap_length() {
    	return get_unmodifiable(max_gap_length, "max_gap_length");
    }

    /**
     * Set the maximumlength of consecutive gaps (insertions or deletions).
     * @param value Maximum gap length.
     * @author Tobias Loka
     */
    void set_max_gap_length(CountType value) {
    	set_unmodifiable(max_gap_length, value, "max_gap_length");
    }

    /**
     * Get the score for a match.
     * @return Score for a match.
     * @author Tobias Loka
     */
    CountType get_match_score() {
  	  return get_unmodifiable(match_score, "match_score");
    }

    /**
     * Set the score for a match.
     * @param value Score for a match.
     * @author Tobias Loka
     */
    void set_match_score(CountType value) {
  	  set_unmodifiable(match_score, value, "match_score");
    }

    /**
     * Get the penalty for a mismatch.
     * @return Penalty for a mismatch.
     * @author Tobias Loka
     */
    CountType get_mismatch_penalty() {
  	  return get_unmodifiable(mismatch_penalty, "mismatch_penalty");
    }

    /**
     * Set the penalty for a mismatch.
     * @param value Penalty for a mismatch.
     * @author Tobias Loka
     */
    void set_mismatch_penalty(CountType value) {
  	  set_unmodifiable(mismatch_penalty, value, "mismatch_penalty");
    }

    /**
     * Get the penalty for an insertion opening.
     * @return Penalty for an insertion opening.
     * @author Tobias Loka
     */
    CountType get_insertion_opening_penalty() {
  	  return get_unmodifiable(insertion_opening_penalty, "insertion_opening_penalty");
    }

    /**
     * Set the penalty for an insertion opening.
     * @param value Penalty for an insertion opening.
     * @author Tobias Loka
     */
    void set_insertion_opening_penalty(CountType value) {
  	  set_unmodifiable(insertion_opening_penalty, value, "insertion_opening_penalty");
    }

    /**
     * Get the penalty for a deletion opening.
     * @return Penalty for a deletion opening.
     * @author Tobias Loka
     */
    CountType get_deletion_opening_penalty() {
  	  return get_unmodifiable(deletion_opening_penalty, "deletion_opening_penalty");
    }

    /**
     * Set the penalty for a deletion opening.
     * @param value Penalty for a deletion opening.
     * @author Tobias Loka
     */
    void set_deletion_opening_penalty(CountType value) {
  	  set_unmodifiable(deletion_opening_penalty, value, "deletion_opening_penalty");
    }

    /**
     * Get the penalty for an insertion extension.
     * @return Penalty for an insertion extension.
     * @author Tobias Loka
     */
    CountType get_insertion_extension_penalty() {
  	  return get_unmodifiable(insertion_extension_penalty, "insertion_extension_penalty");
    }

    /**
     * Set the penalty for an insertion extension.
     * @param value Penalty for an insertion extension.
     * @author Tobias Loka
     */
    void set_insertion_extension_penalty(CountType value) {
  	  set_unmodifiable(insertion_extension_penalty, value, "insertion_extension_penalty");
    }

    /**
     * Get the penalty for a deletion extension.
     * @return Penalty for a deletion extension.
     * @author Tobias Loka
     */
    CountType get_deletion_extension_penalty() {
  	  return get_unmodifiable(deletion_extension_penalty, "deletion_extension_penalty");
    }

    /**
     * Set the penalty for a deletion extension.
     * @param value Penalty for a deletion extension.
     * @author Tobias Loka
     */
    void set_deletion_extension_penalty(CountType value) {
  	  set_unmodifiable(deletion_extension_penalty, value, "deletion_extension_penalty");
    }

    /**
     * Get the penalty for a softclip extension. This is only used for the final output, not for filtering during the alignment algorithm.
     * @return Penalty for a softclip extension.
     * @author Tobias Loka
     */
    float get_softclip_extension_penalty() {
  	  return get_unmodifiable(softclip_extension_penalty, "softclip_extension_penalty");
    }

    /**
     * Set the penalty for a softclip extension.
     * @param value Penalty for a softclip extension.
     * @author Tobias Loka
     */
    void set_softclip_extension_penalty(float value) {
  	  set_unmodifiable(softclip_extension_penalty, value, "softclip_extension_penalty");
    }

    /**
     * Get the penalty for a softclip opening. This is only used for the final output, not for filtering during the alignment algorithm.
     * @return Penalty for a softclip opening.
     * @author Tobias Loka
     */
    float get_softclip_opening_penalty() {
  	  return get_unmodifiable(softclip_opening_penalty, "softclip_opening_penalty");
    }

    /**
     * Set the penalty for a softclip opening.
     * @param value Penalty for a softclip opening.
     * @author Tobias Loka
     */
    void set_softclip_opening_penalty(float value) {
  	  set_unmodifiable(softclip_opening_penalty, value, "softclip_opening_penalty");
    }

    /**
     * Get the maximal relative length of the front softclip.
     * @return Maximal softclip ratio.
     * @author Tobias Loka
     */
    float get_max_softclip_ratio() {
    	return get_unmodifiable(max_softclip_ratio, "max_softclip_ratio");
    }

    /**
     * Set the maximal relative length of the front softclip.
     * @param value The maximal length of the front softclip.
     * @author Tobias Loka
     */
    void set_max_softclip_ratio(float value) {
    	set_unmodifiable(max_softclip_ratio, value, "max_softclip_ratio");
    }



    /**
     * Get the anchor length.
     * @return The anchor length.
     * @author Tobias Loka
     */
    CountType get_anchor_length() {
  	  return get_unmodifiable(anchorLength, "anchorLength");
    }

    /**
     * Set the anchor length.
     * @param value The anchor length.
     * @author Tobias Loka
     */
    void set_anchor_length( CountType value ) {
  	  set_unmodifiable(anchorLength, value, "anchorLength");
    }

    /**
     * Get the frequency to increase the tolerated number of errors during the alingment algorithm.
     * @return Frequency to tolerate more errors during the alignment algorithm.
     * @author Tobias Loka
     */
    CountType get_error_rate() {
  	  return get_unmodifiable(errorRate, "errorRate");
    }

    /**
     * Set the frequency to increase the tolerated number of errors during the alingment algorithm.
     * @param value Frequency to tolerate more errors during the alignment algorithm.
     * @author Tobias Loka
     */
    void set_error_rate(CountType value) {
    	set_unmodifiable(errorRate, value, "errorRate");
    }

    /**
     * Get the interval to create new seeds.
     * @return The interval to create new seeds.
     * @author Tobias Loka
     */
    CountType get_seeding_interval() {
    	return get_unmodifiable(seeding_interval, "seeding_interval");
    }

    /**
     * Get the interval to create new seeds.
     * @param value The interval to create new seeds.
     * @author Tobias Loka
     */
    void set_seeding_interval(CountType value) {

    	// interval must be at least 1
    	value = value == 0 ? 1 : value;

    	set_unmodifiable(seeding_interval, value, "seeding_interval");
    }


  bool get_force_resort() {
	  return get_unmodifiable(force_resort, "force_resort");
  }

  void set_force_resort(bool value) {
	  set_unmodifiable(force_resort, value, "force_resort");
  }


};

#endif
