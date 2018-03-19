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

	//////////////////////////////
	////////// Run Info //////////
	//////////////////////////////

	/** Unique ID of the run. */
	Immutable<std::string> run_id;

	/** Unique ID of the flowcell. */
	Immutable<std::string> flowcell_id;

	/** Unique ID of the instrument. */
	Immutable<std::string> instrument_id;

	///////////////////////////////////////
	////////// Technical options //////////
	///////////////////////////////////////

	/** Memory block size for the input and output buffer in the streamed alignment. */
	Immutable<uint64_t> block_size;

	/** Compression format for alignment files. */
	Immutable<uint8_t> compression_format;

	/** Number of threads to use. */
	Immutable<CountType> num_threads;

	// SWITCH: Keep the old alignment files of previous cycles
	Immutable<std::vector<uint16_t>> keep_aln_files;


	////////////////////////////////////
	////////// Output options //////////
	////////////////////////////////////

	/** Value N for best_n output mode. */
	Immutable<CountType> best_n;

	/** Switch: Write BAM instead of SAM output. */
	Immutable<OutputFormat> output_format;

	/** Cycles for (intermediate) SAM/BAM output. */
	Immutable<std::vector<uint16_t>> output_cycles;

	/** Switch: Activate extended CIGAR annotation. */
	Immutable<bool> extended_cigar;

	/** Output mode. */
	Immutable<AlignmentMode> mode;

	/** Switch: store all sequences. */
	Immutable<bool> keep_all_sequences;

	/** Switch: Report unmapped reads. */
	Immutable<bool> report_unmapped;


	////////////////////////////////////
	////////// Scoring scheme //////////
	////////////////////////////////////

	/** Minimal alignment score to keep an alignment. */
	Immutable<ScoreType> min_as;

	/** Score for an alignment match. */
	Immutable<CountType> match_score;

	/** Penalty for an alignment mismatch. */
	Immutable<CountType> mismatch_penalty;

	/** Penalty for opening an insertion region. */
	Immutable<CountType> insertion_opening_penalty;

	/** Penalty for opening a deletion region. */
	Immutable<CountType> deletion_opening_penalty;

	/** Penalty for extending an insertion region. */
	Immutable<CountType> insertion_extension_penalty;

	/** Penalty for extending a deletion region. */
	Immutable<CountType> deletion_extension_penalty;

	/** Penalty for opening a softclip region (only for output, not for filtering). */
	Immutable<float> softclip_opening_penalty;

	/** Penalty for extending a softclip region (only for output, not for filtering). */
	Immutable<float> softclip_extension_penalty;

	/** Maximal consecutive number of insertions or deletions. */
	Immutable<CountType> max_gap_length;

	// PARAMETER: max. amount of threads used for output
	Immutable<CountType> num_out_threads;

	// PARAMETER: the first cycle to handle. Should be 1 by default.
	Immutable<CountType> start_cycle;

	/** Maximal relative length of softclip region. */
	Immutable<float> max_softclip_ratio;


	///////////////////////////////////////
	////////// Alignment options //////////
	///////////////////////////////////////

	/** Length of the alignment anchor. */
	Immutable<CountType> anchor_length;

	/** Interval to increase the number of tolerated errors. */
	Immutable<CountType> error_rate;

	/** Interval to create new seeds. */
	Immutable<CountType> seeding_interval;

	/** Minimal quality for a base call to be valid. */
	Immutable<CountType> min_qual;

	/** Sort the alignment files for output even if sorted files already exist.*/
	Immutable<bool> force_resort;

	/** List of lanes to process. */
	Immutable<std::vector<uint16_t>> lanes;

	/** List of tiles to process. */
	Immutable<std::vector<uint16_t>> tiles;

	/** Total length of the sequencing procedure. */
	Immutable<CountType> cycles;

	/** List of all sequence fragments. */
	Immutable<std::vector<SequenceElement>> seqs;

	/** Total number of mates. */
	Immutable<uint16_t> mates;


	////////////////////////////////////////////
	////////// Demultiplexing options //////////
	////////////////////////////////////////////

	/** List of barcode sequences (split in an inner vector for multibarcodes). */
	Immutable<std::vector<std::vector<std::string>>> barcode_vector;

	/** Number of tolerated mismatches for barcodes. */
	Immutable<std::vector<uint16_t>> barcode_errors;

	/** Switch: Keep all barcodes (disables barcode filtering). */
	Immutable<bool> keep_all_barcodes;


	//////////////////////////////////////////
	////////// File/Directory paths //////////
	//////////////////////////////////////////

	/** Temporary directory for the streamed alignment. */
	Immutable<std::string> temp_dir;

	/** Root directory of hilive run. */
	Immutable<std::string> root;

	/** Index file(s) prefix. */
	Immutable<std::string> index_fname;

	/** Output directory. */
	Immutable<std::string> out_dir;


	///////////////////////////////////////
	////////// Private functions //////////
	///////////////////////////////////////

	/**
	 * Parse the seq vector to a string vector that can be used for XML output.
	 * @return Vector with one string for each element in format [0-9]*[BR]
	 * @author Tobias Loka
	 */
	std::vector<std::string> seqs_to_string_vector() {
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
	 * Set the barcode vector. Can only be set once (immutable object).
	 * @param value The barcode split in a 2D-vector, where first dimension are the different barcodes and second dimension the fragments of multi-barcodes.
	 * @author Tobias Loka
	 */
	void set_barcode_vector(std::vector<std::vector<std::string> > value) {
		set_immutable(barcode_vector, value);
	}

	/**
	 * Set the seq vector. Can only be set once (immutable object).
	 * @param value Vector of SequenceElement objects that contain all information about the read structure.
	 * @author Tobias Loka
	 */
	void set_seqs(std::vector<SequenceElement> value) {
		set_immutable(seqs, value);
	}

	/**
	 * Set the block size. Can only be set once (immutable object).
	 * @param value The block size as unsigned integer.
	 * @author Tobias Loka
	 */
	void set_block_size(uint64_t value) {
		set_immutable(block_size, value);
	}

	/**
	 * Set the number of mates. Can only be set once (immutable object).
	 * @param value The number of mates as unsigned integer.
	 * @author Tobias Loka
	 */
	void set_mates(uint16_t value) {
		set_immutable(mates, value);
	}

	/**
	 * Set the output mode. Can only be set once (immutable object).
	 * @param value The output mode as AlignmentMode object.
	 * @param bestn (optional) Defines the N for the best_n mode (0 if not defined).
	 * @author Tobias Loka
	 */
	void set_mode(AlignmentMode value, CountType bestn = 0) {
		set_immutable(mode, value);
		set_immutable(best_n, bestn);
	}


public:

	/**
	 * Format a barcode from the internal structure to output style.
	 * @param unformatted_barcode A barcode as single string without delimiter
	 * @return Barcode with a "-" char between dual barcodes. Lengths of the single barcodes is retrieved from the internal data structure.
	 */
	std::string format_barcode(std::string unformatted_barcode) const {
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


	/**
	 * Create a property tree representing the AlignmentSettings object that is filled with all (relevant) settings.
	 * @return Property tree containing all relevant settings.
	 * @author Tobias Loka
	 */
	boost::property_tree::ptree to_ptree() {

		boost::property_tree::ptree ptree;

		// Sequencing options
		putConfigNode ( ptree, "bcl-dir", get_root() );
		putConfigNode ( ptree, "lanes", join( get_lanes() ) );
		putConfigNode ( ptree, "tiles", join( get_tiles() ) );
		putConfigNode ( ptree, "reads", join( seqs_to_string_vector() ) );
		putConfigNode ( ptree, "barcodes", join( get_barcode_string_vector() ) );
		putConfigNode ( ptree, "run-id", get_run_id() );
		putConfigNode ( ptree, "flowcell-id", get_flowcell_id() );
		putConfigNode ( ptree, "instrument-id", get_instrument_id() );

		// Report options
		putConfigNode ( ptree, "out-dir", get_out_dir() );
		putConfigNode ( ptree, "out-format", to_string( get_output_format() ) );
		putConfigNode ( ptree, "out-cycles", join( get_output_cycles() ) );
		putConfigNode ( ptree, "out-mode", to_string ( get_mode(), get_best_n() ) );
		putConfigNode ( ptree, "report-unmapped", get_report_unmapped() );
		putConfigNode ( ptree, "extended-cigar", get_extended_cigar() );
		putConfigNode ( ptree, "force-resort", get_force_resort() );
		putConfigNode ( ptree, "max-softclip-ratio", get_max_softclip_ratio() );

		// Alignment options
		putConfigNode( ptree, "index",  get_index_fname() );
		// Align mode not needed since the implied parameters are stored.
		putConfigNode( ptree, "anchor-length",  get_anchor_length() );
		putConfigNode( ptree, "error-interval",  get_error_rate() );
		putConfigNode( ptree, "seeding-interval",  get_seeding_interval() );
		putConfigNode ( ptree, "barcode-errors", join( get_barcode_errors() ) );
		putConfigNode( ptree, "align-undetermined-barcodes",  get_keep_all_barcodes() );
		putConfigNode( ptree, "min-basecall-quality", get_min_qual() );
		putConfigNode( ptree, "keep-invalid-sequences",  get_keep_all_sequences() );

		// Scoring options
		putConfigNode( ptree, "min-as",  get_min_as() );
		putConfigNode( ptree, "match-score",  get_match_score() );
		putConfigNode( ptree, "mismatch-penalty",  get_mismatch_penalty() );
		putConfigNode( ptree, "insertion-opening-penalty",  get_insertion_opening_penalty() );
		putConfigNode( ptree, "insertion-extension-penalty",  get_insertion_extension_penalty() );
		putConfigNode( ptree, "deletion-opening-penalty",  get_deletion_opening_penalty() );
		putConfigNode( ptree, "deletion-extension-penalty",  get_deletion_extension_penalty() );
		putConfigNode( ptree, "max-gap-length",  get_max_gap_length() );
		putConfigNode( ptree, "softclip-opening-penalty",  get_softclip_opening_penalty() );
		putConfigNode( ptree, "softclip-extension-penalty",  get_softclip_extension_penalty() );

		// Technical options
		putConfigNode( ptree, "temp-dir",  get_temp_dir() );
		putConfigNode ( ptree, "keep-files", join( get_keep_aln_files() ) );
		// Keep all files not necessary ==> stored as keep-files
		putConfigNode( ptree, "block-size",  get_block_size() );
		putConfigNode( ptree, "compression",  get_compression_format() );
		putConfigNode( ptree, "num-threads",  get_num_threads() );
		putConfigNode( ptree, "num-out-threads",  get_num_out_threads() );

		return ptree;
	}


	//////////////////////////////////////////////////////////
	////////// Getter / Setter (Alphabetical order) //////////
	//////////////////////////////////////////////////////////

	/**
	 * Get the anchor length.
	 * @return The anchor length.
	 * @author Tobias Loka
	 */
	CountType get_anchor_length() const {
		return get_immutable(anchor_length);
	}

	/**
	 * Set the anchor length.
	 * @param value The anchor length.
	 * @author Tobias Loka
	 */
	void set_anchor_length( CountType value ) {
		set_immutable(anchor_length, value);
	}

	/**
	 * Create the internal barcode structure.
	 * @param barcodeArg A vector of strings as given via the command line (dual barcodes being split by a "-" char).
	 */
	void set_barcodes ( std::vector< std::string > barcodeArg ) {

		// Get the barcode length(s) from the seqs vector
		std::vector<uint16_t> barcode_lengths;
		for ( uint16_t seq_num = 0; seq_num < get_seqs().size(); seq_num++ ) {
			if ( get_seq_by_id(seq_num).isBarcode() )
				barcode_lengths.push_back( get_seq_by_id(seq_num).length );
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
			split(fragments, *barcode, "-");

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

		set_barcode_vector(barcodeVector);

	}

	/**
	 * Create the internal barcode structure.
	 * @param barcodeArg A comma-separated string (dual barcodes being split by a "-" char).
	 */
	void set_barcodes ( std::string barcodeArg ) {
		std::vector<std::string> barcodes;
		split( barcodes, barcodeArg );
		set_barcodes( barcodes );
	}

	/**
	 * Get the number of errors allowed for each barcode.
	 * @param value The number of errors for each barcode as a vector of ints.
	 */
	std::vector<uint16_t> get_barcode_errors() const {
		return get_immutable(barcode_errors);
	}

	/**
	 * Set the number of errors allowed for each barcode.
	 * @param value The number of errors for each barcode as a vector of ints.
	 */
	void set_barcode_errors(std::vector<uint16_t> value) {
		set_immutable(barcode_errors, value);
	}

	void set_barcode_errors(std::string value) {
		std::vector<uint16_t> barcode_errors;
		split( barcode_errors, value );
		set_barcode_errors(barcode_errors);
	}

	/**
	 * Get a barcode from the internal list of barcodes as string.
	 * @param index Index of the requested barcode.
	 * @return The requested barcode as formatted string (including "-" delimiting between single barcodes).
	 */
	std::string get_barcode_string(CountType index) const {

		// invalid index
		if ( index >= get_barcode_vector().size() ) {
			return "";
		}

		else {

			std::vector<std::string> bc_vec = get_barcode_vector()[index];

			std::stringstream ss;
			for ( auto fragment : bc_vec ) {
				ss << fragment;
			}
			std::string barcode_string = ss.str();
			return format_barcode(barcode_string);
		}
	}

	/**
	 * Parse the barcode vector to a string vector that can be used for XML output.
	 * @return Vector with one string for each barcode. Subsequences are delimited with a "-" character.
	 * @author Tobias Loka
	 */
	std::vector<std::string> get_barcode_string_vector() {
		std::vector<std::string> bc_strings;
		for ( CountType i = 0; i < get_barcode_vector().size(); i++ ) {
			bc_strings.push_back( get_barcode_string(i) );
		}
		return bc_strings;
	}

	/**
	 * Get the barcodes as stored internally.
	 * @return Barcodes as a 2-dimensional vector of strings.
	 */
	std::vector<std::vector<std::string> > get_barcode_vector() const {
		return get_immutable(barcode_vector);
	}

	/**
	 * Get the N value for the "best N scores" mode.
	 * @return Number of "best" scores to print alignments for
	 */
	CountType get_best_n() const {
		return get_immutable(best_n);
	}

	/**
	 * Get the block size for output.
	 * @param value The block size as uint64_t (in bytes)
	 */
	uint64_t get_block_size() const {
		return get_immutable(block_size);
	}

	/**
	 * Set the block size for output.
	 * @param value The block size as String as given by the user input. Format: [0-9]+[BKM]?
	 */
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

	/**
	 * Get the compression format.
	 * @return Compression format as uint8_t.
	 */
	uint8_t get_compression_format() const {
		return get_immutable(compression_format);
	}

	/**
	 * Set the compression format.
	 * @param value Compression format as an uint16_t [0-2]
	 */
	void set_compression_format(uint16_t value) {
		if ( value > 2 )
			value = 2;
		uint8_t one_byte_value = value;
		set_immutable(compression_format, one_byte_value);
	}

	/**
	 * Get the total number of cycles.
	 * @return The total number of cycles
	 */
	CountType get_cycles() const {
		return get_immutable(cycles);
	}

	/**
	 * Set the total number of cycles.
	 * @param value The total number of cycles (must match the sum of cycles given by -r parameter)
	 */
	void set_cycles(CountType value) {
		set_immutable(cycles, value);
	}

	/**
	 * Get the frequency to increase the tolerated number of errors during the alingment algorithm.
	 * @return Frequency to tolerate more errors during the alignment algorithm.
	 * @author Tobias Loka
	 */
	CountType get_error_rate() const {
		return get_immutable(error_rate);
	}

	/**
	 * Set the frequency to increase the tolerated number of errors during the alingment algorithm.
	 * @param value Frequency to tolerate more errors during the alignment algorithm.
	 * @author Tobias Loka
	 */
	void set_error_rate(CountType value) {
		set_immutable(error_rate, value);
	}

	/**
	 * Check to use extended CIGAR format for alignment output.
	 * @return true if extended CIGAR format is activated, false is standard CIGAR is used.
	 */
	bool get_extended_cigar() const {
		return get_immutable(extended_cigar);
	}

	/**
	 * Set to use extended CIGAR format for alignment output.
	 * @param value true to activate extended CIGAR format, false to use standard CIGAR.
	 */
	void set_extended_cigar(bool value) {
		set_immutable(extended_cigar, value);
	}

	/**
	 * Get the flowcell ID.
	 * @return The flowcell ID.
	 */
	std::string get_flowcell_id() {
		return get_immutable(flowcell_id);
	}

	/**
	 * Set the flowcell ID.
	 * @param value The flowcell ID.
	 */
	void set_flowcell_id(std::string value) {
		set_immutable(flowcell_id, value);
	}

	/**
	 * Check if resorting of alignment files before writing output is activated.
	 * @return true if resorting is activated, false otherwise
	 */
	bool get_force_resort() const {
		return get_immutable(force_resort);
	}

	/**
	 * Set resorting of alignment files before writing output is activated or not.
	 * @return true to activate resorting, false otherwise
	 */
	void set_force_resort(bool value) {
		set_immutable(force_resort, value);
	}

	/**
	 * Get the path to the index (including the file prefix)
	 * @return Path to the index as string.
	 */
	std::string get_index_fname() const {
		return get_immutable(index_fname);
	}

	/**
	 * Set the path to the index (including the file prefix)
	 * @param value Path to the index as string.
	 */
	void set_index_fname(std::string value) {
		set_immutable(index_fname, value);
	}

	/**
	 * Get the instument ID.
	 * @return The instrument ID.
	 */
	std::string get_instrument_id() {
		return get_immutable(instrument_id);
	}

	/**
	 * Set the instrument ID.
	 * @param value The instrument ID.
	 */
	void set_instrument_id(std::string value) {
		set_immutable(instrument_id, value);
	}

	/**
	 * Activate keeping the temporary alignment files for all cycles.
	 */
	void set_keep_all_aln_files() {
		std::vector<CountType>keep_all_files (get_cycles());
		std::iota(keep_all_files.begin(), keep_all_files.end(), 1);
		set_keep_aln_files(keep_all_files);
	}

	/**
	 * Get the cycles to keep the alignment files for.
	 * @return Kept cycles as a vector of Integers as given as user input.
	 */
	std::vector<uint16_t> get_keep_aln_files() const {
		return get_immutable(keep_aln_files);
	}

	/**
	 * Get whether reads without a matching barcode are written to an "undetermined" file.
	 * @return true if reads without a matching barcode are written, false otherwise
	 */
	bool get_keep_all_barcodes() const {
		if ( get_barcode_vector().size() == 0 )
			return true;
		return get_immutable(keep_all_barcodes);
	}

	/**
	 * Set to write all reads that don't match a barcode to an "undetermined" file.
	 * @param value true to write reads without a matching barcode, false to not write them
	 */
	void set_keep_all_barcodes(bool value) {
		set_immutable(keep_all_barcodes, value);
	}

	/**
	 * Check to keep all sequences (also those of unmapped reads)
	 * @return true if all sequences are kept
	 */
	bool get_keep_all_sequences() const {
		return get_immutable(keep_all_sequences);
	}

	/**
	 * Set to keep all sequences (also those of unmapped reads)
	 * @param value true to keep all sequences, false otherwise
	 */
	void set_keep_all_sequences(bool value) {
		set_immutable(keep_all_sequences, value);
	}

	/**
	 * Set the cycles to keep the alignment files for.
	 * @rparam value Kept cycles as a vector of Integers as given as user input.
	 */
	void set_keep_aln_files(std::vector<uint16_t> value) {
		set_immutable(keep_aln_files, value);
	}

	/**
	 * Set the cycles to keep the alignment files for.
	 * @rparam value Kept cycles as a comma-separated std::string.
	 */
	void set_keep_aln_files(std::string value) {
		std::vector<uint16_t> keep_aln_files;
		split( keep_aln_files, value );
		set_keep_aln_files(keep_aln_files);
	}

	/**
	 * Check a cycle to be a kept cycle.
	 * @param cycle Requested cycle to check for keeping temporary files.
	 * @return true if requested cycle is a kept cycle. false otherwise
	 */
	bool is_keep_aln_files_cycle(CountType cycle) const {
		auto aln_files_cycles = get_keep_aln_files();
		if ( std::find(aln_files_cycles.begin(), aln_files_cycles.end(), cycle) == aln_files_cycles.end() )
			return false;
		return true;
	}

	/**
	 * Get the lanes to work on.
	 * @return The lanes considered for alignment as vector of Integers.
	 */
	std::vector<uint16_t> get_lanes() const {
		return get_immutable(lanes);
	}

	/**
	 * Set the lanes to work on.
	 * @param value The lanes to consider for alignment as vector of Integers.
	 */
	void set_lanes(std::vector<uint16_t> value) {
		std::sort( value.begin(), value.end() );
		value.erase( std::unique( value.begin(), value.end() ), value.end() );
		set_immutable(lanes, value);
	}

	/**
	 * Set the lanes to work on.
	 * @param value The lanes to consider for alignment as a comma-separated std::string.
	 */
	void set_lanes(std::string value) {
		std::vector<uint16_t> lanes;
		split( lanes, value );
		set_lanes(lanes);
	}

	/**
	 * Get the number of mates.
	 * @return The number of mates.
	 */
	uint16_t get_mates() const {
		return get_immutable(mates);
	}

	/**
	 * Set the tiles to work on by giving only the highest tile number.
	 * All surface, swath and tile counts lower/equal the respective part of the input value will be considered.
	 * @param value The maximum tile number to consider.
	 */
	void set_max_tile( CountType value ) {
		set_tiles ( maxTile_to_tiles( value ) );
	}

	/**
	 * Get the minimum quality for a valid base call.
	 * @return The minimum quality as a number (0-41 for Illumina PHRED scores)
	 */
	CountType get_min_qual() const {
		return get_immutable(min_qual);
	}

	/**
	 * Set the minimum quality for a valid base call.
	 * @param value The minimum quality as a number (0-41 for Illumina PHRED scores)
	 */
	void set_min_qual(CountType value) {
		set_immutable(min_qual, value);
	}

	/**
	 * Get the output mode.
	 * @return The output mode as AlignmentMode object.
	 * @author Tobias Loka
	 */
	AlignmentMode get_mode() const{
		return get_immutable(mode);
	}

	/**
	 * Check if a certain output mode is set.
	 */
	bool is_mode ( AlignmentMode alignment_mode ) const {
		return ( get_immutable(mode) == alignment_mode ) ;
	}

	/**
	 * Set the alignment mode for sequencing [(ALL)(UNIQUE)(BESTN[0-9]+)(ALLBEST)(ANYBEST)AU(N[0-9]+)HB].
	 * @param value Alignment mode as string.
	 */
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

	/**
	 * Get the number of threads used for output.
	 * @return The number of output threads. This number of threads is also used for output if new alignment jobs would be available.
	 */
	CountType get_num_out_threads() const {

		return get_immutable(num_out_threads);
	}

	/**
	 * Set the number of threads used for output.
	 * @param value The number of output threads. This number of threads is also used for output if new alignment jobs would be available.
	 */
	void set_num_out_threads(CountType value) {
		set_immutable(num_out_threads, value);
	}

	/**
	 * Get the total number of threads.
	 * @return Number of threads.
	 */
	CountType get_num_threads() const {
		return get_immutable(num_threads);
	}

	/**
	 * Set the total number of threads.
	 * @param value Number of threads.
	 */
	void set_num_threads(CountType value) {
		set_immutable(num_threads, value);
	}

	/**
	 * Check a cycle to be an output cycle.
	 * @param cycle Requested cycle to check for output.
	 * @return true if requested cycle is an output cycle. false otherwise
	 */
	bool is_output_cycle(CountType cycle) const {
		auto out_cycles = get_output_cycles();
		if ( std::find(out_cycles.begin(), out_cycles.end(), cycle) == out_cycles.end() )
			return false;
		return true;
	}

	/**
	 * Get the cycles to write output.
	 * @return Output cycles as a vector of Integers as given as user input.
	 */
	std::vector<uint16_t> get_output_cycles() const {
		return get_immutable(output_cycles);
	}

	/**
	 * Set the cycles to write output.
	 * @param value Output cycles as a vector of Integers as given as user input.
	 */
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
		set_immutable(output_cycles, the_cycles);
	}

	/**
	 * Set the cycles to write output.
	 * @param value Output cycles as a comma-separated std::string.
	 */
	void set_output_cycles(std::string value) {
		std::vector<uint16_t> output_cycles;
		split( output_cycles, value );
		set_output_cycles(output_cycles);
	}

	/**
	 * Get the output format.
	 * @return The output format set by the user.
	 */
	OutputFormat get_output_format() const {
		return get_immutable(output_format);
	}

	/**
	 * Switch output format by giving a string.
	 * @param value Output format as std::string, e.g. "BAM"
	 */
	void set_output_format(std::string value) {
		if ( value == "SAM" || value == "S" )
			set_immutable(output_format, OutputFormat::SAM);
		else if ( value == "BAM" || value == "B")
			set_immutable(output_format, OutputFormat::BAM);
		else
			throw std::runtime_error("Invalid output format: " + value + ".");
	}

	/**
	 * Get the directory to write output files.
	 * @return Path to the output directory
	 */
	std::string get_out_dir() const {
		return get_immutable(out_dir);
	}

	/**
	 * Set the directory to write output files.
	 * @param value Path to the output directory (no need to exist).
	 */
	void set_out_dir(std::string value) {
		set_immutable(out_dir, value);
	}

	/**
	 * Set the read structure of the sequencing run.
	 * @param read_argument A vector of strings in user input-like format (e.g. 101R 8B 8B 101R)
	 */
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
		set_cycles(lenSum);

	}

	void set_read_structure ( std::string read_argument ) {
		std::vector<std::string> read_structure;
		split( read_structure, read_argument );
		set_read_structure ( read_structure );
	}

	/**
	 * Check to report unmapped reads
	 * @return true if all sequences are kept
	 */
	bool get_report_unmapped() const {
		return get_immutable(report_unmapped);
	}

	/**
	 * Set to keep all sequences (also those of unmapped reads)
	 * @param value true to keep all sequences, false otherwise
	 */
	void set_report_unmapped(bool value) {
		set_immutable(report_unmapped, value);
	}

	/**
	 * Get the root directory (i.e. the path of the BaseCalls directory)
	 * @param value the path to the root directory as string
	 */
	std::string get_root() const {
		return get_immutable(root);
	}

	/**
	 * Set the root directory (i.e. the path of the BaseCalls directory)
	 * @param value the path to the root directory as string
	 */
	void set_root(std::string value) {
		set_immutable(root, value);
	}

	/**
	 * Get the run ID.
	 * @return The run ID.
	 */
	std::string get_run_id() {
		return get_immutable(run_id);
	}

	/**
	 * Set the run ID.
	 * @param value The run ID.
	 */
	void set_run_id(std::string value) {
		set_immutable(run_id, value);
	}

	/**
	 * Get the interval to create new seeds.
	 * @return The interval to create new seeds.
	 * @author Tobias Loka
	 */
	CountType get_seeding_interval() const {
		return get_immutable(seeding_interval);
	}

	/**
	 * Get the interval to create new seeds.
	 * @param value The interval to create new seeds.
	 * @author Tobias Loka
	 */
	void set_seeding_interval(CountType value) {

		// interval must be at least 1
		value = value == 0 ? 1 : value;

		set_immutable(seeding_interval, value);
	}

	/**
	 * Get the sequence elements.
	 * @return Vector of sequence elements. This includes read and barcode elements.
	 */
	std::vector<SequenceElement> get_seqs() const {
		return get_immutable(seqs);
	}

	/**
	 * Get a SequenceElement object from the seqs vector by using the mate number
	 * @param id The mate number of the SequenceElement.
	 * @return The respective SequenceElement object for the given mate number. NULLSEQ if mate==0 (barcodes).
	 * @author Tobias Loka
	 */
	SequenceElement get_seq_by_mate(CountType mate) const {
		if ( mate == 0 ) return NULLSEQ;
		auto the_seq = seqs.get();
		for (uint16_t i = 0; i != the_seq.size(); i++) {
			if(the_seq[i].mate == mate) return the_seq[i];
		}
		return NULLSEQ;
	}

	/**
	 * Get a SequenceElement object from the seqs vector by using the id
	 * @param id The id of the SequenceElement.
	 * @return The respective SequenceElement object for the given id.
	 * @author Tobias Loka
	 */
	SequenceElement get_seq_by_id(CountType id) const {
		return seqs.get()[id];
	}

	/**
	 * Get the cycle to start computations (used for --continue)
	 * @return The cycle to start computations
	 */
	CountType get_start_cycle() const {
		CountType ret = get_immutable(start_cycle);

		// Value if not set.
		if ( ret == 0 )
			return 1;

		return ret;
	}

	/**
	 * Set the cycle to start computations (used for --continue)
	 * @param value The cycle to start computations
	 */
	void set_start_cycle(CountType value) {
		set_immutable(start_cycle, value);
	}

	/**
	 * Get the directory to store temporary data.
	 * @return Path to the directory as string.
	 */
	std::string get_temp_dir() const {
		std::string dir = get_immutable(temp_dir);
		dir = dir!="" ? dir : get_immutable(root);
		return dir;
	}

	/**
	 * Set the directory to store temporary data.
	 * @param value Path to the directory as string.
	 */
	void set_temp_dir(std::string value) {
		set_immutable(temp_dir, value);
	}

	/**
	 * Get the tiles to work on.
	 * @return The tiles considered for alignment as vector of Integers.
	 */
	std::vector<uint16_t> get_tiles() const {
		return get_immutable(tiles);
	}

	/**
	 * Set the tiles to work on.
	 * @param value The tiles to consider for alignment as vector of Integers.
	 */
	void set_tiles(std::vector<uint16_t> value) {
		std::sort( value.begin(), value.end() );
		value.erase( std::unique( value.begin(), value.end() ), value.end() );
		set_immutable(tiles, value);
	}

	/**
	 * Set the tiles to work on.
	 * @param value The tiles to consider for alignment as a comma-separated std::string.
	 */
	void set_tiles(std::string value) {
		std::vector<uint16_t> tiles;
		split ( tiles, value );
		set_tiles( tiles );
	}


	/////////////////////////////////////
	////////// Scoring scheme  //////////
	/////////////////////////////////////

	/**
	 * Get the minimum alignment score for an alignment to be considered.
	 * @return Minimum score for an alignment.
	 * @author Tobias Loka
	 */
	ScoreType get_min_as() const {
		return get_immutable(min_as);
	}

	/**
	 * Set the minimum alignment score for an alignment to be considered.
	 * @param value Minimum score for an alignment.
	 * @author Tobias Loka
	 */
	void set_min_as(ScoreType value) {
		set_immutable(min_as, value);
	}

	/**
	 * Get the maximumlength of consecutive gaps (insertions or deletions).
	 * @return Maximum gap length.
	 * @author Tobias Loka
	 */
	CountType get_max_gap_length() const {
		return get_immutable(max_gap_length);
	}

	/**
	 * Set the maximumlength of consecutive gaps (insertions or deletions).
	 * @param value Maximum gap length.
	 * @author Tobias Loka
	 */
	void set_max_gap_length(CountType value) {
		set_immutable(max_gap_length, value);
	}

	/**
	 * Get the score for a match.
	 * @return Score for a match.
	 * @author Tobias Loka
	 */
	CountType get_match_score() const {
		return get_immutable(match_score);
	}

	/**
	 * Set the score for a match.
	 * @param value Score for a match.
	 * @author Tobias Loka
	 */
	void set_match_score(CountType value) {
		set_immutable(match_score, value);
	}

	/**
	 * Get the penalty for a mismatch.
	 * @return Penalty for a mismatch.
	 * @author Tobias Loka
	 */
	CountType get_mismatch_penalty() const {
		return get_immutable(mismatch_penalty);
	}

	/**
	 * Set the penalty for a mismatch.
	 * @param value Penalty for a mismatch.
	 * @author Tobias Loka
	 */
	void set_mismatch_penalty(CountType value) {
		set_immutable(mismatch_penalty, value);
	}

	/**
	 * Get the penalty for an insertion opening.
	 * @return Penalty for an insertion opening.
	 * @author Tobias Loka
	 */
	CountType get_insertion_opening_penalty() const {
		return get_immutable(insertion_opening_penalty);
	}

	/**
	 * Set the penalty for an insertion opening.
	 * @param value Penalty for an insertion opening.
	 * @author Tobias Loka
	 */
	void set_insertion_opening_penalty(CountType value) {
		set_immutable(insertion_opening_penalty, value);
	}

	/**
	 * Get the penalty for a deletion opening.
	 * @return Penalty for a deletion opening.
	 * @author Tobias Loka
	 */
	CountType get_deletion_opening_penalty() const {
		return get_immutable(deletion_opening_penalty);
	}

	/**
	 * Set the penalty for a deletion opening.
	 * @param value Penalty for a deletion opening.
	 * @author Tobias Loka
	 */
	void set_deletion_opening_penalty(CountType value) {
		set_immutable(deletion_opening_penalty, value);
	}

	/**
	 * Get the penalty for an insertion extension.
	 * @return Penalty for an insertion extension.
	 * @author Tobias Loka
	 */
	CountType get_insertion_extension_penalty() const {
		return get_immutable(insertion_extension_penalty);
	}

	/**
	 * Set the penalty for an insertion extension.
	 * @param value Penalty for an insertion extension.
	 * @author Tobias Loka
	 */
	void set_insertion_extension_penalty(CountType value) {
		set_immutable(insertion_extension_penalty, value);
	}

	/**
	 * Get the penalty for a deletion extension.
	 * @return Penalty for a deletion extension.
	 * @author Tobias Loka
	 */
	CountType get_deletion_extension_penalty() const {
		return get_immutable(deletion_extension_penalty);
	}

	/**
	 * Set the penalty for a deletion extension.
	 * @param value Penalty for a deletion extension.
	 * @author Tobias Loka
	 */
	void set_deletion_extension_penalty(CountType value) {
		set_immutable(deletion_extension_penalty, value);
	}

	/**
	 * Get the penalty for a softclip extension. This is only used for the final output, not for filtering during the alignment algorithm.
	 * @return Penalty for a softclip extension.
	 * @author Tobias Loka
	 */
	float get_softclip_extension_penalty() const {
		return get_immutable(softclip_extension_penalty);
	}

	/**
	 * Set the penalty for a softclip extension.
	 * @param value Penalty for a softclip extension.
	 * @author Tobias Loka
	 */
	void set_softclip_extension_penalty(float value) {
		set_immutable(softclip_extension_penalty, value);
	}

	/**
	 * Get the penalty for a softclip opening. This is only used for the final output, not for filtering during the alignment algorithm.
	 * @return Penalty for a softclip opening.
	 * @author Tobias Loka
	 */
	float get_softclip_opening_penalty() const {
		return get_immutable(softclip_opening_penalty);
	}

	/**
	 * Set the penalty for a softclip opening.
	 * @param value Penalty for a softclip opening.
	 * @author Tobias Loka
	 */
	void set_softclip_opening_penalty(float value) {
		set_immutable(softclip_opening_penalty, value);
	}

	/**
	 * Get the maximal relative length of the front softclip.
	 * @return Maximal softclip ratio.
	 * @author Tobias Loka
	 */
	float get_max_softclip_ratio() const {
		return get_immutable(max_softclip_ratio);
	}

	/**
	 * Set the maximal relative length of the front softclip.
	 * @param value The maximal length of the front softclip.
	 * @author Tobias Loka
	 */
	void set_max_softclip_ratio(float value) {
		set_immutable(max_softclip_ratio, value);
	}

};

#endif
