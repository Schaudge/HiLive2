#ifndef ALNSTREAM_H
#define ALNSTREAM_H

#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "kindex.h"
#include "tools.h"
#include "alnread.h"
#include "illumina_parsers.h"

/** Output alignment stream. Write alignments to file one by one. */
class oAlnStream {

	// --- HEADER INFORMATION ---

	/** Number of the sequencing lane.*/
	uint16_t lane;
	/** Number of the sequencing tile.*/
	uint16_t tile;
	/** Current alignment cycle.*/
	uint16_t cycle;
	/** Root directory.*/
	std::string root;
	/** Total length of the alignment.*/
	CountType rlen;
	/** Number of reads on the current lane.tile.*/
	uint32_t num_reads;


	// --- TECHNICAL SETTINGS ---

	/** Number of written reads. */
	uint32_t num_written;
	/** Data buffer.*/
	std::vector<char> buffer;
	/** Size of the data buffer.*/
	uint64_t buf_size;
	/** Current position in the buffer.*/
	uint64_t buf_pos;

	/** Compression of the output.
	 * 0: no compression
	 * 1: zlib compression (level 1)
	 * 11: lz4 compression (level 1)
	 */
	uint8_t format;

	/** Output file stream. */
	FILE* ofile;
	/** Compression file stream. */
	gzFile ozfile;

	/** Function for writing lz4 compression. */
	uint64_t lz4write(const char* buf, uint64_t size);

 public:
	/** Constructor for output alignment stream.
	 * @param ln The lane number.
	 * @param tl The tile number.
	 * @param cl The current cycle.
	 * @param rt The root directory.
	 * @param rl Total alignment length.
	 * @param nr Number of reads on this lane.tile.
	 * @param bs Buffer size.
	 * @param fmt Compression format
	 */
	oAlnStream(uint16_t ln, uint16_t tl, uint16_t cl, std::string rt, CountType rl, uint32_t nr, uint64_t bs, uint8_t fmt);

	/** Open the output file and write the header.
	 * @param fname Output file name
	 * @return Bytes written (header)
	 */
	uint64_t open(std::string fname);

	/** Write read alignment to output file. Buffer is handled internally.
	 * @param al The read alignment to write to the output file
	 * @return Bytes written
	 */
	uint64_t write_alignment(ReadAlignment * al);

  
	/**
	 * Close the output stream.
	 * @return true, if successful and correct number of bytes written. False otherwise.
	 */
	bool close();
};



/** Input alignment stream. Load alignments to file one by one. */
class iAlnStream {

	// --- HEADER INFORMATION ---

	/** Number of the sequencing lane.*/
	uint16_t lane;
	/** Number of the sequencing tile.*/
	uint16_t tile;
	/** Current alignment cycle.*/
	uint16_t cycle;
	/** Root directory.*/
	std::string root;
	/** Total length of the alignment.*/
	CountType rlen;
	/** Number of reads on the current lane.tile.*/
	uint32_t num_reads;


	// --- TECHNICAL SETTINGS ---

	/** Number of loaded reads. */
	uint32_t num_loaded;
	/** Data buffer.*/
	std::vector<char> buffer;
	/** Size of the data buffer.*/
	uint64_t buf_size;
	/** Current position in the buffer.*/
	uint64_t buf_pos;

	/** Compression of the output.
	 * 0: no compression
	 * 1: zlib compression (level 1)
	 * 11: lz4 compression (level 1)
	 */
	uint8_t format;

	/** Input file stream. */
	FILE* ifile;
	/** Compression file stream. */
	gzFile izfile;

	/**
	 * Function for loading lz4 compression. Loads one block of data to buffer.
	 * @return Size of loaded data in bytes
	 */
	uint64_t lz4read_block();


 public:

	/**
	 * Constructor. Initialize block size and file format.
	 * @param bs Block size
	 * @param fmt Compression format
	 */
	iAlnStream(uint64_t bs, uint8_t fmt);

	/**
	 * Open alignment input stream and load the header.
	 * @param Name of the input alignment file.
	 * @return Number of bytes loaded.
	 */
	uint64_t open(std::string fname);

	/** Load read alignment from input file. Buffer is handled internally.
	 * @return Pointer to a new ReadAlignment object.
	 */
	ReadAlignment* get_alignment();

	/**
	 * Close the input stream.
	 * @return true, if successful and correct number of bytes loaded. False otherwise.
	 */
	bool close();


	// --- GETTERS ---

	/** Get lane number */
	inline uint16_t get_lane() {return lane;};
	/** Get tile number */
	inline uint16_t get_tile() {return tile;};
	/** Get cycle number */
	inline uint16_t get_cycle() {return cycle;};
	/** Get root directory */
	inline std::string get_root() {return root;};
	/** Get total length of the alignment */
	inline CountType get_rlen() {return rlen;};
	/** Get number of reads */
	inline uint32_t get_num_reads() {return num_reads;};
	/** Get number of loaded reads */
	inline uint32_t get_num_loaded() {return num_loaded;};

};


//-------------------------------------------------------------------//
//------  The StreamedAlignment class  ------------------------------//
//-------------------------------------------------------------------//

/** Object organizing a single alignment cycle for one <lane>.<tile> */
class StreamedAlignment {
  
	/** Number of the sequencing lane.*/
  	uint16_t lane;
  	/** Number of the sequencing tile.*/
  	uint16_t tile;
  	/** Root directory.*/
  	std::string root;
  	/** Total length of the alignment.*/
  	CountType rlen;

  	/**
  	 * Fetch the next read from the input stream.
  	 * @return Next read alignment obtained from the input stream
  	 */
  	ReadAlignment get_next_read();

 	/**
  	 * Write an alignment to the output stream
  	 * @param ral Read alignment to write to file
  	 * @return Number of bytes written.
  	 */
  	uint64_t write_alignment(ReadAlignment& ral);

  	/**
  	 * Get the path to the bcl file of a given cycle.
  	 * @param cycle The current cycle
  	 * @param settings The alignment settings
  	 * @param read_number The fragment number of the current sequence (mate, barcode, ..)
  	 * @return Path to the BCL file for the given cycle
  	 */
  	std::string get_bcl_file(uint16_t cycle, uint16_t read_number);

  	/**
  	 * Get the path to the alignment file
  	 * @param cycle The current cycle
  	 * @param mate The current mate
  	 * @param base Path to the temp directory. If not set (no temp directory), it is determined automatically.
  	 * @return Path to the alignment file for the given cycle
  	 */
  	std::string get_alignment_file(uint16_t cycle, uint16_t mate, std::string base = "");

  	/**
  	 * Get the path to the illumina filter file
  	 * @return Path to the illumina filter file
  	 */
  	std::string get_filter_file();

 public:

 	/**
  	 * Constructor for the alignment stream.
	 * @param ln The lane number.
	 * @param tl The tile number.
	 * @param rt The root directory.
	 * @param rl Total alignment length.
	 */
  	StreamedAlignment(uint16_t ln, uint16_t tl, std::string rt, CountType rl): lane(ln), tile(tl), root(rt), rlen(rl) {};

  	/**
  	 * Define when two alignment streams are equal.
  	 */
  	StreamedAlignment& operator=(const StreamedAlignment& other);

  	/**
  	 * Create directories required to store the alignment files (only if not stored in root)
  	 * @param settings The alignment settings
  	 */
  	void create_directories();

  	/**
  	 * Initialize empty alignment. Creates files for a virtual Cycle 0.
  	 * @param mate The current mate number
  	 * @param settings The alignment settings
  	 */
  	void init_alignment(uint16_t mate);

  	/**
  	 * Extend an existing alignment from cycle <cycle-1> to <cycle>
  	 * @param cycle The last cycle
  	 * @param read_no The number of reads
  	 * @param mate The current mate numer
  	 * @param index The alignment index
  	 * @param settings The alignment settings
  	 * @return Number of seeds after the current cycle
  	 */
  	uint64_t extend_alignment(uint16_t cycle, uint16_t read_no, uint16_t mate, KixRun* index);

  	/**
  	 * Extend the barcode for all reads with the information of the current sequencing cycle.
  	 * @param bc_cycle The cycle of the barcode read.
  	 * @param read_cycle The last handled cycle for the respective mate (should always be 0 or the full length)
  	 * @param read_no The number of the sequence read for which the barcode will be extended (:= index in AlignmentSettings::seqs).
  	 * @param mate The read mate to extend the barcode.
  	 * @param settings Object containing the program settings.
  	 * @author Tobias Loka
  	 */
  	void extend_barcode(uint16_t bc_cycle, uint16_t read_cycle, uint16_t read_no, uint16_t mate);


}; /* END class StreamedAlignment */




//-------------------------------------------------------------------//
//------  Streamed SAM generation -----------------------------------//
//-------------------------------------------------------------------//

/**
 * Write alignments to a SAM output file
 * @param ln The lane number
 * @param tl The tile number
 * @param rt Path tothe root directory
 * @param rl Total length of the alignment
 * @param mate The mate number
 * @param index The alignment index
 * @param settings The alignment settings
 * @return 0 if successful, different number otherwise
 */
uint64_t alignments_to_sam(std::vector<uint16_t> lns, std::vector<uint16_t> tls, KixRun* index, CountType cycle);

#endif /* ALNSTREAM_H */
