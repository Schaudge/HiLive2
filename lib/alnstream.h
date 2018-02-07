#ifndef ALNSTREAM_H
#define ALNSTREAM_H

#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "kindex.h"
#include "tools.h"
#include "alnread.h"
#include "illumina_parsers.h"
#include "parallel.h"


/**
 * Output stream to write temporary .align files.
 */
class oAlnStream {

	/** Lane for the output. */
	uint16_t lane;

	/** Tile for the output. */
	uint16_t tile;

	/** Output cycle. */
	uint16_t cycle;

	/** Total read length. */
	CountType rlen;

	/** Total number of reads for this lane/tile. */
	uint32_t num_reads;

	/** Number of reads written to file. */
	uint32_t num_written;

	/** Data buffer. */
	std::vector<char> buffer;

	/** Size of the data buffer. */
	uint64_t buf_size;

	/** Current position in the buffer. */
	uint64_t buf_pos;

	/** Output file compression [0: None; 1: zlib (lvl1); 2: lz4 (lvl1)] */
	uint8_t format;

	/** Standard file handler. */
	FILE* fstream;

	/** zlib file handler. */
	gzFile zfstream;

	/** Name of the file that is currently streamed. */
	std::string fname;

	/** File lock flag. True, if the file was locked by this alignment stream. */
	bool flocked;

	/**
	 * Write with lz4 compression.
	 * @param buf Pointer to the buffer data.
	 * @param size Size of the buffer data.
	 * @return New size of the buffer data.
	 */
	uint64_t lz4write(const char* buf, uint64_t size);

	/** Lock the file (globally in the program). */
	void flock();

	/** Unlock the file (globally in the program). */
	void funlock();


public:

	/**
	 * Constructor.
	 * @param ln Lane for the output.
	 * @param tl Tile for the output.
	 * @param cl Output cycle.
	 * @param rl Total read length.
	 * @param nr Total number of reads.
	 * @param bs Buffer size.
	 * @param fmt Compression format (0: None, 1: zlib, 2: lz4)
	 */
	oAlnStream(uint16_t ln, uint16_t tl, uint16_t cl, CountType rl, uint32_t nr, uint64_t bs, uint8_t fmt);

	/** Default destructor. Unlocks the global fileLock if it was locked by this stream. */
	~oAlnStream();

	/**
	 * Open alignment stream for a file and write the header.
	 * The file will be locked globally.
	 * @param f_name Name of the file to be written (will be overridden if already exists).
	 * @return Number of written bytes.
	 */
	uint64_t open(std::string f_name);

	/**
	 * Write a read alignment to the current output stream.
	 * @param al Pointer to the read alignment object that is written to the output stream.
	 * @return Number of written bytes.
	 */
	uint64_t write_alignment(ReadAlignment * al);

	/**
	 * Close the file stream if all alignments were written.
	 * Unlocks the global file lock.
	 * @return true, if file stream was closed successfully.
	 */
	bool close();
};


/**
 * Input stream to read temporary .align files.
 */
class iAlnStream {

	/** Lane for the output. */
	uint16_t lane;

	/** Tile for the output. */
	uint16_t tile;

	/** Output cycle. */
	uint16_t cycle;

	/** Total read length. */
	CountType rlen;

	/** Total number of reads for this lane/tile. */
	uint32_t num_reads;

	/** Number of reads loaded from the file. */
	uint32_t num_loaded;

	/** Data buffer. */
	std::vector<char> buffer;

	/** Size of the data buffer. */
	uint64_t buf_size;

	/** Current position in the buffer. */
	uint64_t buf_pos;

	/** Output file compression [0: None; 1: zlib (lvl1); 2: lz4 (lvl1)] */
	uint8_t format;

	/** Standard file handler. */
	FILE* fstream;

	/** zlib file handler. */
	gzFile zfstream;

	/** Name of the file that is currently streamed. */
	std::string fname;

	/** File lock flag. True, if the file was locked by this alignment stream. */
	bool flocked;

	/**
	 * Load a lz4-compressed block to the buffer.
	 * @return The new buffer size.
	 */
	uint64_t lz4read_block();

	/** Lock the file (globally in the program). */
	void flock();

	/** Unlock the file (globally in the program). */
	void funlock();

public:

	/**
	 * Constructor.
	 * @param bs Buffer size.
	 * @param fmt Compression format (0: None, 1: zlib, 2: lz4)
	 */
	iAlnStream(uint64_t bs, uint8_t fmt);

	/** Default destructor. Unlocks the global fileLock if it was locked by this stream. */
	~iAlnStream();

	/**
	 * Open alignment stream for a file and load the header.
	 * The file will be locked globally.
	 * @param f_name Name of the file to be loaded.
	 * @return Number of loaded bytes.
	 */
	uint64_t open(std::string f_name);

	/**
	 * Loasd a read alignment from the current input stream.
	 * @return The read alignment object that was loaded from the file.
	 */
	ReadAlignment* get_alignment();

	/**
	 * Close the file stream if all alignments were loaded.
	 * Unlocks the global file lock.
	 * @return true, if file stream was closed successfully.
	 */
	bool close();

	// Getter
	uint16_t get_lane() {return lane;};
	uint16_t get_tile() {return tile;};
	uint16_t get_cycle() {return cycle;};
	CountType get_rlen() {return rlen;};
	uint32_t get_num_reads() {return num_reads;};
	uint32_t get_num_loaded() {return num_loaded;};
};


/**
 * Streamer for new base calls to the alignment algorithm.
 */
class StreamedAlignment {
  
  /** The lane to be handled. */
  uint16_t lane;

  /** The tile to be handled. */
  uint16_t tile;

  /** Total read length. */
  CountType rlen;

  /**
   * Get the path to the bcl file of a given cycle.
   * @param cycle The current read cycle.
   * @param mate Number of the current mate.
   * @return Path to the bcl file.
   */
  std::string get_bcl_file(uint16_t cycle, uint16_t mate);

  /**
   * Get the path to the align file.
   * @param cycle The current read cycle.
   * @param mate Number of the current mate.
   * @param base Base of the path to the align files.
   * @return Path to the align file.
   */
  std::string get_alignment_file(uint16_t cycle, uint16_t mate, std::string base = "");

  /**
   * Get the path to the filter file.
   * @return Path to the filter file.
   */
  std::string get_filter_file();

 public:

  /**
   * Constructor.
   * @param ln The lane to be handled.
   * @param tl The tile to be handled.
   * @param rl Total read length.
   */
  StreamedAlignment(uint16_t ln, uint16_t tl, CountType rl): lane(ln), tile(tl), rlen(rl) {};  

  /**
   * Create the underlying directories of the align files.
   */
  void create_directories();

  /**
   * Initialize empty alignments for the current mate (stored as output of a virtual cycle 0).
   * @param mate Number of the current mate.
   */
  void init_alignment(uint16_t mate);
  
  /**
   * Extend the alignments for all reads of the specified lane and tile by one cycle.
   * @param cycle Current cycle, i.e. the cycle that will be extended.
   * @param read_no Total number of reads.
   * @param mate Number of the current mate.
   * @param index Pointer to the reference index.
   * @return Total number of seeds (for all reads).
   */
  uint64_t extend_alignment(uint16_t cycle, uint16_t read_no, uint16_t mate, KixRun* index);

  /**
   * Extend the barcode for all reads with the information of the current sequencing cycle.
   * @param bc_cycle The cycle of the barcode read.
   * @param read_cycle The last handled cycle for the respective mate (should always be 0 or the full length)
   * @param read_no The number of the sequence read for which the barcode will be extended (:= index in globalAlignmentSettings.seqs).
   * @param mate The read mate to extend the barcode.
   */
  void extend_barcode(uint16_t bc_cycle, uint16_t read_cycle, uint16_t read_no, uint16_t mate);

  StreamedAlignment& operator=(const StreamedAlignment& other);

};

#endif /* ALNSTREAM_H */
