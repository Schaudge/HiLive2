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




/**
 * Extension of SeqAn's BamFileOut data type supporting multithreading via an atomic flag.
 * @author Tobias Loka
 */
struct Atomic_bfo {

private:

	/** The BamFileOut stream.*/
	seqan::BamFileOut bfo;

	/** The atomic flag to perform a spinlock while writing.*/
	std::atomic_flag flag = ATOMIC_FLAG_INIT;

	/**
	 * Lock the atomic flag.
	 */
	void lock() {
		while ( flag.test_and_set(std::memory_order_acquire))
			; // spin
	}

	/**
	 * Unlock the atomic flag.
	 */
	void unlock() {
		flag.clear(std::memory_order_release);
	}

public:

	/**
	 * Default constructor that also initializes the BamFileOut stream.
	 * @param f_name Output file name.
	 * */
	Atomic_bfo( std::string f_name) : bfo(f_name.c_str()) { }


	/** Destructor. */
	~Atomic_bfo ( ) {	}

	/**
	 * Set the context of the BamFileOut stream.
	 * @param context The new context.
	 */
	void setContext ( seqan::BamIOContext<seqan::StringSet<seqan::CharString> > & context ) {
		bfo.context = context;
	}

	/**
	 * Write the header to the output file.
	 * @param header The header for the output file.
	 */
	void writeHeader ( seqan::BamHeader & header ) {
		seqan::writeHeader( bfo, header );
	}

	/**
	 * Write records to the output file in a "thread-safe" manner.
	 * @param records Reference to a vector containing a set of records.
	 */
	void writeRecords ( std::vector<seqan::BamAlignmentRecord> & records ) {

		if ( records.size() == 0 )
			return;

		lock();

		seqan::writeRecords(bfo, records);

		unlock();

	}

};

/**
 * Extends a deque of atomic BamFileOut streams.
 * Store the context, refNames and refNamesCache such that it exist only once and will not be destructed as long as it is needed.
 * @author Tobias Loka
 */
class BamFileOutDeque {

	/** Deque of bfos. */
	std::deque<Atomic_bfo> bfos;

	/** All fields needed for the context of the BamFileOut streams. */
	seqan::BamIOContext<seqan::StringSet<seqan::CharString> > context;
	seqan::StringSet<seqan::CharString> refNames;
	seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > refNamesCache;

public:

	/**
	 * Set the context of the list.
	 * @param seq_names Names of all sequences in the database
	 * @param seq_length Lengths of all sequences in the database
	 */
	void set_context(StringListType & seq_names, std::vector<uint32_t> seq_lengths) {
		seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > rnc(refNames);
		refNamesCache = rnc;
		seqan::BamIOContext<seqan::StringSet<seqan::CharString> > cxt(refNames, refNamesCache);
		context = cxt;
		seqan::contigNames(context) = seq_names;
		seqan::contigLengths(context) = seq_lengths;
	}

	/**
	 * Add a Bam Output stream for the given file name.
	 * The stream will be created in-place.
	 * @param f_name File name for the Bam output stream.
	 */
	void emplace_back ( std::string f_name ) {
		bfos.emplace_back ( f_name.c_str() );
		bfos.back().setContext(context);
	}

	/**
	 * Get a reference to the last atomic bfo of the deque.
	 * @return Reference to the last atomic bfo of the deque.
	 */
	Atomic_bfo & back() { return bfos.back(); }

	/**
	 * Get the Atomic_bfo at a certain position of the deque.
	 * @param i Index of the Atomic_bfo in the deque.
	 * @return reference to the Atomic_bfo at a certain position of the deque.
	 */
	Atomic_bfo & operator [](int i) { return bfos[i]; }

	/**
	 * Clear the deque of Atomic bfos.
	 */
	void clear () {
		bfos.clear();
	}
};




/**
 * Class to organize the output for a specific cycle.
 * @author Tobias Loka
 */
class AlnOut {

private:

	/** Map of tasks and the status. */
	std::map<Task, ItemStatus> tasks;

	/** Mutex to lock tasks when their status is getting modified. */
	std::mutex tasks_lock;

	/** True if the output was finalized. */
	bool finalized = false;

	/** True if the output was initialized. */
	bool initialized = false;

	/** Mutex to lock the initialization and finalization. */
	std::mutex if_lock;

	/** Cycle for the output. */
	CountType cycle;

	/** Deque of output streams. */
	BamFileOutDeque bfos;

	/** Vector of barcodes. */
	std::vector<std::string> barcodes;

	/** Vector containing the current cycle of all mates. */
	std::vector<CountType> mateCycles;

	/** Minimal alignment score for a certain cycle to print an alignment. */
	std::vector<CountType> min_as_scores;

	/** The underlying index for the output. */
	KixRun* index;

	/**
	 * Set the status of a task (only if the task exists).
	 * @param t The task.
	 * @param status The new status.
	 * @return true, if the task exists and the status was successfully changed. false otherwise.
	 */
	bool set_task_status( Task t, ItemStatus status );

	/**
	 * Change the status of a task if it had a specified status before.
	 * @param t The task.
	 * @param oldStatus The previous status of the task.
	 * @param newStatus The new status of the task.
	 * @return true, if the task exists, has the correct previous status and the status was successfully changed. false otherwise.
	 */
	bool set_task_status_from_to( Task t, ItemStatus oldStatus, ItemStatus newStatus );

	/**
	 * Change the status of the next task with a given status.
	 * @param getStatus Status to be changed.
	 * @param setToStatus Status the task is set to.
	 * @return The task for that the status was changed.
	 */
	Task get_next( ItemStatus getStatus, ItemStatus setToStatus );

	/**
	 * Add a new task with a certain status. If the task already exists, the status will not be changed.
	 * @param t The new task.
	 * @param status The status.
	 * @return true, if the task didn't exist before and was successfully created.  false otherwise.
	 */
	bool add_task( Task t, ItemStatus status );

	/**
	 * Create a temporary align file that is sorted by score.
	 * @param ln The lane.
	 * @param tl The tile.
	 * @param mate The mate.
	 * @param cycle The mate cycle (not the sequencing cycle).
	 * @param overwrite If true and a sorted file already exists, it will be sorted again and the old file will be overridden (default: false)
	 * @return true, if sorting was successful.
	 */
	bool sort_tile ( CountType ln, CountType tl, CountType mate, CountType cycle, bool overwrite = false );

	/**
	 * Calls __write_tile_to_bam__(Task t) to start the output of a task with handled exceptions.
	 * @param t Task that contains the information about lane and tile.
	 */
	void write_tile_to_bam ( Task t );

	/**
	 * Start the output of a task. Should be called by write_tile_to_bam(Task t) to handle exceptions in an appropriate manner.
	 * @param t Task that contains the information about lane and tile.
	 */
	void __write_tile_to_bam__ ( Task t );

	/**
	 * Initalize the output streams.
	 */
	void init();

public:

	/**
	 * Constructor.
	 * Includes opening all output streams and writing the header.
	 * @param lns The lanes to consider.
	 * @param tls The tiles to consider.
	 * @param cycl The sequencing cycle.
	 * @param idx The underlying index.
	 */
	AlnOut (std::vector<CountType> lns, std::vector<CountType> tls, CountType cycl, KixRun* idx);

	/**
	 * Destructor.
	 * The destructor will wait for all running threads to be finished.
	 * The temporary output files are moved to their final locations.
	 */
	~AlnOut ();

	/**
	 * Check if the output writing is finished.
	 * @return true, if all tasks are finished or failed. false otherwise.
	 * TODO: maybe change the return value to the number of failed tasks!?
	 */
	bool is_finished() {
		if ( is_finalized() )
			return true;
		return ( get_task_status_num( FINISHED ) + get_task_status_num( FAILED ) ) == tasks.size();
	};

	/**
	 * Set that the alignment file of a task is available.
	 * @param t The task containing information about lane and tile.
	 * @return true if the task was set to available status.
	 */
	bool set_task_available ( Task t ) {
		return set_task_status_from_to( t, WAITING, BCL_AVAILABLE );
	}

	/**
	 * Check if a task is contained in the list of tasks.
	 * @param t The task of interest.
	 * @return true, if the task is in the list of tasks (the status is not considered). False if not.
	 */
	bool has_task ( Task t ) {
		if ( tasks.find(t) != tasks.end() )
			return true;
		return false;
	}

	/**
	 * Check if the status of a task is FINISHED.
	 * @param t The task of interest.
	 * @return true, if the task has status FINISHED. false, if the task has a status other than FINISHED or does not exist.
	 */
	bool is_task_finished ( Task t ) {
		if ( !has_task(t) )
			return false;
		return tasks[t] == FINISHED;
	}

	/**
	 * Check if the status of a task is FAILED.
	 * @param t The task of interest.
	 * @return true, if the task has status FAILED. false, if the task has a status other than FAILED or does not exist.
	 */
	bool is_task_failed ( Task t ) {
		if ( !has_task(t) )
			return false;
		return tasks[t] == FAILED;
	}

	/**
	 * Write the next available task. This function is implemented in a thread-safe manner.
	 * @return The task that was written. NO_TASK if no task was written.
	 */
	Task write_next ( );

	/**
	 * Get the number of tasks with a given status.
	 * @param getStatus Requested status.
	 * @return Number of tasks with the requested status.
	 */
	CountType get_task_status_num ( ItemStatus getStatus );

	/**
	 * Renames the temporary sam files if all tasks are finished.
	 * @return true, if finished with success.
	 */
	bool finalize();

	/**
	 * Check if this output controller was already finalized.
	 * @return true, if output controller was already finalized.
	 */
	bool is_finalized () {
		return finalized;
	}

	/**
	 * Check if the output controller was already initialized (i.e., the output streams are created).
	 * @return true, if the output controller was already initialized.
	 */
	bool is_initialized(){
		return initialized;
	}

	/**
	 * Get the cycle number for this output controller.
	 * @return The sequencing cycle.
	 */
	CountType get_cycle() {
		return cycle;
	}

	/**
	 * Get the cycle for a certain mate.
	 * @param Mate number.
	 * @return Mate cycle.
	 */
	CountType get_mate_cycle( CountType mate ) {
		return mateCycles[mate-1];
	}

};

#endif /* ALNSTREAM_H */
