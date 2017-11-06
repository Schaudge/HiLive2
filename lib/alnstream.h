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

// Output alignment stream: write alignments to file one by one
class oAlnStream {
  // dataset information for the header
  uint16_t lane;
  uint16_t tile;
  uint16_t cycle;
  CountType rlen;
  uint32_t num_reads;

  // number of reads written to file
  uint32_t num_written;

  // data buffer (don't write everything at once)
  std::vector<char> buffer;

  // size of the buffer
  uint64_t buf_size;

  // current position in the buffer
  uint64_t buf_pos;

  // output file format
  // 0: no compression
  // 1: zlib compression (level 1)
  // 11: lz4 compression (level 1)
  uint8_t format;

  // file handles
  FILE* ofile;
  gzFile ozfile;

  // write function for lz4 compression
  uint64_t lz4write(const char* buf, uint64_t size);

 public:
  // constructor initializes all member variables
  oAlnStream(uint16_t ln, uint16_t tl, uint16_t cl, CountType rl, uint32_t nr, uint64_t bs, uint8_t fmt);

  // open Alignment stream file and write header
  uint64_t open(std::string fname);

  // writes a read alignment to the output Alignment file. 
  // Buffering is handled internally
  uint64_t write_alignment(ReadAlignment * al);
  
  // checks if the correct number of alignments was written and closes the Alignment file
  bool close();
};



// Input alignment stream: loads read alignments from a file one by one
class iAlnStream {
  // dataset information for the header
  uint16_t lane;
  uint16_t tile;
  uint16_t cycle;
  CountType rlen;
  uint32_t num_reads;

  // number of reads loaded from file
  uint32_t num_loaded;

  // data buffer (read blocks of data)
  std::vector<char> buffer;

  // size of the buffer
  uint64_t buf_size;

  // current position in the buffer
  uint64_t buf_pos;

  // output file format
  // 0: no compression
  // 1: zlib compression (level 1)
  // 11: lz4 compression (level 1)
  uint8_t format;

  // file pointer
  FILE* ifile;
  gzFile izfile;

  // read function for LZ4 compression. Reads one block of data to buffer
  uint64_t lz4read_block();

 public:
  // constructor initializes only block size and file format
  iAlnStream(uint64_t bs, uint8_t fmt);

  // open Alignment stream file and load header
  uint64_t open(std::string fname);

  // Try to open Alignment stream file. If not successful, create a copy of the file and work on the copy.
  uint64_t open_c(std::string fname, std::string c_ext_name);

  // loads a read alignment from the input Alignment file. 
  // Buffering is handled internally
  ReadAlignment* get_alignment();
  
  // checks if the correct number of alignments was loaded and closes the Alignment file
  bool close();

  // get dataset information
  inline uint16_t get_lane() {return lane;};
  inline uint16_t get_tile() {return tile;};
  inline uint16_t get_cycle() {return cycle;};
  inline CountType get_rlen() {return rlen;};
  inline uint32_t get_num_reads() {return num_reads;};
  inline uint32_t get_num_loaded() {return num_loaded;};
};


//-------------------------------------------------------------------//
//------  The StreamedAlignment class  ------------------------------//
//-------------------------------------------------------------------//

class StreamedAlignment {
  
  // dataset information
  uint16_t lane;
  uint16_t tile;
  CountType rlen;

  // fetch the next read from the input stream
  ReadAlignment get_next_read();
  
  // write an alignment to the output stream
  uint64_t write_alignment(ReadAlignment& ral);

  // get the path to the bcl file of a given cycle
  std::string get_bcl_file(uint16_t cycle, uint16_t read_number);

  // get the path to the alignment file. The alignment file is located in
  // <base>/L00<lane>/s_<lane>_<tile>.<cycle>.align
  // if base == "": base = globalAlignmentSettings.get_root()
  std::string get_alignment_file(uint16_t cycle, uint16_t mate, std::string base = "");

  // get the path to the filter file. The illumina filter information is located in
  // <root>/L00<lane>/s_<lane>_<tile>.filter
  std::string get_filter_file();

 public:
  StreamedAlignment(uint16_t ln, uint16_t tl, CountType rl): lane(ln), tile(tl), rlen(rl) {};  

  StreamedAlignment& operator=(const StreamedAlignment& other);
  
  // create directories required to store the alignment files (only if not stored in root)
  void create_directories();

  // initialize empty alignment. Creates files for a virtual Cycle 0
  void init_alignment(uint16_t mate);
  
  // extend an existing alignment from cycle <cycle-1> to <cycle>
  uint64_t extend_alignment(uint16_t cycle, uint16_t read_no, uint16_t mate, KixRun* index);

  /**
   * Extend the barcode for all reads with the information of the current sequencing cycle.
   * @param bc_cycle The cycle of the barcode read.
   * @param read_cycle The last handled cycle for the respective mate (should always be 0 or the full length)
   * @param read_no The number of the sequence read for which the barcode will be extended (:= index in globalAlignmentSettings.seqs).
   * @param mate The read mate to extend the barcode.
   * @return
   * @author Tobias Loka
   */
  void extend_barcode(uint16_t bc_cycle, uint16_t read_cycle, uint16_t read_no, uint16_t mate);

}; /* END class StreamedAlignment */




//-------------------------------------------------------------------//
//------  Streamed SAM generation -----------------------------------//
//-------------------------------------------------------------------//

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
	std::mutex tasks_mutex;

	/** Mutex to lock the finalization. */
	std::mutex finalizing;

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

	bool finalized = false;

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
	bool add_task( Task t, ItemStatus status ) {
		std::lock_guard<std::mutex> lock(tasks_mutex);
		if ( tasks.find(t) != tasks.end() )
			return false;
		tasks[t] = status;
		return true;
	}

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
