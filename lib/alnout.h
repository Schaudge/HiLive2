#ifndef ALNOUT_H
#define ALNOUT_H

#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "kindex.h"
#include "tools.h"
#include "alnstream.h"
#include "alnread.h"
#include "illumina_parsers.h"
#include "parallel.h"

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

	/**
	 * Write records to the output file in a "thread-safe" manner.
	 * @param records Reference to a vector containing an other vector of records.
	 */
	void writeRecords ( std::vector<std::vector<seqan::BamAlignmentRecord>> & records ) {

		if ( records.size() == 0 )
			return;

		lock();

		for ( auto & subvector : records)
		seqan::writeRecords(bfo, subvector);

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
	void set_context(StringListType & seq_names, std::vector<uint32_t> & seq_lengths) {
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
	 * Set the SAM fields according to the information about the other mates
	 * @param reference to a vecor of vectors for SAM records of each mate.
	 */
	void setMateSAMFlags( std::vector<std::vector<seqan::BamAlignmentRecord>> & records );

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
	 * Open an input alignment stream of the sorted align file of a mate.
	 * If a sorted align file does not exist yet or resort if forced, the sort_tile function is called to create this file.
	 * @param lane The lane.
	 * @param tile The tile.
	 * @param mateCycle Cycle of the desired mate. In general, this is NOT the sequencing cycle.
	 * @param mate The current mate.
	 * @param istream The input stream to open the file.
	 * @return 0 on success. 1 on invalid mate. 2 on matecycle equals 0. 3 on sorted align file not available.
	 */
	CountType openiAlnStream( CountType lane, CountType tile, CountType mateCycle, CountType mate, iAlnStream* istream);

	/**
	 * Open the input alignment streams of the sorted align files of all mates.
	 * @param lane The lane.
	 * @param tile The tile.
	 * @param filter_exist true, if a filter file exist.
	 * @param filter_size Number of alignments in the filter file. This value is ignored if filter_exist==false.
	 * @return Vector containing pointer to all opened input streams.
	 */
	std::vector<iAlnStream*> openiAlnStreams( CountType lane, CountType tile, bool filter_exist, unsigned filter_size);


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
