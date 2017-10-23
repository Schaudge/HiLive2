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
	~Atomic_bfo ( ) {  }

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

	void set_context(StringListType & seq_names, std::vector<uint32_t> seq_lengths) {
		seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > rnc(refNames);
		refNamesCache = rnc;
		seqan::BamIOContext<seqan::StringSet<seqan::CharString> > cxt(refNames, refNamesCache);
		context = cxt;
		seqan::contigNames(context) = seq_names;
		seqan::contigLengths(context) = seq_lengths;
	}

	void emplace_back ( std::string f_name ) {
		bfos.emplace_back ( f_name.c_str() );
		bfos.back().setContext(context);
	}

	Atomic_bfo & back() { return bfos.back(); }

	Atomic_bfo & operator [](int i) { return bfos[i]; }

};


class AlnOut {

// private:
public:

	std::map<Task, ItemStatus> tasks;
	std::mutex tasks_mutex;
	CountType cycle;

	std::deque<std::thread> threads;

	BamFileOutDeque bfos;

//	std::deque<seqan::BamFileOut> bfos;
//	std::deque<std::mutex> bfo_mutexes;
//	seqan::BamIOContext<seqan::StringSet<seqan::CharString> > bamIOContext;
//	seqan::StringSet<seqan::CharString> referenceNames;
//	seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > referenceNamesCache;


	std::vector<std::string> barcodes;
	std::vector<CountType> mateCycles;
	std::vector<CountType> min_as_scores;

	std::vector<std::vector<seqan::BamAlignmentRecord>> records_buffer;


	KixRun* index;

	Task get_next( ItemStatus getStatus, ItemStatus setToStatus ) {
		std::lock_guard<std::mutex> lock(tasks_mutex);
		for ( auto it = tasks.begin(); it != tasks.end(); ++it ) {
			if ( it->second == getStatus ) {
				tasks[it->first] = setToStatus;
				return it->first;
			}
		}
		return NO_TASK;
	}

	CountType get_task_status_num ( ItemStatus getStatus ) {
		CountType num = 0;
		std::lock_guard<std::mutex> lock(tasks_mutex);
		for ( auto it = tasks.begin(); it != tasks.end(); ++it ) {
			if ( it->second == getStatus ) {
				num += 1;
			}
		}
		return num;
	}

	bool set_task_status( Task t, ItemStatus status ) {
		std::lock_guard<std::mutex> lock(tasks_mutex);
		if ( tasks.find(t) == tasks.end() )
			return false;
		tasks[t] = status;
		return true;
	}

	bool set_task_status_from_to( Task t, ItemStatus oldStatus, ItemStatus newStatus ) {
		std::lock_guard<std::mutex> lock(tasks_mutex);
		if ( tasks.find(t) == tasks.end() )
			return false;
		if ( tasks[t] == oldStatus ) {
			tasks[t] = newStatus;
			return true;
		}
		return false;
	}

	bool add_task( Task t, ItemStatus status ) {
		std::lock_guard<std::mutex> lock(tasks_mutex);
		if ( tasks.find(t) != tasks.end() )
			return false;
		tasks[t] = status;
		return true;
	}

	std::string getOutFileName( CountType barcodeIndex, CountType cycle );
	std::string getTempOutFileName( CountType barcodeIndex, CountType cycle );
	std::string getTileOutFileName ( CountType ln, CountType tl, CountType mate, CountType cycle );
	std::string getTempTileOutFileName ( CountType ln, CountType tl, CountType mate, CountType cycle );

	void write_tile_to_bam ( Task t ) {
		try {
			__write_tile_to_bam__ (t);
			set_task_status( t, FINISHED );
		} catch ( const std::exception& e) {
			set_task_status( t, FAILED );
			std::cerr << "Writing of task " << t << " failed: " << e.what() << std::endl;
		}
	}

	void __write_tile_to_bam__ ( Task t );


public:

	AlnOut (std::vector<CountType> lns, std::vector<CountType> tls, CountType cycl, KixRun* idx);
	~AlnOut ();
	bool is_finished() {
		return ( get_task_status_num( FINISHED ) + get_task_status_num( FAILED ) ) == tasks.size();
	};
	bool sort_tile ( CountType ln, CountType tl, CountType mate, CountType cycle, bool overwrite = false );

	bool task_available ( Task t ) {
		return set_task_status_from_to( t, WAITING, BCL_AVAILABLE );
	}

	Task write_next ( ) {

		Task t = get_next ( BCL_AVAILABLE, RUNNING );
		if ( t != NO_TASK )
			threads.emplace_back(&AlnOut::write_tile_to_bam, this, t);

		return t;
	};

	CountType get_num_threads() { return threads.size(); }

//	void write_records(CountType min_num_records = 1) {
//		for ( CountType i = 0; i < records_buffer.size(); i++ ) {
//			std::vector<seqan::BamAlignmentRecord> buffer;
//			if ( records_buffer[i].size() >= min_num_records ) {
//				buffer.swap(records_buffer[i]);
//				bfos[i].writeRecords(buffer);
//			}
//		}
//	}

	void join() {
		for ( auto& t : threads )
			t.join();
	}

};

//uint64_t alignments_to_sam(std::vector<uint16_t> lns, std::vector<uint16_t> tls, KixRun* index, CountType cycle);

#endif /* ALNSTREAM_H */
