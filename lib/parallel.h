#ifndef PARALLEL_H
#define PARALLEL_H

#include "headers.h"
#include "definitions.h"
#include "tools.h"
#include "kindex.h"

//------ Threading tools --------------------------------------------//

/** Task data structure. Contains all information for a thread to process a BCL file.
 * @author Martin Lindner
 */
struct Task {
  // dataset information
  uint16_t lane;
  uint16_t tile;
  SequenceElement seqEl;
//  uint16_t read_number; // Number of the read fragment (:= index in settings.seqLength array)
  uint16_t cycle; // Current cycle of the read fragment (must be <= value in settings.seqLength array)
//  uint16_t read_length; // Number of bases for the current read fragment (:= value in settings.seqLength array)
//  uint16_t mate; // Read mate. Normally, there exist only 1 (single-end) or 1 and 2 (paired end). Barcodes are marked as 0.
  std::string root;

  Task() : lane(255), tile(0), seqEl(NULLSEQ), cycle(0), root("") {};

  // constructor initializes all member variables
 Task(uint16_t ln, uint16_t tl, SequenceElement seq, uint16_t cl, std::string rt):
	 lane(ln), tile(tl), seqEl(seq), cycle(cl), root(rt) {};

  // overload << operator for printing
  friend std::ostream& operator<<(std::ostream& os, const Task& t);
};

inline bool operator==(const Task& l, const Task& r){ return (r.lane==l.lane)&&(r.tile==l.tile)&&(r.cycle==l.cycle)&&(r.seqEl==l.seqEl)&&(r.root==l.root); }
inline bool operator!=(const Task& l, const Task& r){ return !(l==r); }

const Task NO_TASK (255,0,NULLSEQ,0,"");

// Task queue data structure. Manages a list of task objects in a thread safe way.
class TaskQueue {
  // the internal queue
  std::queue<Task> tasks;
  
  // mutex to ensure that only one process can access the queue at once
  std::mutex m;

 public:
  // Add element to the task list
  void push(Task t);
  
  // Get element from the task list
  Task pop();

  // return the size of the queue
  uint64_t size();
};


// Agenda item status
typedef uint8_t ItemStatus;
const ItemStatus WAITING = 0;
const ItemStatus BCL_AVAILABLE = 1;
const ItemStatus RUNNING = 2;
const ItemStatus FINISHED = 3;
const ItemStatus RETRY = 4;
const ItemStatus FAILED = 5;
const ItemStatus ERROR = std::numeric_limits<ItemStatus>::max();


// Agenda: monitors the sequencing process and manages the alignment
// - Monitors BCL files
// - generates new tasks
// - receive finished/fail signals
class Agenda {
  // list of items on the agenda. items[lane][tile][cycle]
  std::vector< std::vector< std::vector<ItemStatus> > > items;
  
  // dataset information
  std::string root;
  uint16_t rlen;
  std::vector<uint16_t> lanes;
  std::vector<uint16_t> tiles;
  AlignmentSettings* settings;

 public:
  // initialize agenda with root directory and read length only (all lanes, all tiles)
  Agenda (std::string rt, uint16_t rl, AlignmentSettings* set);

  // initialize agenda with root directory, read length, and lanes (all tiles)
  Agenda (std::string rt, uint16_t rl, std::vector<uint16_t> ln, AlignmentSettings* set);

  // initialize agenda with root directory, read length, lanes, and tiles
  Agenda (std::string rt, uint16_t rl, std::vector<uint16_t> ln, std::vector<uint16_t> tl, AlignmentSettings* set);

  // check for BCL files and update item status
  void update_status();

  // generate a new task from the agenda
  Task get_task();

  // set the status of a task
  void set_status(Task t, ItemStatus status, AlignmentSettings* set);

  // get the status of a task
  ItemStatus get_status(Task t);

  // check if all items of the agenda were processed, if possible
  bool finished();

  // the total number of tasks on the agenda
  uint32_t task_count();

  // the total number of finished tasks on the agenda
  uint32_t tasks_finished();

  // generate a complete TaskQueue with tasks to generate SAM files
  // SAM files can only be generated for tiles where all cycles are completed
  std::vector<Task> get_SAM_tasks();
};



// create a vector with all lane numbers
std::vector<uint16_t> all_lanes();

// create a vector with one lane number
std::vector<uint16_t> one_lane(uint16_t l);

// create a vector with all tile numbers
std::vector<uint16_t> all_tiles();

// create a vector with one tile number
std::vector<uint16_t> one_tile(uint16_t t);



#endif /* PARALLEL_H */
