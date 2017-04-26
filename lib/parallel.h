#ifndef PARALLEL_H
#define PARALLEL_H

#include "headers.h"
#include "definitions.h"
#include "tools.h"
#include "kindex.h"

//------ Threading tools --------------------------------------------//

/**
 * Task data structure. Contains all information for a thread to process a BCL file.
 * @author Martin Lindner
 */
struct Task {
  /** The lane of the task. */
  uint16_t lane;
  /** The tile of the task. */
  uint16_t tile;
  /** Struct containing the read properties (Barcode vs. sequence; length; mate). */
  SequenceElement seqEl;
  /** Current cycle of the particular read (in general, this does NOT equal the sequencing cycle!). Must be <=seqEl.length. */
  uint16_t cycle;
  /** Base call root directory */
  std::string root;

  /**
   * Constructor for a NULL task.
   * @author Tobias Loka
   */
  Task() : lane(255), tile(0), seqEl(NULLSEQ), cycle(0), root("") {};

  /**
   * Constructor for a valid task.
   * @param ln The lane number.
   * @param tl The tile number.
   * @param seq The respective seqEl element for the current read containing information about length, type (barcode vs. sequence), mate number ...
   * @param cl The cycle of the current read (in general, this does NOT equal the sequencing cycle!). Must be <=seqEl.length.
   * @param rt Base call root directory.
   * @author Martin Lindner
   */
 Task(uint16_t ln, uint16_t tl, SequenceElement seq, uint16_t cl, std::string rt):
	 lane(ln), tile(tl), seqEl(seq), cycle(cl), root(rt) {};

  /**
   * Overload of the << operator. Defines the cout form of a task.
   * @author Martin Lindner
   */
  friend std::ostream& operator<<(std::ostream& os, const Task& t);
};

/**
 * Overload of the == operator.
 * @return true, if all fields/variables of the compared tasks equal.
 * @author Martin Lindner
 */
inline bool operator==(const Task& l, const Task& r){ return (r.lane==l.lane)&&(r.tile==l.tile)&&(r.cycle==l.cycle)&&(r.seqEl==l.seqEl)&&(r.root==l.root); }

/**
 * Overload of the != operator.
 * @return true, if at least one field/variable of the compared tasks is different.
 * @author Martin Lindner
 */
inline bool operator!=(const Task& l, const Task& r){ return !(l==r); }

/**
 * Definition of a NULL task.
 * @author Martin Lindner
 */
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

 public:
  // initialize agenda with root directory and read length only (all lanes, all tiles)
  Agenda (std::string rt, uint16_t rl);

  // initialize agenda with root directory, read length, and lanes (all tiles)
  Agenda (std::string rt, uint16_t rl, std::vector<uint16_t> ln);

  // initialize agenda with root directory, read length, lanes, and tiles
  Agenda (std::string rt, uint16_t rl, std::vector<uint16_t> ln, std::vector<uint16_t> tl);

  // check for BCL files and update item status
  void update_status();

  // generate a new task from the agenda
  Task get_task();

  // set the status of a task
  void set_status(Task t, ItemStatus status);

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
