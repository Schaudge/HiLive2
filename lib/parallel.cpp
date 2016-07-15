#include "parallel.h"



std::ostream& operator<<(std::ostream& os, const Task& t)
{
  os << "Lane " << t.lane << " Tile " << t.tile << " Cycle " << t.cycle;
  return os;
}


// Add element to the task list
void TaskQueue::push(Task t) {
  std::lock_guard<std::mutex> lk(m);
  tasks.push(t);
}
  
// Get element from the task list. If TaskList is empty, return NO_TASK.
Task TaskQueue::pop() {
  std::lock_guard<std::mutex> lk(m);
  if (!tasks.empty()) {
    Task t = tasks.front();
    tasks.pop();
    return t;
  }
  else {
    return NO_TASK;
  }
}

// return the size of the queue
uint64_t TaskQueue::size() {
  std::lock_guard<std::mutex> lk(m);
  return tasks.size();
}



// create a vector with all lane numbers
std::vector<uint16_t> all_lanes() {
  std::vector<uint16_t> ln;
  for (uint16_t l=0; l < 8; l++)
    ln.push_back(l+1);
  return ln;
}

// create a vector with one lane number
std::vector<uint16_t> one_lane(uint16_t l) {
  return std::vector<uint16_t> (1,l);
}


// create a vector with all tile numbers
std::vector<uint16_t> all_tiles() {
  std::vector<uint16_t> tl;
  for (uint16_t l = 0; l < 2; l++) {
    for (uint16_t s = 0; s < 3; s++) {
      for (uint16_t t = 0; t < 16; t++) {
	// construct tile number
	tl.push_back( (l+1)*1000 + (s+1)*100 + (t+1) );
      }      
    }
  }
  return tl;
}

// create a vector with one tile number
std::vector<uint16_t> one_tile(uint16_t t) {
  return std::vector<uint16_t> (1,t);
}


// initialize agenda with root directory and read length only (all lanes, all tiles)
Agenda::Agenda (std::string rt, uint16_t rl) {
  
  // add lanes 1-8 to the list
  std::vector<uint16_t> ln = all_lanes();
  
  // call the tiles constructor
  Agenda(rt, rl, ln);

}

// initialize agenda with root directory, read length, and lanes (all tiles)
Agenda::Agenda (std::string rt, uint16_t rl, std::vector<uint16_t> ln) {

  // add all tiles to the list
  std::vector<uint16_t> tl = all_tiles();
  
  // call the full constructor
  Agenda (rt, rl, ln, tl);

}

// initialize agenda with root directory, read length, lanes, and tiles
Agenda::Agenda (std::string rt, uint16_t rl, std::vector<uint16_t> ln, std::vector<uint16_t> tl) {

  root = rt;
  rlen = rl;
  lanes = ln;
  tiles = tl;

  // set up the agenda
  items.clear();
  for (uint16_t ln_id = 0; ln_id < lanes.size(); ln_id++) {
    std::vector<std::vector<ItemStatus> > lane_status;
    for (uint16_t tl_id = 0; tl_id < tiles.size(); tl_id++) {
      std::vector<ItemStatus> tile_status (rlen, WAITING);
      lane_status.push_back(tile_status);
    } 
    items.push_back(lane_status);    
  }

}

// check for BCL files and update item status
void Agenda::update_status () {

  // iterate over lanes
  for (uint16_t ln_id = 0; ln_id < items.size(); ++ln_id) {

    // iterate over all tiles
    for (uint16_t tl_id = 0; tl_id < items[ln_id].size(); ++tl_id) {

      // get the first cycle that is not in the FINISHED status
      uint16_t first_unfinished = 0;
      while ( (first_unfinished < items[ln_id][tl_id].size()) && (items[ln_id][tl_id][first_unfinished] == FINISHED)) {
	first_unfinished++;
      }

      // if there is one, check if there is a BCL file available
      if ((first_unfinished != items[ln_id][tl_id].size()) && (items[ln_id][tl_id][first_unfinished] == WAITING)) {
	std::string this_fname = bcl_name(root, lanes[ln_id], tiles[tl_id], first_unfinished+1);
	// only change the status if the file exists
	if ( file_exists(this_fname) ) {
	  // TODO: probably find a way to check if the machine currently writes to that file
	  items[ln_id][tl_id][first_unfinished] = BCL_AVAILABLE;
	}
      }

    }

  }
  
}


// generate a new task from the agenda
Task Agenda::get_task(){
  // iterate over lanes
  for (uint16_t ln_id = 0; ln_id < items.size(); ++ln_id) {

    // iterate over all tiles
    for (uint16_t tl_id = 0; tl_id < items[ln_id].size(); ++tl_id) {

      // check if there is a cycle with an unprocessed BCL file
      uint16_t unprocessed = 0;
      while ( (unprocessed < items[ln_id][tl_id].size()) && (items[ln_id][tl_id][unprocessed] != BCL_AVAILABLE)) {
	unprocessed++;
      }

      // generate a new task if there is an unprocessed BCL file
      if ( unprocessed != items[ln_id][tl_id].size() ) {
	Task t (lanes[ln_id], tiles[tl_id], unprocessed+1, rlen, root);
	return t;
      }

    }

  }
  // return indicator that no new task could be created
  return NO_TASK;
}


// set a status
void Agenda::set_status(Task t, ItemStatus status) {
  // get the lane index
  uint64_t diff = std::find(lanes.begin(), lanes.end(), t.lane) - lanes.begin();
  if ( diff >= lanes.size() ) {
    throw std::out_of_range("Lane ID out of range.");
  }
  uint16_t ln_id = diff;
  
  // get the tile index
  diff = std::find(tiles.begin(), tiles.end(), t.tile) - tiles.begin();
  if ( diff >= tiles.size() ) {
    throw std::out_of_range("Tile ID out of range.");
  }
  uint16_t tl_id = diff;

  // get the cycle index
  if ( (t.cycle > rlen) || (t.cycle == 0) ) {
    throw std::out_of_range("Cycle out of range.");
  }
  uint16_t cl_id = t.cycle -1;

  items[ln_id][tl_id][cl_id] = status;
}


// get the status of a task
ItemStatus Agenda::get_status(Task t) {
  // get the lane index
  uint64_t diff = std::find(lanes.begin(), lanes.end(), t.lane) - lanes.begin();
  if ( diff >= lanes.size() ) {
    throw std::out_of_range("Lane ID out of range.");
  }
  uint16_t ln_id = diff;
  
  // get the tile index
  diff = std::find(tiles.begin(), tiles.end(), t.tile) - tiles.begin();
  if ( diff >= tiles.size() ) {
    throw std::out_of_range("Tile ID out of range.");
  }
  uint16_t tl_id = diff;

  // get the cycle index
  if ( (t.cycle > rlen) || (t.cycle == 0) ) {
    throw std::out_of_range("Cycle out of range.");
  }
  uint16_t cl_id = t.cycle -1;

  return items[ln_id][tl_id][cl_id];
}


// check if all items of the agenda were processed, if possible
bool Agenda::finished() {
  // check for each tile if either all cycles are finished OR there is a failed status item
  for (uint16_t ln_id = 0; ln_id < items.size(); ++ln_id) {
    for (uint16_t tl_id = 0; tl_id < items[ln_id].size(); ++tl_id) {
      for (uint16_t cl_id = 0; cl_id < items[ln_id][tl_id].size(); ++cl_id) {
	ItemStatus s = items[ln_id][tl_id][cl_id];
	if ( s == FAILED ) {
	  // the rest of the tile is "allowed" to be unprocessed --> skip
	  continue;
	}
	else if (s != FINISHED) {
	  // otherwise any other status means that the agenda is not finished
	  return false;
	}
      }
    }
  }
  return true;
}


// the total number of tasks on the agenda
uint32_t Agenda::task_count() {
  return lanes.size() * tiles.size() * rlen;
}


// the total number of finished tasks on the agenda
uint32_t Agenda::tasks_finished() {
  uint32_t num_finished = 0;
  // iterate over all items and count the finished tasks
  for (uint16_t ln_id = 0; ln_id < items.size(); ++ln_id) {
    for (uint16_t tl_id = 0; tl_id < items[ln_id].size(); ++tl_id) {
      for (uint16_t cl_id = 0; cl_id < items[ln_id][tl_id].size(); ++cl_id) {
	if (items[ln_id][tl_id][cl_id] == FINISHED) {
	  num_finished++;
	}
      }
    }
  }
  return num_finished;
}


// generate a complete TaskQueue with tasks to generate SAM files
// SAM files can only be generated for tiles where all cycles are completed
std::vector<Task> Agenda::get_SAM_tasks() {
  std::vector<Task> tv;
  // find tiles that are completely mapped
  for (uint16_t ln_id = 0; ln_id < items.size(); ++ln_id) {
    for (uint16_t tl_id = 0; tl_id < items[ln_id].size(); ++tl_id) {
      if ( items[ln_id][tl_id][rlen-1] == FINISHED ) {
	tv.push_back(Task(lanes[ln_id],tiles[tl_id],rlen,rlen,root));
      }
    }
  }
  
  return tv;
}
