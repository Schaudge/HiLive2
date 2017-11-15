#include "parallel.h"

std::ostream& operator<<(std::ostream& os, const Task& t)
{
  std::string mate = t.seqEl.mate == 0 ? "b" : std::to_string(t.seqEl.mate);
  os << "Lane " << t.lane << " Tile " << t.tile << " Cycle " << mate << "." << t.cycle;
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

// initialize agenda with read length only (all lanes, all tiles)
Agenda::Agenda (uint16_t rl) {
  
  // add lanes 1-8 to the list
  std::vector<uint16_t> ln = all_lanes();
  
  // call the tiles constructor
  Agenda(rl, ln);

}

// initialize agenda with read length and lanes (all tiles)
Agenda::Agenda (uint16_t rl, std::vector<uint16_t> ln) {

  // add all tiles to the list
  std::vector<uint16_t> tl = all_tiles();
  
  // call the full constructor
  Agenda (rl, ln, tl);

}

// initialize agenda with read length, lanes, and tiles
Agenda::Agenda (uint16_t rl, std::vector<uint16_t> ln, std::vector<uint16_t> tl) {

	Agenda(rl, ln, tl, 1);

}

// initialize agenda with read length, lanes, and tiles
Agenda::Agenda (uint16_t rl, std::vector<uint16_t> ln, std::vector<uint16_t> tl, CountType start_cycle) {

  rlen = rl;
  lanes = ln;
  tiles = tl;

  // set up the agenda
  items.clear();
  for (uint16_t ln_id = 0; ln_id < lanes.size(); ln_id++) {
    std::vector<std::vector<ItemStatus> > lane_status;
    for (uint16_t tl_id = 0; tl_id < tiles.size(); tl_id++) {

    	// Status for finished cycles if "--continue" was used.
    	std::vector<ItemStatus> tile_status (start_cycle-1, FINISHED);

    	// Waiting cycles
    	std::vector<ItemStatus> waiting_status (rlen-(start_cycle-1), WAITING);

    	// Merge vectors
    	tile_status.insert(tile_status.end(), waiting_status.begin(), waiting_status.end());

    	// Push back
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

//			std::cout << ln_id << ";" << tl_id << ";" << first_unfinished << std::endl;

			// if there is one, check if there is a BCL file available
			if ((first_unfinished != items[ln_id][tl_id].size()) && (items[ln_id][tl_id][first_unfinished] == WAITING)) {
				std::string this_fname = bcl_name(lanes[ln_id], tiles[tl_id], first_unfinished+1);
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
				uint16_t cycle = unprocessed + 1;
				uint16_t read_no = 0;
				while ( cycle > globalAlignmentSettings.getSeqById(read_no).length) {
					cycle -= globalAlignmentSettings.getSeqById(read_no).length;
					read_no += 1;
				}
				Task t (lanes[ln_id], tiles[tl_id], globalAlignmentSettings.getSeqById(read_no), cycle);
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
  uint16_t cl_id = getSeqCycle(t.cycle,t.seqEl.id) -1;

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

// check if all items of the agenda were processed, if possible
bool Agenda::finished( CountType cycle ) {
	// check for each tile if either all cycles are finished OR there is a failed status item
	for (uint16_t ln_id = 0; ln_id < items.size(); ++ln_id) {
		for (uint16_t tl_id = 0; tl_id < items[ln_id].size(); ++tl_id) {
			for (uint16_t cl_id = 0; cl_id < cycle; ++cl_id) {
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

bool Agenda::cycle_available( CountType cycle ) {

	if ( cycle == 0 || cycle > rlen )
		return false;

	for (uint16_t ln_id = 0; ln_id < items.size(); ++ln_id) {
		for (uint16_t tl_id = 0; tl_id < items[ln_id].size(); ++tl_id) {
			if ( items[ln_id][tl_id][cycle-1] == WAITING )
				return false;
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
