#include <iostream>
#include <boost/program_options.hpp>

#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/bamout.h"


namespace po = boost::program_options;

std::string license =
"Copyright (c) 2015, Martin S. Lindner, marzin at mail-lindner.de\n"
"All rights reserved.\n"
"\n"
"Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:\n"
"\n"
"1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.\n"
"\n"
"2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.\n"
"\n"
"3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.\n"
"\n"
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.";




// the worker function for the threads
void worker (TaskQueue & tasks, TaskQueue & finished, TaskQueue & failed, BAMOut * bamfile, AlignmentSettings* settings, KixRun* idx, bool & surrender ) {

  // loop that keeps on running until the surrender flag is set
  while ( !surrender ) {
    // try to obtain a new task
    Task t = tasks.pop();
    if ( t != NO_TASK ) {
      // Execute the task
      bool success = true;
      try {
        StreamedAlignment s (t.lane, t.tile, t.root, t.rlen);
        uint64_t num_seeds;
        if (t.cycle < t.rlen) {
            num_seeds = s.extend_alignment(t.cycle,idx,settings);
        }
        else {
            num_seeds = s.extend_last_alignment(t.cycle,idx,settings,bamfile);
        }
        
        std::cout << "Task [" << t << "]: Found " << num_seeds << " seeds." << std::endl;
      }
      catch (const std::exception &e) {
        std::cerr << "Failed to finish task [" << t << "]: " << e.what() << std::endl;
        success = false;
      }
      if (success) {
        finished.push(t);
      }
      else {
        failed.push(t);
      }

    }
    else {
      // send this thread to sleep for a second
      std::this_thread::sleep_for (std::chrono::milliseconds(100));
    }
  }  

}


// create a special SAM worker, that writes out a SAM file for a tile
void sam_worker (TaskQueue & tasks, AlignmentSettings* settings, KixRun* idx) {

  // loop that keeps on running until the surrender flag is set
  while ( true ) {
    // try to obtain a new task
    Task t = tasks.pop();
    if ( t != NO_TASK ) {
      // Execute the task
      alignments_to_sam(t.lane,t.tile,t.root,t.rlen, idx,settings);
    }
    else {
      return;
    }
  }  

}



int main(int argc, char* argv[]) {
  time_t t_start = time(NULL);

  // setting up the command line interface
  po::options_description general("General");
  general.add_options()
    ("help,h", "Print this help message and exit")
    ("license", "Print licensing information and exit");

  po::options_description parameters("Parameters");
  parameters.add_options()
    ("BC_DIR", po::value<std::string>()->required(), "Illumina BaseCalls directory")
    ("INDEX", po::value<std::string>()->required(), "Path to k-mer index")
    ("CYCLES", po::value<CountType>()->required(), "Number of cycles")
    ("OUT", po::value<std::string>()->required(), "Output BAM file name");

  po::options_description io_settings("IO settings");
  io_settings.add_options()
    ("temp", po::value<std::string>(), "Temporary directory for the alignment files [Default: use BaseCalls directory]")
    ("sam,S", po::value<std::string>(), "Create SAM files for each tile. [Default: no SAM files]")
    ("keep-files,k", "Keep intermediate alignment files [Default: false]")
    ("lanes,l", po::value< std::vector<uint16_t> >(), "Select lane [Default: all lanes]")
    ("tiles,t", po::value< std::vector<uint16_t> >(), "Select tile numbers [Default: all tiles]");

  po::options_description alignment("Alignment settings");
  alignment.add_options()
    ("min-errors,e", po::value<CountType>()->default_value(2), "Number of errors tolerated in read alignment")
    ("best-hit,H", "Report only the best alignmnet(s) for each read (default)")
    ("best-n,N", po::value<CountType>()->default_value(2), "Report the N best alignmnets for each read")
    ("all-hits,A", "Report all valid alignments for each read")
    ("disable-ohw-filter", "Disable the One-Hit Wonder filter")
    ("start-ohw", po::value<CountType>()->default_value(K+5), "First cycle to apply One-Hit Wonder filter")
    ("window,w", po::value<DiffType>()->default_value(5), "Set the window size to search for alignment continuation, i.e. maximum insertion/deletion size")
    ("min-quality", po::value<uint16_t>()->default_value(1), "Minimum allowed basecall quality");

  po::options_description technical("Technical settings");
  technical.add_options()
    ("block-size", po::value<uint64_t>(), "Block size for the alignment input/output stream in Bytes. Use -K or -M to specify in Kilobytes or Megabytes")
    (",K", "Interpret the block-size argument as Kilobytes instead of Bytes")
    (",M", "Interpret the block-size argument as Megabytes instead of Bytes")
    ("compression,c", po::value<uint16_t>()->default_value(2), "Compress alignment files. 0: no compression (default) 1: Deflate (smaller) 2: LZ4 (faster)")
    ("num-threads,n", po::value<int>()->default_value(1), "Number of threads to spawn")
    ;

  po::options_description cmdline_options;
  cmdline_options.add(general).add(parameters).add(io_settings).add(alignment).add(technical);

  po::options_description visible_options;
  visible_options.add(general).add(io_settings).add(alignment).add(technical);
  
  std::stringstream help_message;
  help_message << "HiLive v"<< HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR << " - Realtime Alignment of Illumina Reads" << std::endl;
  help_message << "Copyright (c) 2015, Martin S. Lindner" << std::endl;
  help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;
  help_message << "Fixed k-mer size: " << K << std::endl << std::endl;
  help_message << "Usage: " << std::string(argv[0]) << " [options] BC_DIR INDEX CYCLES" << std::endl;
  help_message << "  BC_DIR       Illumina BaseCalls directory of the sequencing run to analyze" << std::endl;
  help_message << "  INDEX        Path to k-mer index file (*.kix)" << std::endl;
  help_message << "  CYCLES       Total number of cycles for read 1" << std::endl;
  help_message << "  OUT          Output BAM file name" << std::endl << std::endl;
  help_message << visible_options << std::endl;


  po::positional_options_description p;
  p.add("BC_DIR", 1);
  p.add("INDEX", 1);
  p.add("CYCLES", 1);
  p.add("OUT", 1);

  po::variables_map vm;
  try {
    // parse arguments
    po::store(po::command_line_parser(argc, argv).
	      options(cmdline_options).positional(p).run(), vm);
    // first check if -h or --help was called
    if (vm.count("help")) {
      std::cout << help_message.str();
      return 1;
    }
    // first check if --license was called
    if (vm.count("license")) {
      std::cout << license << std::endl;
      return 1;
    }
    
    // then check arguments
    po::notify(vm);  
  }
  catch ( boost::program_options::required_option& e ) {
    std::cerr << "Missing Parameter: " << e.what() << std::endl << std::endl;
    std::cout << help_message.str();
    return -1;  
  }
  catch( boost::program_options::error& e) { 
    std::cerr << "Error while parsing command line options: " << e.what() << std::endl << std::endl; 
    std::cout << help_message.str();
    return -1;  
  } 


  // variables to set
  std::string root = vm["BC_DIR"].as<std::string>();
  std::string index_fname = vm["INDEX"].as<std::string>();
  CountType rlen = vm["CYCLES"].as<CountType>();
  std::string outfile = vm["OUT"].as<std::string>();

  std::vector<uint16_t> lanes;
  std::vector<uint16_t> tiles;

  AlignmentSettings settings;
  
  bool write_sam;
  std::string sam_dir;

  if (vm.count("temp"))
    settings.temp_dir = vm["temp"].as<std::string>();
  else
    settings.temp_dir = "";

  if (vm.count("sam")) {
    write_sam = true;
    settings.sam_dir = vm["sam"].as<std::string>();
    if (settings.sam_dir == "") {
      if (settings.temp_dir == "") 
	settings.sam_dir = root;
      else
	settings.sam_dir = settings.temp_dir;
    }
  }
  else
    write_sam = false;


  settings.keep_aln_files = (vm.count("keep-files") > 0);

  if (vm.count("lanes"))
    lanes = vm["lanes"].as< std::vector<uint16_t> >();
  else
    lanes = all_lanes();
 
  if (vm.count("tiles"))
    tiles = vm["tiles"].as< std::vector<uint16_t> >();
  else
    tiles = all_tiles();
  
  if (vm.count("min-errors")) 
    settings.min_errors = vm["min-errors"].as<uint16_t>();
  else
    settings.min_errors = 2;

  if (vm.count("all-hits")) {
    // all hits: disable other modes
    settings.best_hit_mode = false;
    settings.best_n_mode = false;
  } else if (vm.count("best-n")) {
    // enable best-n mode and get parameter
    settings.best_n_mode = true;
    settings.best_n = vm["best-n"].as<CountType>();
  }

  if (vm.count("best-hit")) {
    // enable best-hit mode in all other cases, overwrite other parameters
    settings.best_hit_mode = true;
    settings.best_n_mode = false;
  }
  
  settings.discard_ohw = (vm.count("disable-ohw-filter") == 0);
  
  if (vm.count("start-ohw")) 
    settings.start_ohw = vm["start-ohw"].as<CountType>();
  else
    settings.start_ohw = K+5;

  if (vm.count("window")) 
    settings.window = vm["window"].as<DiffType>();
  else
    settings.window = 5;

  if (vm.count("min-quality")) 
    settings.min_qual = vm["min-quality"].as<uint16_t>();
  else
    settings.min_qual = K+5;

  if (vm.count("block-size")) {
    if (vm.count("-M"))
      settings.block_size = vm["block-size"].as<uint64_t>()*1024*1024;
    else if (vm.count("-K"))
      settings.block_size = vm["block-size"].as<uint64_t>()*1024;
    else
      settings.block_size = vm["block-size"].as<uint64_t>();
  }
  else
    settings.block_size = 64*1024*1024;

  if (vm.count("compression")) 
    settings.compression_format = (uint8_t) vm["compression"].as<uint16_t>();
  else
    settings.compression_format = 2;

  int num_threads = 1;
  if (vm.count("num-threads")) 
    num_threads = vm["num-threads"].as<int>();
  else {
    uint32_t n_cpu = std::thread::hardware_concurrency();
    if (n_cpu > 1){
      num_threads = n_cpu;
    }
  }


  // check paths and file names
  if (!file_exists(index_fname)){
    std::cerr << "Input error: Could not find k-mer index file " << index_fname << std::endl;
    return -1;
  }

  std::size_t found = root.find("BaseCalls");
  if (!(found != std::string::npos && found >= root.size()-10)) {
    std::cerr << "Warning: BaseCalls directory seems to be invalid: " << root << std::endl;
  } 

  if (!is_directory(root)){
    std::cerr << "Input error: Could not find BaseCalls directory " << root << std::endl;
    return -1;
  }

  for ( uint16_t ln : lanes ) {
    std::string ln_dir = root;
    if ( ln < 10 )
      ln_dir += "/L00";
    else if ( ln < 100 )
      ln_dir += "/L0";
    else
      ln_dir += "/L";
    ln_dir += std::to_string(ln);
    if (!is_directory(ln_dir)){
      std::cerr << "Input error: Could not find location of Lane " << ln << ": " << ln_dir << std::endl;
      return -1;
    }
  }
  

  // Report the basic settings
  std::cout << "Running HiLive with " << num_threads << " thread(s)." << std::endl;
  std::cout << "BaseCalls directory:   " << root << std::endl;
  if (settings.temp_dir != "") {
    std::cout << "Temporary directory:   " << settings.temp_dir << std::endl;
  }
  if (write_sam) {
    std::cout << "SAM output directory:  " << settings.sam_dir << std::endl;
  }
  std::cout << "Lanes:                 ";
  for ( uint16_t ln : lanes )
    std::cout << ln << " ";
  std::cout << std::endl;
  std::cout << "K-mer index:           " << index_fname << std::endl;
  std::cout << "Read length:           " << rlen << std::endl;
  std::cout << "Mapping error:         " << settings.min_errors << std::endl;
  if (settings.best_hit_mode) 
    std::cout << "Mapping mode:          Best-Hit-Mode" << std::endl;
  else if (settings.best_n_mode) 
    std::cout << "Mapping mode:          Best-N-Mode" << std::endl;
  else
    std::cout << "Mapping mode:          All-Hits-Mode" << std::endl;
  std::cout << std::endl;


  // load the index
  std::cout << "Loading Index" << std::endl;
  KixRun* index = new KixRun();
  index->deserialize_file(index_fname);


  // Create the overall agenda
  Agenda agenda (root, rlen, lanes, tiles);


  // prepare the alignment
  std::cout << "Initializing Alignment files. Waiting for the first cycle to finish." << std::endl;
  bool first_cycle_available = false;

  // wait for the first cycle to be written. Attention - this loop will wait infinitely long if no first cycle is found
  while ( !first_cycle_available ) {
    // check for new BCL files and update the agenda status
    agenda.update_status();

    // check if the first cycle is available for all tiles
    first_cycle_available = true;
    for ( auto ln : lanes ) {
      for ( auto tl : tiles ) {
	if ( agenda.get_status(Task(ln,tl,1,rlen,"")) != BCL_AVAILABLE) {
	  first_cycle_available = false;
	}
      }
    }
    
    // take a small break
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }

  std::cout << "First cycle complete. Starting alignment." << std::endl;
  
  // write empty alignment file for each tile
  for (uint16_t ln : lanes) {
    for (uint16_t tl : tiles) {
      StreamedAlignment s (ln, tl, root, rlen);
      s.create_directories(&settings);
      s.init_alignment(&settings);
    }
  }

  if (settings.temp_dir != "" && !is_directory(settings.temp_dir)){
    std::cerr << "Error: Could not find temporary directory " << settings.temp_dir << std::endl;
    return -1;
  }

  // Set up the queues
  TaskQueue toDoQ;
  TaskQueue finishedQ;
  TaskQueue failedQ;
  
  // Open the BAM file for output
  BAMOut bamfile (outfile, index);

  // Create the threads
  std::cout << "Creating " << num_threads << " threads." << std::endl;
  bool surrender = false;
  std::vector<std::thread> workers;
  for (int i = 0; i < num_threads; i++) {
    workers.push_back(std::thread(worker, std::ref(toDoQ), std::ref(finishedQ), std::ref(failedQ), &bamfile, &settings, index, std::ref(surrender)));
  }
  
  // Process all tasks on the agenda
  while ( !agenda.finished() ) {
    // check for new BCL files and update the agenda status
    agenda.update_status();
    
    // fill the ToDo queue with tasks from the agenda
    while(true) {
      Task t = agenda.get_task();
      if (t == NO_TASK)
        break;
      toDoQ.push(t);
      agenda.set_status(t,RUNNING);
    }

    // take a look in the finished queue and process finished tasks
    while(true) {
      Task t = finishedQ.pop();
      if (t == NO_TASK)
        break;
      agenda.set_status(t,FINISHED);
    }

    // take a look in the failed queue and process failed tasks
    while(true) {
      Task t = failedQ.pop();
      if (t == NO_TASK)
        break;
      if (agenda.get_status(t) == RUNNING) {
        // give it one more chance
        agenda.set_status(t,RETRY);
        toDoQ.push(t);
      }
      else {
        agenda.set_status(t,FAILED);
        std::cout << "Task failed! " << t << std::endl;
      }
    }

    // take a small break
    std::this_thread::sleep_for (std::chrono::milliseconds(100));
  }  
  
  // Halt the threads
  surrender = true;
  for (auto& w : workers) {
    w.join();
  }
  std::cout << "All threads joined." << std::endl;

  if (write_sam) {
    std::cout << "Writing SAM files." << std::endl;
    // Create individual SAM files for every tile
    TaskQueue sam_tasks; 
    std::vector<Task> tv = agenda.get_SAM_tasks();
    for ( auto t: tv ) {
      sam_tasks.push(t);
    }
    
    workers.clear();
    for (int i = 0; i < num_threads; i++) {
      workers.push_back(std::thread(sam_worker, std::ref(sam_tasks), &settings, index));
    }
    
    for (auto& w : workers) {
      w.join();
    }
  }

  std::cout << "Total run time: " << time(NULL) - t_start << " s" << std::endl;

  delete index;

  return 1;
}
