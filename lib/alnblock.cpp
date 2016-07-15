#include "alnblock.h"


uint16_t DatasetAlignment::get_cycle() {
  return cycle;
}

uint16_t DatasetAlignment::get_lane() {
  return lane;
}

uint16_t DatasetAlignment::get_tile() {
  return tile;
}

std::string DatasetAlignment::get_root() {
  return root;
}
  
CountType DatasetAlignment::get_read_length() {
  return rlen;
}

std::string DatasetAlignment::get_bcl_file(uint16_t cl /* = 0 */) {
  std::ostringstream path_stream;
  uint16_t thiscycle = cl;
  if ( cl == 0 ) {
    // only allow cycle numbers > 0
    assert(cycle > 0);
    thiscycle = cycle;
  }
  path_stream << root << "/L00" << int(lane) << "/C" << int(thiscycle) << ".1/s_"<< int(lane) <<"_" << int(tile) << ".bcl";
  return path_stream.str();
}


std::string DatasetAlignment::get_alignment_file() {
  std::ostringstream path_stream;
  path_stream << root << "/L00" << int(lane) << "/s_"<< int(lane) << "_" << int(tile) << ".align";
  return path_stream.str();
}


std::string DatasetAlignment::get_filter_file() {
  std::ostringstream path_stream;
  path_stream << root << "/L00" << int(lane) << "/s_"<< int(lane) << "_" << int(tile) << ".filter";
  return path_stream.str();
}


std::vector<char> DatasetAlignment::serialize() {
  // calculate total size first
  unsigned long int total_size = 0;

  total_size += sizeof(uint16_t); // lane
  total_size += sizeof(uint16_t); // tile
  total_size += sizeof(CountType); // cycle

  // root directory name
  uint16_t root_size = root.size();
  total_size += sizeof(uint16_t);
  total_size += root_size;

  // number of reads
  unsigned long int num_reads = reads.size();
  total_size += sizeof(uint64_t);

  // read length
  total_size += sizeof(CountType);

  // sum up the sizes of the single alignments
  for (auto it = reads.begin(); it != reads.end(); ++it) {
    total_size += sizeof(uint32_t);
    total_size += (*it).serialize_size();
  }

  // create the vector to store the data
  std::vector<char> data (total_size);
  char* d = data.data();

  // write the lane
  memcpy(d,&lane,sizeof(uint16_t));
  d += sizeof(uint16_t);

  // write the tile
  memcpy(d,&tile,sizeof(uint16_t));
  d += sizeof(uint16_t);

  // write the cycle
  memcpy(d,&cycle,sizeof(CountType));
  d += sizeof(CountType);

  // root directory name
  memcpy(d,&root_size,sizeof(uint16_t));
  d += sizeof(uint16_t);

  memcpy(d,root.c_str(),root_size);
  d += root_size;

  // write the read length
  memcpy(d,&rlen,sizeof(CountType));
  d += sizeof(CountType);

  // write the number of reads
  memcpy(d,&num_reads,sizeof(uint64_t));
  d += sizeof(uint64_t);

  // write the reads
  int i = 0;
  for (auto it = reads.begin(); it != reads.end(); ++it,++i) {

    std::vector<char> sread = it->serialize();
    uint32_t size = sread.size();

    memcpy(d,&size,sizeof(uint32_t));
    d += sizeof(uint32_t);
    memcpy(d,sread.data(),size);
    d += size;
    
  }

  return data;
}


uint64_t DatasetAlignment::serialize_file(std::string f /* = "" */) {
  std::string fname = get_alignment_file();  
  if ( f != "" ) {
    fname = f;
  }

  // serialize data
  std::vector<char> sdata = serialize();

  // write data
  return write_binary_file(fname,sdata);

}



uint64_t DatasetAlignment::deserialize(char* d) {
  // the total number of bytes read
  uint64_t bytes = 0; 

  // read the lane
  memcpy(&lane,d+bytes,sizeof(uint16_t));
  bytes += sizeof(uint16_t);

  // read the tile
  memcpy(&tile,d+bytes,sizeof(uint16_t));
  bytes += sizeof(uint16_t);

  // read the cycle
  memcpy(&cycle,d+bytes,sizeof(unsigned short));
  bytes += sizeof(unsigned short);

  // read the root directory name
  uint16_t root_size;
  memcpy(&root_size,d+bytes,sizeof(uint16_t));
  bytes += sizeof(uint16_t);

  char tmp[root_size+1];
  memcpy(tmp,d+bytes,root_size);
  tmp[root_size] = 0; // make the string null-terminated
  root = tmp;
  bytes += root_size;

  // read the read length
  memcpy(&rlen,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the number of reads
  unsigned long int num_reads;
  memcpy(&num_reads,d+bytes,sizeof(uint64_t));
  bytes += sizeof(uint64_t);

  // reading the reads
  reads.clear();
  reads.reserve(num_reads);
  for (unsigned int i = 0; i < num_reads; ++i) {
    // get the size of the read
    uint32_t size;
    memcpy(&size,d+bytes,sizeof(uint32_t));
    bytes += sizeof(uint32_t);

    // read as a chunk and deserialize
    reads.emplace_back(ReadAlignment(rlen));
    //if (flags[i]) {
    std::vector<char> sread (size);
    memcpy(sread.data(),d+bytes,size);
    reads.back().deserialize(sread.data());
    //}

    bytes += size;
  }

  return bytes;
}


unsigned long int DatasetAlignment::deserialize_file(std::string f /* = "" */) {
  std::string fname = get_alignment_file();
  if ( f != "" ) {
    fname = f;
  }

  // read raw data from file
  std::vector<char> sdata = read_binary_file(fname);
  
  // deserialize data
  deserialize(sdata.data());

  return sdata.size();
}


void DatasetAlignment::initialize(std::string rt, int ln, int tl, int rl) {
  // set the basic stuff
  root = rt;
  lane = ln;
  tile = tl;
  cycle = 0;
  rlen = rl;

  // get the number of reads in this tile by looking in the first bcl file
  std::string fname = get_bcl_file(1);

  // extract the number of reads
  uint32_t num_reads = num_reads_from_bcl(fname);

  // create vector for read alignments
  reads.clear();
  reads.reserve(num_reads);
  for (uint32_t i = 0; i < num_reads; i++){
    reads.emplace_back(ReadAlignment(rlen));
  }
}



void DatasetAlignment::add_nucleotide(KixRun* index, AlignmentSettings* settings) {

  // TODO: check the prerequisites: bcl file must be available

  // ----------------------------------------
  // read the BCL file
  // ----------------------------------------
  std::string fname = get_bcl_file(cycle+1);
  
  std::vector<char> bcl_data = read_binary_file(fname);

  // extract the number of reads from the BCL file
  uint64_t num_reads;
  memcpy(&num_reads,bcl_data.data(),4);

  if (num_reads != reads.size()){
    std::cerr << "Input Error: Number of base calls (" << num_reads << ") does not match the number of reads in this dataset (" << reads.size() << ")." << std::endl;
    return;
  }
  
  // pointer to the beginning of the base call data block
  char* base_calls = bcl_data.data() + 4;

  // ----------------------------------------
  // load the filter flags
  // ----------------------------------------
  std::string filter_fname = get_filter_file();

  std::vector<char> filterdata = read_binary_file(filter_fname);

  // extract the number of reads from the filter file
  unsigned int num_reads_filter;
  memcpy(&num_reads_filter,filterdata.data()+8,4);

  if (num_reads != num_reads_filter){
    std::cerr << "Input Error: Number of reads in filter file (" << num_reads_filter << ") does not match the number of reads in the BCL file (" << num_reads << ")." << std::endl;
    return;
  }

  // pointer to the beginning of the base call data block
  char* filter_flags = filterdata.data() + 12;


  // ----------------------------------------  
  // update the reads
  // ----------------------------------------

  // timing stuff
  std::chrono::high_resolution_clock::duration d_vec = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::high_resolution_clock::duration d_seed = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::high_resolution_clock::duration d_add = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::high_resolution_clock::duration d_rem = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::high_resolution_clock::duration d_sort = std::chrono::high_resolution_clock::duration::zero();


  
  auto rit = reads.begin();
  for (uint64_t i = 0; i < num_reads; ++i, /*++fit,*/ ++rit) {

    // extend the alignment by one nucleotide, if flags are ok
    if (*(filter_flags+i) > 0) {
      std::vector<std::chrono::high_resolution_clock::duration> t = (*rit).extend_alignment(*(base_calls+i), index, settings);
	d_vec += t[0];
	d_seed += t[1];
	d_add += t[2];
	d_rem += t[3];
	d_sort += t[4];
      }else {
      (*rit).seeds.clear();
    }
  }
  
  cycle++;
  
  std::cout << "Time gather seeds:  " << d_vec.count()/1000 << std::endl;
  std::cout << "Time extend seeds:  " << d_seed.count()/1000 << std::endl;
  std::cout << "Time add new seeds: " << d_add.count()/1000 << std::endl;
  std::cout << "Time filter seeds:  " << d_rem.count()/1000 << std::endl;
  std::cout << "Time sorting seeds: " << d_sort.count()/1000 << std::endl;
  

  return;
}



void DatasetAlignment::report_alignments_spam(std::string fname, KixRun* idx) {

  // open output file for writing
  FILE* ofile;
  ofile = fopen(fname.c_str(), "w");
  if (!ofile) {
    std::cerr << "Error opening output file " << fname << std::endl;
    return;
  }

  // write out all seeds of each read
  uint64_t rd = 0;
  for (auto rit = reads.begin(); rit != reads.end(); ++rit, ++rd) {
    // compose a read name
    std::string rname = std::string("read_") + std::to_string(rd);

    for (auto sit = rit->seeds.begin(); sit != rit->seeds.end(); ++sit) {
      PositionType pos = (*sit)->start_pos;
      if (pos < 0) {
	pos = -pos - rlen + K;
      }
      std::stringstream out;
      out << rname << "\t" << idx->get_name((*sit)->gid)  << "\t" << pos  << "\t" << CIGAR(**sit) << std::endl;
      fwrite(out.str().c_str(), 1, out.str().size(), ofile);
    }

  }

  fclose(ofile);
}


