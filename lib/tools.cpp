#include "tools.h"

// reads a binary file from hdd and stores its raw content in a char vector
std::vector<char> read_binary_file(const std::string &fname) {
  
  // get file size
  uint64_t size = get_filesize(fname);

  // open binary file
  FILE* f;
  f = fopen(fname.c_str(), "rb");

  if (!f) {
    std::cerr << "Error reading binary file " << fname << ": Could not open file." << std::endl;
    return std::vector<char>();
  }

  // allocate memory
  std::vector<char> data (size);
  
  // read all data at once
  uint64_t read = fread(data.data(), 1, size, f);

  if (read != size){
    std::cerr << "Error reading binary file " << fname << ": File size: " << size << " bytes. Read: " << read << " bytes." << std::endl;
    return std::vector<char>();
  }
  
  fclose(f);

  return data;
}


// checks if a directory with the given path exists
bool is_directory(const std::string &path) {
  if ( boost::filesystem::exists(path) ) {
    if ( boost::filesystem::is_directory(path) ) {
      return true;
    }
    else {
      return false;
    }
  }
  else { 
    return false;
  }
}

// checks if a file exists
bool file_exists(const std::string &fname) {
  return boost::filesystem::exists(fname);
  /*if (FILE *f = fopen(fname.c_str(), "r")) {
    fclose(f);
    return true;
  }
  else {
    return false;
    }*/ 
}


// writes a char vector into a binary file
uint64_t write_binary_file(const std::string &fname, const std::vector<char> & data) {

  // open binary file
  FILE* ofile;
  ofile = fopen(fname.c_str(), "wb");

  if (!ofile) {
    std::cerr << "Error serializing object to file " << fname << ": Could not open file for writing." << std::endl;
    return 0;
  }

  // write all data
  uint64_t written = fwrite(data.data(), 1, data.size(), ofile);
  
  // close file
  fclose(ofile);

  if (written != data.size()){
    std::cerr << "Error serializing object to file " << fname << ": Total size: " << data.size() << " bytes. Written: " << written << " bytes." << std::endl;
  }

  return written;
}



// extract the number of reads from a BCL file
uint32_t num_reads_from_bcl(std::string bcl) {
  // open BCL file of first cycle
  FILE* ifile;
  ifile = fopen(bcl.c_str(), "rb");

  if (!ifile) {
    std::cerr << "Error reading BCL file " << bcl << ": Could not open file." << std::endl;
    return 0;
  }

  // extract the number of reads
  uint32_t num_reads;
  assert(fread(&num_reads, 1, sizeof(uint32_t), ifile));

  // close file
  fclose (ifile);
  
  return num_reads;
}



// Update an existing kmer by left-shifting all nucleotides and appending new nucleotide
void update_kmer(HashIntoType &kmer, HashIntoType nuc) {
  
  // left-shift the previous k-mer
  kmer = kmer << 2;
  
  // 'or' in the current nucleotide
  // only use the last 2 bits
  kmer |= (nuc & 3);
  
  // mask off the 2 bits we shifted over.
  kmer &= MASK;
  
}


/* get the size (in bytes) of a file */
std::ifstream::pos_type get_filesize(const std::string &fname)
{
  std::ifstream in(fname, std::ios::binary | std::ios::ate);
  return in.tellg(); 
}


/* calculates the first forward and reverse complement k-mer in the 
   string <kmer> and returns the canonical representation. */
HashIntoType hash(const char * kmer, HashIntoType& _h, HashIntoType& _r)
{
  assert(strlen(kmer) >= K);

  HashIntoType h = 0, r = 0;

  h |= twobit_repr(kmer[0]);
  r |= twobit_comp(kmer[K-1]);

  for (unsigned int i = 1, j = K - 2; i < K; i++, j--) {
    h = h << 2;
    r = r << 2;

    h |= twobit_repr(kmer[i]);
    r |= twobit_comp(kmer[j]);
  }

  _h = h;
  _r = r;

  return (h)<(r)?h:r;
}

/* calculates the first forward k-mer in the string <kmer> */
HashIntoType hash_fw(const char * kmer, HashIntoType& _h)
{
  assert(strlen(kmer) >= K);

  HashIntoType h = 0;
  HashIntoType last_invalid = K+1;

  h |= twobit_repr(kmer[0]);

  for (unsigned int i = 1; i < K; i++) {
    h = h << 2;
    h |= twobit_repr(kmer[i]);
    if ( seq_chars.find(kmer[i]) == std::string::npos ) {
      last_invalid = K+1-i;
    }
  }

  _h = h;

  return last_invalid;
}


/* returns the reverse complement of a k-mer */
HashIntoType rc(HashIntoType fw)
{

  HashIntoType rc = 0;

  // Illumina uses
  // A = 0b00
  // C = 0b01
  // G = 0b10
  // T = 0b11
  // so, inverting bits produces the rc: rc(A) = ~A

  for (unsigned int i = 0; i < 2*K; i+=2) {
    rc |= (~(fw >> i) & 3) << (2*K - i - 2);
  }

  return rc;
}



// file name construction functions

// construct BCL file name from: root, lane, tile, cycle
std::string bcl_name(std::string rt, uint16_t ln, uint16_t tl, uint16_t cl) {
  std::ostringstream path_stream;
  path_stream << rt << "/L00" << ln << "/C" << cl << ".1/s_"<< ln <<"_" << tl << ".bcl";
  return path_stream.str();
}


// construct alignment file name from: root, lane, tile, cycle
std::string alignment_name(std::string rt, uint16_t ln, uint16_t tl, uint16_t cl) {
  std::ostringstream path_stream;
  path_stream << rt << "/L00" << ln << "/s_"<< ln << "_" << tl << "." << cl << ".align";
  return path_stream.str();
}

// construct tile-wise SAM file name from: root, lane, tile
std::string sam_tile_name(std::string rt, uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << rt << "/L00" << ln << "/s_"<< ln << "_" << tl << ".sam";
  return path_stream.str();
}

// construct lane-wise SAM file name from: root, lane
std::string sam_lane_name(std::string rt, uint16_t ln) {
  std::ostringstream path_stream;
  path_stream << rt << "/L00" << ln << "/s_"<< ln << ".sam";
  return path_stream.str();
}

// construct filter file name from: root, lane, tile
std::string filter_name(std::string rt, uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << rt << "/L00" << ln << "/s_"<< ln << "_" << tl << ".filter";
  return path_stream.str();
}

// construct position file name from: root, lane, tile
std::string position_name(std::string rt, uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << rt << "../L00" << ln << "/s_"<< ln << "_" << tl << ".clocs";
  return path_stream.str();
}


