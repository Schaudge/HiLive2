#include "tools.h"

// compares two genome positions by position (not by genome id)
bool gp_compare (GenomePosType i,GenomePosType j) { 
  return (i.pos < j.pos); 
}

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
  bool res = fread(&num_reads, 1, sizeof(uint32_t), ifile);
  assert(res); // fread(&num_reads, 1, sizeof(uint32_t), ifile)

  // close file
  fclose (ifile);
  
  return num_reads;
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
  assert(strlen(kmer) >= globalAlignmentSettings.get_kmer_span());

  HashIntoType h = 0, r = 0;

  h |= twobit_repr(kmer[0]);
  r |= twobit_comp(kmer[globalAlignmentSettings.get_kmer_span()-1]);

  for (unsigned int i = 1, j = globalAlignmentSettings.get_kmer_span()-2; i < globalAlignmentSettings.get_kmer_span(); i++, j--) {
    // if i not gap position
    if (std::find(globalAlignmentSettings.get_kmer_gaps().begin(), globalAlignmentSettings.get_kmer_gaps().end(), i+1) == globalAlignmentSettings.get_kmer_gaps().end()) {
      h = h << 2;
      h |= twobit_repr(kmer[i]);
      r = r << 2;
      r |= twobit_comp(kmer[j]);
    }
  }

  _h = h;
  _r = r;

  return (h)<(r)?h:r;
}

/* calculates the first forward k-mer in the string <kmer> */
std::string::const_iterator hash_fw(std::string::const_iterator it, std::string::const_iterator end, HashIntoType& _h)
{
  assert(it+globalAlignmentSettings.get_kmer_span()-1 < end);
  HashIntoType h = 0;
  std::string::const_iterator last_invalid = it-1;

  h |= twobit_repr(*it);

  std::string::const_iterator kmerEnd = it+globalAlignmentSettings.get_kmer_span();
  ++it;
  int positionInKmer = 2;
  for (; it != kmerEnd; ++it, ++positionInKmer) {
    if (std::find(globalAlignmentSettings.get_kmer_gaps().begin(), globalAlignmentSettings.get_kmer_gaps().end(), positionInKmer) != globalAlignmentSettings.get_kmer_gaps().end())
        continue;
    h = h << 2;
    h |= twobit_repr(*it);
    if ( seq_chars.find(*it) == std::string::npos ) {
      last_invalid = it+globalAlignmentSettings.get_kmer_span()-1;
    }
  }

  _h = h;
  return last_invalid;
}


/* returns the sequence of a k-mer */
std::string unhash(HashIntoType myHash, unsigned hashLen)
{
	std::string kmer = "";

	unsigned mask = 3;
	for (unsigned i = 1; i<pow(2,2*hashLen); i *= 4) {
		kmer.push_back(revtwobit_repr(myHash & mask));
		myHash = myHash >> 2;
	}
	std::reverse(kmer.begin(), kmer.end());
	return kmer;
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
std::string sam_tile_name(std::string rt, uint16_t ln, uint16_t tl, bool write_bam) {
  std::ostringstream path_stream;
  if (write_bam)
    path_stream << rt << "/L00" << ln << "/s_"<< ln << "_" << tl << ".bam";
  else
    path_stream << rt << "/L00" << ln << "/s_"<< ln << "_" << tl << ".sam";
  return path_stream.str();
}

// construct lane-wise SAM file name from: root, lane
std::string sam_lane_name(std::string rt, uint16_t ln, bool write_bam) {
  std::ostringstream path_stream;
  if (write_bam)
    path_stream << rt << "/L00" << ln << "/s_"<< ln << ".bam";
  else
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
