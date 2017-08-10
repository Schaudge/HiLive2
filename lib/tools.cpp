#include "tools.h"


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

	  auto gaps_vec = globalAlignmentSettings.get_kmer_gaps();
	  if (std::find(gaps_vec.begin(), gaps_vec.end(), i+1) == gaps_vec.end()) {
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
  if (!(it+globalAlignmentSettings.get_kmer_span()-1 < end)) {
    std::cerr << "Error: hash_fw was called using an begin position which had not at least kmer_span bases behind it." << std::endl;
  }
  HashIntoType h = 0;
  std::string::const_iterator last_invalid = it-1;

  h |= twobit_repr(*it);

  std::string::const_iterator kmerEnd = it+globalAlignmentSettings.get_kmer_span();
  ++it;
  int positionInKmer = 2;
  auto kmer_gaps = globalAlignmentSettings.get_kmer_gaps();
  for (; it != kmerEnd; ++it, ++positionInKmer) {
    if (std::find(kmer_gaps.begin(), kmer_gaps.end(), positionInKmer) != kmer_gaps.end())
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
std::string bcl_name(uint16_t ln, uint16_t tl, uint16_t cl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L00" << ln << "/C" << cl << ".1/s_"<< ln <<"_" << tl << ".bcl";
  return path_stream.str();
}

// construct alignment file name from: root, lane, tile, cycle
std::string alignment_name(uint16_t ln, uint16_t tl, uint16_t cl, uint16_t mt){
  std::ostringstream path_stream;
  std::string base = globalAlignmentSettings.get_temp_dir() != "" ? globalAlignmentSettings.get_temp_dir() : globalAlignmentSettings.get_root();
  path_stream << base << "/L00" << ln << "/s_"<< ln << "_" << tl << "." << mt << "."<< cl << ".align";
  return path_stream.str();
}

// construct tile-wise SAM file name from: root, lane, tile
std::string sam_tile_name(uint16_t ln, uint16_t tl, bool write_bam) {
  std::ostringstream path_stream;
  if (write_bam)
    path_stream << globalAlignmentSettings.get_out_dir().string() << "/L00" << ln << "/s_"<< ln << "_" << tl << ".bam";
  else
    path_stream << globalAlignmentSettings.get_out_dir().string() << "/L00" << ln << "/s_"<< ln << "_" << tl << ".sam";
  return path_stream.str();
}

uint16_t getSeqCycle(uint16_t cycle, uint16_t read_number) {
	uint16_t seq_cycle = cycle;
	for ( int i = 0; i < read_number; i++ )
		seq_cycle += globalAlignmentSettings.getSeqById(i).length;
	return seq_cycle;
}

// construct filter file name from: root, lane, tile
std::string filter_name(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L00" << ln << "/s_"<< ln << "_" << tl << ".filter";
  return path_stream.str();
}

// construct position file name from: root, lane, tile
std::string position_name(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "../L00" << ln << "/s_"<< ln << "_" << tl << ".clocs";
  return path_stream.str();
}

std::string get_xml_out_name() {
	std::ostringstream path_stream;
	std::string base = globalAlignmentSettings.get_temp_dir() != "" ? globalAlignmentSettings.get_temp_dir() : globalAlignmentSettings.get_root();
	path_stream << base << "/hilive_settings.xml";
	return path_stream.str();
}
