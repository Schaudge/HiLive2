#include "tools.h"


///////////////////////////////////
////////// K-mer Hashing //////////
///////////////////////////////////

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

std::string unhash(HashIntoType myHash, unsigned hashLen)
{
	std::string kmer = "";

	unsigned mask = 3;
	for (unsigned i = 1; i<pow(2,2*hashLen); i *= 4) {
		kmer.push_back(revtwobit_repr(myHash & mask));
		myHash = myHash >> 4;
	}
	std::reverse(kmer.begin(), kmer.end());
	return kmer;
}



////////////////////////////////////////////
////////// File name construction //////////
////////////////////////////////////////////

std::string bcl_name(uint16_t ln, uint16_t tl, uint16_t cl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L00" << ln << "/C" << cl << ".1/s_"<< ln <<"_" << tl << ".bcl";
  return path_stream.str();
}

std::string alignment_name(uint16_t ln, uint16_t tl, uint16_t cl, uint16_t mt){
  std::ostringstream path_stream;
  std::string base = globalAlignmentSettings.get_temp_dir() != "" ? globalAlignmentSettings.get_temp_dir() : globalAlignmentSettings.get_root();
  path_stream << base << "/L00" << ln << "/s_"<< ln << "_" << tl << "." << mt << "."<< cl << ".align";
  return path_stream.str();
}

uint16_t getSeqCycle(uint16_t cycle, uint16_t seq_id) {
	uint16_t seq_cycle = cycle;
	for ( int i = 0; i < seq_id; i++ )
		seq_cycle += globalAlignmentSettings.getSeqById(i).length;
	return seq_cycle;
}

uint16_t getMateCycle( uint16_t mate_number, uint16_t seq_cycle ) {

	// Invalid mate
	if ( mate_number == 0 || mate_number > globalAlignmentSettings.get_mates() )
		return 0;

	// Iterate through all sequence elements (including barcodes)
	for ( CountType id = 0; id < globalAlignmentSettings.get_seqs().size(); id++ ) {

		// Current sequence element
		SequenceElement seq = globalAlignmentSettings.getSeqById(id);

		// Seq is mate of interest
		if ( seq.mate == mate_number )
			return ( seq.length > seq_cycle ? seq_cycle : seq.length );

		// Not enough cycles left to reach mate of interest
		else if ( seq.length >= seq_cycle )
			return 0;

		// Reduce number of cycles by the Seq length
		else
			seq_cycle -= seq.length;

	}

	// Should not be reached
	return 0;
}


std::string filter_name(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L00" << ln << "/s_"<< ln << "_" << tl << ".filter";
  return path_stream.str();
}


std::string position_name(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "../L00" << ln << "/s_"<< ln << "_" << tl << ".clocs";
  return path_stream.str();
}


std::string get_settings_name() {
	std::ostringstream path_stream;
	std::string base = globalAlignmentSettings.get_temp_dir() != "" ? globalAlignmentSettings.get_temp_dir() : globalAlignmentSettings.get_root();
	path_stream << base << "/hilive_settings.xml";
	return path_stream.str();
}


std::string get_out_log_name() {
	return ( globalAlignmentSettings.get_out_dir() + "/hilive_out.log" );
}



////////////////////////////////////
////////// SAM/BAM output //////////
////////////////////////////////////

seqan::BamHeader getBamHeader() {
	std::stringstream ss;
	ss.str(std::string());
	ss << HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR;

	seqan::BamHeader header;
	resize(header, 2);

	// @HD header.
	seqan::resize(header[0].tags, 2);
	header[0].type = seqan::BAM_HEADER_FIRST;
	header[0].tags[0].i1 = "VN";
	header[0].tags[0].i2 = "1.5";
	header[0].tags[1].i1 = "GO";
	header[0].tags[1].i2 = "query";

	// @PG header.
	seqan::resize(header[1].tags, 3);
	header[1].type = seqan::BAM_HEADER_PROGRAM;
	header[1].tags[0].i1 = "ID";
	header[1].tags[0].i2 = "hilive";
	header[1].tags[1].i1 = "PN";
	header[1].tags[1].i2 = "HiLive";
	header[1].tags[2].i1 = "VN";
	header[1].tags[2].i2 = ss.str();

	return header;
}


std::string getBamTempFileName(std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_write_bam() ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_out_dir() << "/hilive_out_" << "cycle" << std::to_string(cycle) << "_" << barcode << ".temp" << file_suffix;
	return fname.str();
}


std::string getBamFileName(std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_write_bam() ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_out_dir() << "/hilive_out_" << "cycle" << std::to_string(cycle) << "_" << barcode << file_suffix;
	return fname.str();
}
