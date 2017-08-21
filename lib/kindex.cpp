#include "kindex.h"

//////////////////////////////
////////// KixBuild //////////
//////////////////////////////

int KixBuild::create_index( const std::string &fname, const std::string &iname, bool convert_spaces, bool trim_ids ) {

	// Add all sequences of a (multi-)fasta file to the set of sequences.
	try {
	if ( add_fasta(fname, convert_spaces, trim_ids) == 0 )
		return 1;
	} catch ( std::ifstream::failure &readErr ) {
		std::cerr << readErr.what() << std::endl;
		return 2;
	}

	// save the metadata
	if ( ! save_metadata(iname) ) {
		std::cerr << "Writing metadata failed." << std::endl;
		return 3;
	}

	// modify sequences for index building (reverse and add reverse complement)
	seqan::StringSet<seqan::DnaString> all_seqs;

	for(uint16_t i = 0; i < length(seqs); i++) {
		seqan::ModifiedString<seqan::DnaString, seqan::ModReverse> rev(seqs[i]);
		seqan::appendValue(all_seqs, rev);
		seqan::complement(seqs[i]);
		seqan::appendValue(all_seqs, seqs[i]);
	}

	// Define the index
	FMIndex idx(all_seqs);

	// Force the index to be created
	if( ! seqan::indexCreate(idx) )
		return 4;

	// Save the index to file
	if ( save_fmindex(idx, iname) != 1 )
		return 5;

	return 0;
}

int KixBuild::add_fasta(const std::string &fname, bool convert_spaces, bool trim_ids) {

  std::ios::sync_with_stdio(false);
  std::ifstream::sync_with_stdio(false);

  // open input fasta file
  std::ifstream infile (fname.c_str());
  assert(infile.is_open());

  std::string line;	// current line of fasta file
  std::string seq_name;	// name of the current sequence
  bool startNewSequence = false;	// true if last line was a header / new sequence begins
  seqan::DnaString newSeq = "";	// current DNA sequence

  while(getline(infile, line)) {

	  // ignore empty lines
	  if (line.length() == 0) {continue;};

    
	  // check for Windows newline characters (just in case the fasta comes from Windows --> thanks Simon)
	  if (line[line.length()-1] == '\r'){
		  line.erase(line.length()-1);
	  }

	  // handle header line
	  if (line[0] == '>') {

		  // add last sequence to sequence vector (if exist)
		  if ( newSeq != "" ) {
			  seqan::appendValue(seqs, newSeq);
		  }

		  // trim the header (remove leading and trailing whitespaces)
		  trim(line);

		  // handle sequence name options
		  if (convert_spaces)
			  std::replace( line.begin(), line.end(), ' ', '_');

		  if (trim_ids)
			  seq_name = line.substr(1,line.find(' ')-1);
		  else
			  seq_name = line.substr(1,line.length()-1);

		  // add the sequence name to the sequence name vector
		  seq_names.push_back(seq_name);

		  // init new sequence fields
		  seq_lengths.push_back(0);
		  startNewSequence = true;
		  newSeq = "";
	  }

	  // handle sequence lines
	  else {
		  // TODO: remove leading and trailling Ns

		  // start new sequence
		  if (startNewSequence) {
			  newSeq = line;
			  startNewSequence = false;
		  }

		  // extend existing sequence
		  else {
			  newSeq += line;
		  }
		  *(--seq_lengths.end())+=line.length();
	  }
  }

  infile.close();

  // append the last sequence to the StringSet
  if ( newSeq != "" ) {
	  seqan::appendValue(seqs, newSeq);
  }

  return seq_lengths.size();

}

int KixBuild::save_fmindex( FMIndex & idx, const std::string &iname ) {
	const char * filename = iname.c_str();
	return seqan::save(idx, filename);
}

int KixBuild::save_metadata( const std::string &iname ) {

	bool success = save_seqnames(iname);
	success = success && save_seqlengths(iname);

	return success;

}

int KixBuild::save_seqnames( const std::string &iname ) {
	std::string file_name = iname;
	file_name += ".seqnames";
	// get total amount of bytes needed
	uint32_t data_size = sizeof(uint32_t); // number of sequences
	for ( CountType i = 0; i < seq_names.size(); i++ ) {
		data_size += sizeof(uint32_t); // length of the name
		data_size += seq_names[i].size();
	}
	// init data vector
	std::vector<char> data(data_size);
	char* d = data.data();

	// write number of sequences
	uint32_t size = seq_names.size();
	memcpy(d,&size,sizeof(uint32_t));
	d += sizeof(uint32_t);

	// write all sequence names
	for ( uint32_t i = 0; i < seq_names.size(); i++ ) {
		// write sequence name length
		size = seq_names[i].size();
		memcpy(d,&size,sizeof(uint32_t));
		d += sizeof(uint32_t);
		const char* seqn = seq_names[i].c_str();
			// write sequence name
			memcpy(d,seqn,seq_names[i].size());
			d += seq_names[i].size();
		}

		// write data to file
		uint32_t written = write_binary_file(file_name, data);

		return ( written == data_size );
}

int KixBuild::save_seqlengths( const std::string &iname ) {

		std::string file_name = iname;
		file_name += ".seqlengths";

		// get total amount of bytes needed
		uint32_t data_size = sizeof(uint32_t); // number of sequences
		for ( CountType i = 0; i < seq_lengths.size(); i++ ) {
			data_size += sizeof(uint32_t); // sequence length
		}

		// init data vector
		std::vector<char> data(data_size);
		char* d = data.data();

		// write number of sequences
		uint32_t seqlen = seq_lengths.size();
		memcpy(d,&seqlen,sizeof(uint32_t));
		d += sizeof(uint32_t);

		// write all sequence names
		for ( uint32_t i = 0; i < seq_lengths.size(); i++ ) {

			// write sequence length
			seqlen = seq_lengths[i];
			memcpy(d,&seqlen,sizeof(uint32_t));
			d += sizeof(uint32_t);

		}

		// write data to file
		uint32_t written = write_binary_file(file_name, data);

		return ( written == data_size );
}




////////////////////////////
////////// KixRun //////////
////////////////////////////


int KixRun::load_fmindex( std::string index_name ) {

	const char * filename = index_name.c_str();

	if ( seqan::open(idx, filename) != 1 ) // function returns 1 on success
		return 1;
	return 0;


}

int KixRun::load_metadata( const std::string &iname ) {

	bool success = load_seqnames(iname);
	success = success && load_seqlengths(iname);
	return success;
}

int KixRun::load_seqnames( const std::string &iname ) {

	std::string file_name = iname;
	file_name += ".seqnames";

	std::vector<char> data = read_binary_file(file_name);
	char* d = data.data();

	uint32_t num_seqs;
	memcpy(&num_seqs, d, sizeof(uint32_t));
	d += sizeof(uint32_t);

	uint32_t curr_seqLength;
	for ( uint32_t i=0; i<num_seqs; i++) {

		memcpy(&curr_seqLength, d, sizeof(uint32_t));
		d += sizeof(uint32_t);

		char* seq_name = new char[curr_seqLength];

		memcpy(seq_name, d, curr_seqLength);
		d += curr_seqLength;

		// Don't know why the substring is necessary here, but without it there are strange artifacts in some names.
		seq_names.push_back(std::string(seq_name).substr(0,curr_seqLength));

	}
	return true;
}

int KixRun::load_seqlengths( const std::string &iname ) {
	std::string file_name = iname;
	file_name += ".seqlengths";

	std::vector<char> data = read_binary_file(file_name);
	char* d = data.data();

	uint32_t num_seqs;
	memcpy(&num_seqs, d, sizeof(uint32_t));
	d += sizeof(uint32_t);

	uint32_t seqlength;
	for ( uint32_t i=0; i<num_seqs; i++) {
		memcpy(&seqlength, d, sizeof(uint32_t));
		d += sizeof(uint32_t);

		seq_lengths.push_back(seqlength);

	}
	return true;
}

