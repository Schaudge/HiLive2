#include "kindex.h"



int KixBuild::add_fasta(const std::string &fname) {
  GenomeIdListType added_ids;
  return add_fasta(fname, added_ids);
}

int KixBuild::add_fasta(const std::string &fname, GenomeIdListType &ids) {
  //std::cout << "Sequence " << num_seq << ".  Reading file " << fname << std::endl;
  std::ios::sync_with_stdio(false);
  std::ifstream::sync_with_stdio(false);

  std::ifstream infile (fname.c_str());
  assert(infile.is_open());

  std::string line, seq;
  HashIntoType fw; // forward k-mer
  PositionType pos = 1; // use 1-based positions (to allow for negative positions)
  GenomeIdType seq_id = 0;
  std::string  seq_name;
  GenomeIdListType added_ids;
  unsigned long int last_invalid = K+1; // number of characters since the last invalid character (e.g. N)

  while(getline(infile, line)) {
    if (line.length() == 0) {continue;};

    if (line[0] == '>') { // header line
      // initialize a new sequence
      pos = 1;
      fw = 0;
      num_seq += 1;
      seq_id = num_seq - 1;
      last_invalid = K+1;
      if (line.find(" ") != std::string::npos) {
	seq_name = line.substr(0,line.find(" "));
      } else {
	seq_name = line;
      }
      added_ids.push_back(seq_id);
      seq_names.push_back(seq_name);
      seq_lengths.push_back(0);
      assert(seq_names.size() == num_seq);
    } 
    else { // sequence line
      const char * lstr = line.c_str(); 
      if (pos==1) { // beginning of a sequence
	last_invalid = hash_fw(lstr, fw); // hash the first k-mer
	seq_lengths.back()+=1;

	if (last_invalid > K) {
	  add_kmer(fw,seq_id,pos);
	}

	for (unsigned int i = K; i < line.length(); i++) {
	  pos = i - K + 1 + 1; // use 1-based positions (to allow for negative positions)
	  seq_lengths.back()+=1;

	  if ( seq_chars.find(lstr[i]) == std::string::npos ) {
	    last_invalid = 1;
	    continue;
	  } else {
	    last_invalid += 1;
	  }

	  update_kmer(fw, twobit_repr(lstr[i]));
	
	  // add k-mer to database
	  if (last_invalid > K) {
	    add_kmer(fw,seq_id,pos);
	  }
	}
      }
      else { // continue a sequence
	for (unsigned int i = 0; i < line.length(); i++) {
	  pos++;
	  seq_lengths.back()+=1;
	
	  if ( seq_chars.find(lstr[i]) == std::string::npos ) {
	    last_invalid = 1;
	    continue;
	  } else {
	    last_invalid += 1;
	  }

	  update_kmer(fw, twobit_repr(lstr[i]));

	  // add k-mer to database
	  if (last_invalid > K) {
	    add_kmer(fw,seq_id,pos);
	  }
	}
      }
    }
    
  }
  infile.close();
  ids = added_ids;
  return added_ids.size();
}



GenomeIdType KixBuild::add_sequence(const std::string &s) {
  /* Add all k-mers in a sequence string s to the database.
     A new ID is created for this sequence.
     Return: sequence ID
   */

  assert(seq_names.size() == num_seq);

  // increase the sequence counter
  num_seq += 1; 

  // add a default name for the sequence
  seq_names.push_back(std::string("Sequence_")+std::to_string(num_seq));
  seq_lengths.push_back(0);

  PositionType pos = 1; // use 1-based positions (to allow for negative positions)
  const char * sp = s.c_str(); 
  unsigned int length = s.length();

  // the current k-mers
  HashIntoType fw; // forward k-mer
  HashIntoType last_invalid = hash_fw(sp, fw); // hash the first k-mer
  seq_lengths.back()+=1;
  if (last_invalid > K) {
    add_kmer(fw,num_seq-1,pos);
  }

  for (unsigned int i = K; i < length; i++) {
    pos = i - K + 1 + 1; // use 1-based positions (to allow for negative positions)
    seq_lengths.back()+=1;
    if ( seq_chars.find(sp[i]) == std::string::npos ) {
      last_invalid = 1;
      continue;
    } else {
      last_invalid += 1;
    }

    update_kmer(fw, twobit_repr(sp[i]));

    // add k-mer to database
    if (last_invalid > K) {
      add_kmer(fw,num_seq-1,pos);
    }
  }
  
  return num_seq;
}


int KixBuild::add_kmer(HashIntoType kmer, GenomeIdType id, PositionType pos) {
  assert(kmer < db.size());
  assert(id < num_seq);
  GenomePosType gp;
  gp.gid = id;
  gp.pos = pos;
  db[kmer].push_back(gp);
  return 1;
}

/* Trim the k-mer index: removes all k-mers with more than
   max_count occurrences in the reference genomes. Trimmed
   k-mers are marked by the GenomeIdType TRIMMED (from
   definitions.h). */
uint64_t KixBuild::trim(uint64_t max_count) {

  uint64_t trimmed = 0;
  GenomePosType gp_trimmed (TRIMMED,0);

  for (auto it = db.begin(); it != db.end(); ++it) {
    if ((*it).size() > max_count) {
      trimmed += (*it).size();
      (*it).clear();
      (*it).push_back(gp_trimmed);
    }
  }

  return trimmed;
}



std::string KixBuild::get_name(GenomeIdType id) {
  return seq_names[id];
}


std::vector<char> KixBuild::serialize() {
  // first of all, sort the database entries by position
  for (auto it = db.begin(); it != db.end(); ++it) {
    std::sort(it->begin(), it->end(), gp_compare);
  }

  // calculate total size
  unsigned long int total_size = 0;

  // K itself
  total_size += 1;

  // total number of sequences in database
  total_size += sizeof(GenomeIdType);

  // sequence names
  for (uint32_t i = 0; i < seq_names.size(); i++) {
    uint16_t nm_length = seq_names.size();
    total_size += sizeof(uint16_t);
    total_size += nm_length;
  }

  // reference sequence lengths
  for (uint32_t i = 0; i < seq_lengths.size(); i++) {
    total_size += sizeof(uint64_t);
  }
  
  // database entries
  for (auto it = db.begin(); it != db.end(); ++it) {
    total_size += sizeof(uint32_t); // number of positions

    uint32_t num_positions = (*it).size();
    total_size += num_positions*(sizeof(GenomeIdType) + sizeof(PositionType));
  }
  

  // create the vector to store the data
  std::vector<char> data (total_size);
  char* d = data.data();

  // write K
  uint8_t kk = K;
  memcpy(d,&kk,1);
  d++;

  // total number of sequences in database
  memcpy(d,&num_seq,sizeof(GenomeIdType));
  d += sizeof(GenomeIdType);

  // sequence names
  for (uint32_t i = 0; i < seq_names.size(); i++) {
    uint16_t nm_length = seq_names[i].size();

    memcpy(d,&nm_length,sizeof(uint16_t));
    d += sizeof(uint16_t);
    
    memcpy(d,seq_names[i].c_str(),nm_length);
    d += nm_length;
  }
  for (uint32_t i = 0; i < seq_lengths.size(); i++) {
    uint64_t seq_len = seq_lengths[i];
    memcpy(d,&seq_len,sizeof(uint64_t));
    d += sizeof(uint64_t);
  }
  
  // database entries
  for (auto it = db.begin(); it != db.end(); ++it) {
    // number of positions
    uint32_t num_positions = (*it).size();
    memcpy(d,&num_positions,sizeof(uint32_t));
    d += sizeof(uint32_t);

    // genome ID and position
    for( uint32_t i = 0; i < num_positions; i++) {
      GenomeIdType gid = (*it)[i].gid;
      PositionType pos = (*it)[i].pos;

      memcpy(d,&gid,sizeof(GenomeIdType));
      d += sizeof(GenomeIdType);
      
      memcpy(d,&pos,sizeof(PositionType));
      d += sizeof(PositionType);
    }
  }

  return data;
}

uint64_t KixBuild::serialize_file(std::string f) {
  std::string fname = f;

  // serialize data
  std::vector<char> sdata = serialize();

  // open binary file
  FILE* ofile;
  ofile = fopen(fname.c_str(), "wb");

  if (!ofile) {
    std::cerr << "Error serializing object to file " << fname << ": Could not open file for writing." << std::endl;
    return 0;
  }

  // write all data
  uint64_t written = fwrite(sdata.data(), 1, sdata.size(), ofile);
  
  // close file
  fclose(ofile);

  if (written != sdata.size()){
    std::cerr << "Error serializing object to file " << fname << ": Total size: " << sdata.size() << " bytes. Written: " << written << " bytes." << std::endl;
  }
  
  return written;
}


uint64_t KixBuild::deserialize(char* d) {
  // the total number of bytes read
  uint64_t bytes = 0; 

  // read and check K
  uint8_t kk;
  memcpy(&kk,d+bytes,1);
  bytes++;
  assert(K == kk);

  // read total number of sequences in database
  memcpy(&num_seq,d+bytes,sizeof(GenomeIdType));
  bytes += sizeof(GenomeIdType);


  // sequence names
  seq_names.clear();
  seq_lengths.clear();
  for (uint32_t i = 0; i < num_seq; i++) {
    uint16_t nm_length;
    memcpy(&nm_length,d+bytes,sizeof(uint16_t));
    bytes += sizeof(uint16_t);

    char tmp[nm_length+1];
    memcpy(tmp,d+bytes,nm_length);
    tmp[nm_length] = 0; // make the string null-terminated
    seq_names.push_back(tmp);
    bytes += nm_length;
  }
  for (uint32_t i = 0; i < num_seq; i++) {
    uint64_t seq_len;
    memcpy(&seq_len,d+bytes,sizeof(uint64_t));
    bytes += sizeof(uint64_t);
    seq_lengths.push_back(seq_len);
  }

  
  // database entries
  for (auto it = db.begin(); it != db.end(); ++it) {
    // number of positions
    uint32_t num_positions;
    memcpy(&num_positions,d+bytes,sizeof(uint32_t));
    bytes += sizeof(uint32_t);
    (*it).clear();
    (*it).reserve(num_positions);

    // genome ID and position
    for( uint32_t i = 0; i < num_positions; i++) {
      GenomeIdType gid;
      PositionType pos;
      GenomePosType gp;
      
      memcpy(&gid,d+bytes,sizeof(GenomeIdType));
      bytes += sizeof(GenomeIdType);
      
      memcpy(&pos,d+bytes,sizeof(PositionType));
      bytes += sizeof(PositionType);
      
      gp.gid = gid;
      gp.pos = pos;
      
      (*it).push_back(gp);
    }
  }

  return bytes;
}


uint64_t KixBuild::deserialize_file(std::string f) {
  std::string fname = f;
  
  // obtain file size
  unsigned long int size = get_filesize(fname);

  // open binary file
  FILE* ifile;
  ifile = fopen(fname.c_str(), "rb");

  if (!ifile) {
    std::cerr << "Error reading from file " << fname << ": Could not open file." << std::endl;
    return 0;
  }

  // allocate memory
  std::vector<char> sdata (size);
  
  // read all data
  unsigned long int read = fread(sdata.data(), 1, size, ifile);

  // close file
  fclose (ifile);

  if (read != size){
    std::cerr << "Error reading from file " << fname << ": File size: " << size << " bytes. Read: " << read << " bytes." << std::endl;
    return 0;
  }

  // deserialize data
  deserialize(sdata.data());

  return read;
}






bool gp_compare (GenomePosType i,GenomePosType j) { 
  return (i.pos < j.pos); 
}





uint64_t KixRun::deserialize(char* d) {
  // the total number of bytes read
  uint64_t bytes = 0; 

  // read and check K
  uint8_t kk;
  memcpy(&kk,d+bytes,1);
  bytes++;
  assert(kk == K);

  // read total number of sequences in database
  memcpy(&num_seq,d+bytes,sizeof(GenomeIdType));
  bytes += sizeof(GenomeIdType);

  // sequence names
  seq_names.clear();
  seq_lengths.clear();
  for (uint32_t i = 0; i < num_seq; i++) {
    uint16_t nm_length;
    memcpy(&nm_length,d+bytes,sizeof(uint16_t));
    bytes += sizeof(uint16_t);

    char tmp[nm_length+1];
    memcpy(tmp,d+bytes,nm_length);
    tmp[nm_length] = 0; // make the string null-terminated
    seq_names.push_back(tmp);
    bytes += nm_length;
  }

  // sequence lengths
  for (uint32_t i = 0; i < num_seq; i++) {
    uint64_t seq_len;
    memcpy(&seq_len,d+bytes,sizeof(uint64_t));
    bytes += sizeof(uint64_t);
    seq_lengths.push_back(seq_len);
  }

  
  // database entries
  for (auto it = db.begin(); it != db.end(); ++it) {
    // number of positions
    uint32_t* num_positions = (uint32_t*)(d+bytes);
    
    // total size of data block for this k-mer
    uint32_t total_size = sizeof(uint32_t) + (*num_positions)*GenomePos_size;

    // allocate the memory
    (*it) = d+bytes;

    // increase pointer
    bytes += total_size;
  }

  return bytes;
}



uint64_t KixRun::deserialize_file(std::string f) {
  std::string fname = f;
  
  sdata = read_binary_file(f);
  /*
  // obtain file size
  unsigned long int size = get_filesize(fname);

  // open binary file
  FILE* ifile;
  ifile = fopen(fname.c_str(), "rb");

  if (!ifile) {
    std::cerr << "Error reading from file " << fname << ": Could not open file." << std::endl;
    return 0;
  }

  // allocate memory
  sdata.resize(size,0);
  
  // read all data
  unsigned long int read = fread(sdata.data(), 1, size, ifile);

  // close file
  fclose (ifile);

  if (read != size){
    std::cerr << "Error reading from file " << fname << ": File size: " << size << " bytes. Read: " << read << " bytes." << std::endl;
    return 0;
  }
  */
  // deserialize data
  deserialize(sdata.data());

  return sdata.size();
}



std::string KixRun::get_name(GenomeIdType id) {
  return seq_names[id];
}


/* Retrieve all occurrences (fwd & rc) of kmer in the reference from the index */
GenomePosListType KixRun::retrieve_positions(HashIntoType kmer) {

  // get the reverse complement of the kmer
  HashIntoType kmer_rc = rc(kmer);
  
  // obtain the list of positions for each k-mer
  char* fwd_begin = db[kmer];
  uint32_t fwd_len;
  memcpy(&fwd_len,fwd_begin,sizeof(uint32_t));
  char* rev_begin = db[kmer_rc];
  uint32_t rev_len;
  memcpy(&rev_len,rev_begin,sizeof(uint32_t));

  // the position list: all positions in all genomes, where the current k-mer was found
  GenomePosListType pos;
  pos.reserve(fwd_len+rev_len);
  
  // indicate reverse complement hits by negative position and append in reverse order
  for(uint64_t i = 0; i < rev_len; i++) {
    GenomePosType rev;
    memcpy(&rev.gid,rev_begin+sizeof(uint32_t)+(rev_len-1-i)*GenomePos_size, sizeof(GenomeIdType));
    memcpy(&rev.pos,rev_begin+sizeof(uint32_t)+(rev_len-1-i)*GenomePos_size+sizeof(GenomeIdType), sizeof(PositionType));
    rev.pos = - rev.pos;
    if (rev.gid == TRIMMED) {
      pos.clear();
      pos.push_back(GenomePosType(TRIMMED,0));
      return pos;
    }
    pos.push_back(rev);
  }

  // then add the forward hits with positive position in normal order --> position list is sorted if index was sorted!
  for(uint64_t i = 0; i < fwd_len; i++) {
    GenomePosType fwd;
    memcpy(&fwd.gid,fwd_begin+sizeof(uint32_t)+i*GenomePos_size, sizeof(GenomeIdType));
    memcpy(&fwd.pos,fwd_begin+sizeof(uint32_t)+i*GenomePos_size+sizeof(GenomeIdType), sizeof(PositionType));
    if (fwd.gid == TRIMMED) {
      pos.clear();
      pos.push_back(GenomePosType(TRIMMED,0));
      return pos;
    }
    pos.push_back(fwd);
  }

  return pos;
}


// generate a SAM compliant header string
std::string KixRun::get_SAM_header() {
  std::stringstream s;
  
  // header
  s << "@HD\tVN:1.5\tSO:unsorted" << std::endl;

  // sequence information
  assert(seq_names.size() == seq_lengths.size());
  for ( uint64_t i = 0; i < seq_names.size(); i++  ) {
    s << "@SQ\tSN:" << seq_names[i] << "\tLN:" << seq_lengths[i] << std::endl;
  }

  // program information
  s << "@PG\tHiLive" << std::endl;


  return s.str();
}

