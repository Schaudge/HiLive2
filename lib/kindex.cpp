#include "kindex.h"


int KixBuild::add_fasta(const std::string &fname, bool convert_spaces, bool trim_ids) {
  GenomeIdListType added_ids;
  return add_fasta(fname, added_ids, convert_spaces, trim_ids);
}

int KixBuild::add_fasta(const std::string &fname, GenomeIdListType &ids, bool convert_spaces, bool trim_ids) {
  std::ios::sync_with_stdio(false);
  std::ifstream::sync_with_stdio(false);

  std::ifstream infile (fname.c_str());
  assert(infile.is_open());

  std::string line;
  GenomeIdType seq_id = 0;
  std::string  seq_name;
  GenomeIdListType added_ids;
  PositionType sequencePosition;
  bool startNewSequence = false;
  std::string tailingKmer;

  while(getline(infile, line)) {
    if (line.length() == 0) {continue;};
    
    if (line[line.length()-1] == '\r'){
      // check for Windows newline characters (just in case the fasta comes from Windows --> thanks Simon)
      line.erase(line.length()-1);
    }

    if (line[0] == '>') { // header line
      // initialize a new sequence
      num_seq += 1;
      seq_id = num_seq - 1;

      if (convert_spaces)
        std::replace( line.begin(), line.end(), ' ', '_');
      if (trim_ids)
        seq_name = line.substr(1,line.find(' ')-1);
      else
        seq_name = line.substr(1,line.length()-1);

      added_ids.push_back(seq_id);
      seq_names.push_back(seq_name);
      seq_lengths.push_back(0); // gets later corrected to globalAlignmentSettings.get_kmer_span() - 1
      assert(seq_names.size() == num_seq);
      startNewSequence = true;
    } 
    else { // sequence line
      if (startNewSequence) {
        if (line.length() < globalAlignmentSettings.get_kmer_span())
            continue; // ignore sequences shorter than K
        start_sequence(line, tailingKmer, sequencePosition);
        startNewSequence = false;
      }
      else
        continue_sequence(line, tailingKmer, sequencePosition);
    }
  }
  infile.close();
  ids = added_ids;
  return added_ids.size();
}


/* Start adding all k-mers in a sequence string s to the database.
    A new ID is created for this sequence.
    Return: sequence ID
*/
GenomeIdType KixBuild::start_sequence(const std::string &s, std::string& tailingKmer, PositionType& sequencePosition) {
  assert(seq_names.size() == num_seq);
  assert(s.length() >= globalAlignmentSettings.get_kmer_span());

  // add sequence kmers to index
  sequencePosition = 0; // use 1-based positions (to allow for negative positions)
  std::string::const_iterator it_s = s.begin();
  std::string::const_iterator last_invalid;
  HashIntoType fw; // forward k-mer

  seq_lengths.back()=globalAlignmentSettings.get_kmer_span()-1;
  for (; it_s < s.end()-globalAlignmentSettings.get_kmer_span()+1; ++it_s) {
    ++sequencePosition; // use 1-based positions (to allow for negative positions)
    seq_lengths.back()+=1;
    last_invalid = hash_fw(it_s, s.end(), fw);

    // add k-mer to database
    if (last_invalid < it_s)
      add_kmer(fw,num_seq-1,sequencePosition);
    else {
      unsigned jumplength = std::min(last_invalid - it_s, s.end() - globalAlignmentSettings.get_kmer_span() - it_s);
      sequencePosition += jumplength;
      seq_lengths.back() += jumplength;
      it_s = last_invalid;
    }
  }
  tailingKmer = s.substr(s.length()-globalAlignmentSettings.get_kmer_span());
  return num_seq;
}


/* Continue adding all k-mers in a sequence string s to the database.
    Return: sequence ID
*/
GenomeIdType KixBuild::continue_sequence(const std::string &s, std::string& tailingKmer, PositionType& sequencePosition) {
  assert(seq_names.size() == num_seq);

  std::string concatString = tailingKmer + s;
  // add sequence kmers to index
  std::string::const_iterator it_s = concatString.begin() + 1;
  std::string::const_iterator last_invalid;
  HashIntoType fw; // forward k-mer

  for (; it_s < concatString.end()-globalAlignmentSettings.get_kmer_span()+1; ++it_s) {
    ++sequencePosition; // use 1-based positions (to allow for negative positions)
    seq_lengths.back()+=1;
    last_invalid = hash_fw(it_s, concatString.end(), fw);

    // add k-mer to database
    if (last_invalid < it_s)
      add_kmer(fw,num_seq-1,sequencePosition);
    else {
      unsigned jumplength = std::min(last_invalid - it_s, concatString.end()- globalAlignmentSettings.get_kmer_span() - it_s);
      sequencePosition += jumplength;
      seq_lengths.back() += jumplength;
      it_s = last_invalid;
    }
  }
  tailingKmer = concatString.substr(concatString.length()-globalAlignmentSettings.get_kmer_span());
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

  for (auto it = db.begin(); it != db.end(); ++it)
    if ((*it).size() > max_count) {
      trimmed += (*it).size();
      (*it).clear();
      (*it).push_back(gp_trimmed);
    }
  return trimmed;
}


std::vector<char> KixBuild::serialize() {
  // first of all, sort the database entries by position
  for (auto it = db.begin(); it != db.end(); ++it)
    std::sort(it->begin(), it->end(), gp_compare);

  // calculate total size
  unsigned long int total_size = 0;

  // K_HiLive itself
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
  uint8_t kk = K_HiLive;
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
    uint32_t num_positions = (*it).size(); //db is an array of GenomePosType structs
    memcpy(d,&num_positions,sizeof(uint32_t));
    d += sizeof(uint32_t);

    // genome ID and position
    for(GenomePosListIt entry=(*it).begin(); entry!=(*it).end(); ++entry) {
      GenomeIdType gid = (*entry).gid;
      PositionType pos = (*entry).pos;

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
  assert(K_HiLive == kk);

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

    char * tmp = new char[nm_length+1];
    memcpy(tmp,d+bytes,nm_length);
    tmp[nm_length] = 0; // make the string null-terminated
    seq_names.push_back(tmp);
    delete tmp;
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


uint64_t KixRun::deserialize(char* d) {
  // the total number of bytes read
  uint64_t bytes = 0; 

  // read and check K
  uint8_t kk;
  memcpy(&kk,d+bytes,1);
  bytes++;
  assert(kk == K_HiLive);

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

    char * tmp = new char[nm_length+1];
    memcpy(tmp,d+bytes,nm_length);
    tmp[nm_length] = 0; // make the string null-terminated
    seq_names.push_back(tmp);
    delete tmp;
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
  deserialize(sdata.data());

  return sdata.size();
}


/* Retrieve all occurrences (fwd & rc) of kmer in the reference from the index */
GenomePosListType KixRun::retrieve_positions(std::string kmerSpan) {

  // get the reverse complement of the kmer
  HashIntoType fwHashValue;
  HashIntoType rcHashValue;
  hash(kmerSpan.c_str(), fwHashValue, rcHashValue);
  
  // obtain the list of positions for each k-mer
  char* fwd_begin = db[fwHashValue];
  uint32_t fwd_len;
  memcpy(&fwd_len,fwd_begin,sizeof(uint32_t));
  char* rev_begin = db[rcHashValue];
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
