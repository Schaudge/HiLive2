#include "alnstream.h"


//-------------------------------------------------------------------//
//------  The output Alignment Stream class  ------------------------//
//-------------------------------------------------------------------//

// new output Alignment Stream class
oAlnStream::oAlnStream(uint16_t ln, uint16_t tl, uint16_t cl, std::string rt, CountType rl, uint32_t nr, uint64_t bs, uint8_t fmt):
  lane(ln), tile(tl), cycle(cl), root(rt), rlen(rl), num_reads(nr), num_written(0), buffer(bs,0), buf_size(bs), buf_pos(0), format(fmt), ofile(NULL), ozfile(Z_NULL) {}


// write function for lz4 compression
uint64_t oAlnStream::lz4write(const char* source, uint64_t size) {
  // allocate buffer for the compressed data
  std::vector<char> buf (LZ4_COMPRESSBOUND(size),0);
  
  // compress the data
  uint32_t compressed_size = LZ4_compress (source, buf.data(), size);
  if (!compressed_size)
    throw std::runtime_error("Error compressing data with LZ4.");
  
  // write the block size
  if ( !fwrite(&compressed_size, 1, sizeof(uint32_t), ofile) )
    throw std::runtime_error("Error writing block size to file while compressing data with LZ4.");

  // write the data chunk
  if ( !fwrite(buf.data(), 1, compressed_size, ofile) )
    throw std::runtime_error("Error writing data to file while compressing with LZ4.");
  
  return size;
}


uint64_t oAlnStream::open(std::string fname) {
  // open the new Alignment file
  switch (format) {
  case 0: case 2: 
    ofile = fopen(fname.c_str(), "wb");
    if (!ofile) {
      std::cerr << "Could not open file " << fname << " for writing." << std::endl;
      return 0;
    }
    break;
  case 1:
    ozfile = gzopen(fname.c_str(), "wb1"); //Don't compress too much, not enough bang for the buck
    if (ozfile == Z_NULL) {
      std::cerr << "Could not open file " << fname << " for writing." << std::endl;
      return 0;
    }
    break;
  default:
    throw std::invalid_argument("Output file format not recognized");
  }

  // write the header:

  // calculate total size first
  unsigned long int total_size = 0;

  total_size += sizeof(uint16_t); // lane
  total_size += sizeof(uint16_t); // tile
  total_size += sizeof(CountType); // cycle

  // root directory name
  uint16_t root_size = root.size();
  total_size += sizeof(uint16_t);
  total_size += root_size;

  // read length
  total_size += sizeof(CountType);

  // number of reads
  total_size += sizeof(uint32_t);

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
  memcpy(d,&num_reads,sizeof(uint32_t));
  d += sizeof(int32_t);

  // write all data
  uint64_t written = 0;
  switch (format) {
  case 0: case 2:  written = fwrite(data.data(), 1, data.size(), ofile); break;
  case 1: written = gzwrite(ozfile, data.data(), data.size()); break;
  }
  
  return written;
}


// writes a read alignment to the output Alignment file. 
// Buffering is handled internally
uint64_t oAlnStream::write_alignment(ReadAlignment * al) {
  if ( (!ofile && (format == 0 || format == 2)) || (ozfile == Z_NULL && format == 1) ){
    throw std::runtime_error("Could not write alignment to file. File handle not valid.");
  }
  if (num_written >= num_reads) {
    throw std::length_error("Could not write alignment to file. All alignments were already written.");
  }

  std::vector<char> data = al->serialize();
  uint64_t al_size = data.size();

  // first, write the size of the serialized alignment (uint32_t = 4 bytes)
  if (buf_pos+sizeof(uint32_t) <= buf_size) {
    // directly copy if all 4 bytes have space in the buffer (should be almost always the case)
    memcpy(buffer.data()+buf_pos,&al_size,sizeof(uint32_t));
    buf_pos += sizeof(uint32_t);
  } 
  else {
    // copy the first bytes into temporary buffer to compose the alignment size
    std::vector<char> temp (sizeof(uint32_t),0);
    memcpy(temp.data(),&al_size,sizeof(uint32_t));

    uint64_t first_part = buf_size-buf_pos;
    memcpy(buffer.data()+buf_pos,temp.data(),first_part);

    // write out buffer
    uint64_t written = 0;
    switch (format) {
    case 0: written = fwrite(buffer.data(), 1, buffer.size(), ofile); break;
    case 1: written = gzwrite(ozfile, buffer.data(), buffer.size()); break;
    case 2: written = lz4write(buffer.data(), buffer.size()); break;
    }
    assert(written == buf_size);

    // copy remaining data
    memcpy(buffer.data(),temp.data()+first_part,sizeof(uint32_t)-first_part);
    buf_pos = sizeof(uint32_t)-first_part;
  }

  // finally, write the serialized data
  uint64_t copied = 0;
  while (copied < al_size) {
    uint64_t to_copy = std::min(al_size-copied,buf_size-buf_pos);
    memcpy(buffer.data()+buf_pos, data.data()+copied, to_copy);
    buf_pos += to_copy;
    copied += to_copy;

    // write buffer to disk if full
    if(buf_pos >= buf_size){
      uint64_t written = 0;
      switch (format) {
      case 0: written = fwrite(buffer.data(), 1, buffer.size(), ofile); break;
      case 1: written = gzwrite(ozfile, buffer.data(), buffer.size()); break;
      case 2: written = lz4write(buffer.data(), buffer.size()); break;
      }
      assert(written == buf_size);
      buf_pos = 0;
    }
  }

  num_written++;

  return num_written;  
}


bool oAlnStream::close() {
  if ( ((format == 0 || format == 2) && ofile) || (format == 1 && ozfile != Z_NULL) ) {
    // write remaining buffer content to file
    uint64_t written = 0;
    switch (format) {
    case 0: written = fwrite(buffer.data(), 1, buf_pos, ofile); break;
    case 1: written = gzwrite(ozfile, buffer.data(), buf_pos); break;
    case 2: written = lz4write(buffer.data(), buf_pos); break;
    }
    assert(written == buf_pos);
    buf_pos = 0;
    if (num_written == num_reads) {
      switch (format) {
      case 0: case 2: fclose(ofile); break;
      case 1: gzclose(ozfile); break;
      }
      return true;
    }
    else {
      std::cerr << "Error: Could not close output alignment file! "<< num_reads - num_written <<" alignments missing." << std::endl;
      return false;
    }
  }
  else {
    std::cerr << "Error: Could not close output alignment file! File handle not valid." << std::endl;
    return false;
  }
}



//-------------------------------------------------------------------//
//------  The input Alignment Stream class  -------------------------//
//-------------------------------------------------------------------//

// new Alignment Stream class
iAlnStream::iAlnStream(uint64_t bs, uint8_t fmt):
  lane(0), tile(0), cycle(0), root(""), rlen(0), num_reads(0), num_loaded(0), buffer(bs,0), buf_size(bs), buf_pos(bs), format(fmt), ifile(NULL), izfile(Z_NULL) {}


// read function for lz4 decompression, reads one block of data
uint64_t iAlnStream::lz4read_block() {
  // get the size of the next block
  uint32_t compressed_size = 0;
  if ( !fread(&compressed_size,sizeof(uint32_t),1,ifile) )
    return 0;
  
  // allocate buffer for the compressed data
  std::vector<char> cbuf (compressed_size,0);

  // read the data
  if ( !fread(cbuf.data(),compressed_size,1,ifile) )
    throw std::runtime_error("Malformed input file. Could not read next block.");
  
  // decompress the data
  int64_t r_size = LZ4_decompress_safe (cbuf.data(), buffer.data(), compressed_size, buffer.size());
  if ( r_size < 0 )
    throw std::runtime_error("Error while decompressing LZ4 compressed block.");
  
  // update the current buffer size
  buf_size = r_size;

  return (uint64_t)r_size;
}


uint64_t iAlnStream::open(std::string fname) {
  // open the new Alignment file
  switch (format) {
  case 0: case 2:
    ifile = fopen(fname.c_str(), "rb");
    if (!ifile) {
      std::cerr << "Error opening file " << fname << " for reading." << std::endl;
      return 0;
    }
    break;
  case 1:
    izfile = gzopen(fname.c_str(), "rb");
    if (izfile == Z_NULL) {
      std::cerr << "Error opening file " << fname << " for reading." << std::endl;
      return 0;
    }
    break;
  default:
    throw std::invalid_argument("Input file format not recognized.");
  }

  // load the header:

  uint64_t bytes = 0;
  uint16_t root_size;
  switch (format) {
  case 0: case 2:
    {
      // read the lane
      bytes += fread(&lane,sizeof(uint16_t),1,ifile);
      // read the tile
      bytes += fread(&tile,sizeof(uint16_t),1,ifile);
      // read the cycle
      bytes += fread(&cycle,sizeof(CountType),1,ifile);
      // root directory name size
      bytes += fread(&root_size,sizeof(uint16_t),1,ifile);
      // root name
      char * tmp = new char[root_size+1];
      bytes += fread(tmp,1,root_size,ifile);
      tmp[root_size] = 0; // make the string null-terminated
      root = tmp;
      delete tmp;
      // read the read length
      bytes += fread(&rlen,sizeof(CountType),1,ifile);
      // read the number of reads
      bytes += fread(&num_reads,sizeof(uint32_t),1,ifile);
      break;
    }
  case 1:
    {
      // read the lane
      bytes += gzread(izfile,&lane,sizeof(uint16_t));
      // read the tile
      bytes += gzread(izfile,&tile,sizeof(uint16_t));
      // read the cycle
      bytes += gzread(izfile,&cycle,sizeof(CountType));
      // root directory name size
      bytes += gzread(izfile,&root_size,sizeof(uint16_t));
      // root name
      char * tmp = new char[root_size+1];
      bytes += gzread(izfile,tmp,root_size);
      tmp[root_size] = 0; // make the string null-terminated
      root = tmp;
      delete tmp;
      // read the read length
      bytes += gzread(izfile,&rlen,sizeof(CountType));
      // read the number of reads
      bytes += gzread(izfile,&num_reads,sizeof(uint32_t));
      break;
    }
  }

  return bytes;
}



ReadAlignment* iAlnStream::get_alignment() {

  if ( (format==0 && !ifile) || (format==1 && izfile == Z_NULL) ){
    throw std::runtime_error("Could not load alignment from file. File handle not valid.");
  }
  if (num_loaded >= num_reads) {
    throw std::length_error("Could not load alignment from file. All alignments were already loaded.");
  }

  // first, get the size of the serialized alignment (uint32_t = 4 bytes)
  uint32_t al_size = 0;
  if (buf_pos+sizeof(uint32_t) <= buf_size) {
    // directly copy if all 4 bytes are in the buffer (should be almost always the case)
    memcpy(&al_size,buffer.data()+buf_pos,sizeof(uint32_t));
    buf_pos += sizeof(uint32_t);
  } 
  else {
    // copy the first bytes into temporary buffer to compose the alignment size
    std::vector<char> temp (sizeof(uint32_t),0);
    uint64_t first_part = buf_size-buf_pos;
    memcpy(temp.data(),buffer.data()+buf_pos,first_part);

    // load new buffer
    uint64_t loaded;
    switch (format) {
    case 0:
      loaded = fread(buffer.data(),1,buf_size,ifile);
      assert( (loaded == buf_size) || feof(ifile) );
      break;
    case 1:
      loaded = gzread(izfile,buffer.data(),buf_size);    
      assert( (loaded == buf_size) || gzeof(izfile) );
      break;
    case 2:
      loaded = lz4read_block();    
      assert( loaded>0 || feof(ifile) );
      break;
    }

    // copy remaining data and copy to variable
    memcpy(temp.data()+first_part,buffer.data(),sizeof(uint32_t)-first_part);
    buf_pos = sizeof(uint32_t)-first_part;
    memcpy(&al_size,temp.data(),sizeof(uint32_t));
  }
  
  // then, copy the content to the data vector
  std::vector<char> data(al_size,0);
  uint64_t copied = 0;
  while (copied < al_size) {
    uint64_t to_copy = std::min(al_size-copied,buf_size-buf_pos);
    memcpy(data.data()+copied, buffer.data()+buf_pos, to_copy);
    buf_pos += to_copy;
    copied += to_copy;

    // read new buffer from disk if necessary
    if(buf_pos >= buf_size){
      uint64_t loaded;
      switch (format) {
      case 0:
	loaded = fread(buffer.data(),1,buf_size,ifile);
	assert( (loaded == buf_size) || feof(ifile) );
	break;
      case 1:
	loaded = gzread(izfile,buffer.data(),buf_size);    
	assert( (loaded == buf_size) || gzeof(izfile) );
	break;
      case 2:
	loaded = lz4read_block();    
	assert( loaded>0 || feof(ifile) );
	break;
      }
      buf_pos = 0;
    }
  }

  // finally, deserialize the alignment
  ReadAlignment* ra = new ReadAlignment();
  ra->set_rlen(rlen);
  ra->deserialize(data.data());
  
  num_loaded++;

  return ra;
}


bool iAlnStream::close() {
  //if (ifile) {
  if ( ((format==0 || format==2) && ifile) || (format==1 && izfile != Z_NULL)) {
    if (num_loaded == num_reads) {
      switch (format) {
      case 0: case 2: fclose(ifile); break;
      case 1: gzclose(izfile); break;
      }
      return true;
    }
    else {
      std::cerr << "Error: Could not close alignment file! "<< num_reads - num_loaded <<" alignments missing." << std::endl;
      return false;
    }
  }
  else {
    throw std::runtime_error("Could not close alignment file. File handle not valid.");
  }
}




//-------------------------------------------------------------------//
//------  The StreamedAlignment class  ------------------------------//
//-------------------------------------------------------------------//

StreamedAlignment& StreamedAlignment::operator=(const StreamedAlignment& other) {
  if(&other == this)
    return *this;
  
  lane = other.lane;
  tile = other.tile;
  root = other.root;
  rlen = other.rlen;
  
  return *this;
}

std::string StreamedAlignment::get_bcl_file(uint16_t cycle) {
  std::ostringstream path_stream;

  path_stream << root << "/L00" << lane << "/C" << cycle << ".1/s_"<< lane <<"_" << tile << ".bcl";
  return path_stream.str();
}


// get the path to the alignment file. The alignment file is located in
// <base>/L00<lane>/s_<lane>_<tile>.<cycle>.align
// if base == "": base = root
std::string StreamedAlignment::get_alignment_file(uint16_t cycle, std::string base){
  if (base == "") {
    base = root;
  }
  std::ostringstream path_stream;
  path_stream << base << "/L00" << lane << "/s_"<< lane << "_" << tile << "." << cycle << ".align";
  return path_stream.str();
}


std::string StreamedAlignment::get_filter_file() {
  std::ostringstream path_stream;
  path_stream << root << "/L00" << lane << "/s_"<< lane << "_" << tile << ".filter";
  return path_stream.str();
}


// create directories required to store the alignment files (only if not stored in root)
void StreamedAlignment::create_directories(AlignmentSettings* settings) {
  std::ostringstream path_stream;
  if (settings->temp_dir == "") {
    path_stream << root;
  }
  else {
    path_stream << settings->temp_dir;
  }
  path_stream << "/L00" << lane;

  boost::filesystem::create_directories(path_stream.str());

  std::ostringstream sam_stream;
  sam_stream << settings->out_dir;
  sam_stream << "/L00" << lane;
  
  boost::filesystem::create_directories(sam_stream.str());
}


// initialize empty alignment. Creates alignment files for a virtual Cycle 0
void StreamedAlignment::init_alignment(AlignmentSettings* settings) {
  std::string out_fname = get_alignment_file(0, settings->temp_dir);

  // get the number of reads in this tile by looking in the first bcl file
  std::string first_cycle = get_bcl_file(1);

  // extract the number of reads
  uint32_t num_reads = num_reads_from_bcl(first_cycle);
  
  // open output alignment stream
  oAlnStream output (lane, tile, 0, root, rlen, num_reads, settings->block_size, settings->compression_format);
  output.open(out_fname);

  // write empty read alignments for each read
  for (uint32_t i = 0; i < num_reads; ++i) {
    ReadAlignment * ra = new ReadAlignment();
    ra->set_rlen(rlen);
    output.write_alignment(ra);
    delete ra;
  }

  if(!output.close()) {
    std::cerr << "Error: Could not create initial alignment file." << std::endl;
  }
} 



// extend an existing alignment from cycle <cycle-1> to <cycle>. returns the number of seeds
uint64_t StreamedAlignment::extend_alignment(uint16_t cycle, KixRun* index, AlignmentSettings* settings) {

  // 1. Open the input file
  //-----------------------
  std::string in_fname = get_alignment_file(cycle-1, settings->temp_dir);
  std::string bcl_fname = get_bcl_file(cycle);
  std::string filter_fname = get_filter_file();

  iAlnStream input ( settings->block_size, settings->compression_format );
  input.open(in_fname);
  assert(input.get_cycle() == cycle-1);
  assert(input.get_lane() == lane);
  assert(input.get_tile() == tile);
  assert(input.get_root() == root);
  assert(input.get_rlen() == rlen);

  uint32_t num_reads = input.get_num_reads();


  // 2. Open output stream
  //----------------------------------------------------------
  std::string out_fname = get_alignment_file(cycle, settings->temp_dir);
  oAlnStream output (lane, tile, cycle, root, rlen, num_reads, settings->block_size, settings->compression_format);
  output.open(out_fname);


  
  // 3. Read the full BCL file (this is not too much)
  //-------------------------------------------------
  BclParser basecalls;
  basecalls.open(bcl_fname);
  
  // extract the number of reads from the BCL file
  uint32_t num_base_calls = basecalls.size();
  assert(num_base_calls == num_reads);


  // 4. Load the filter flags if filter file is available
  // ----------------------------------------------------
  FilterParser filters;
  if (file_exists(filter_fname)) {
    filters.open(filter_fname);
    // extract the number of reads from the filter file
    uint32_t num_reads_filter = filters.size();
    
    if (num_reads != num_reads_filter){
      std::string msg = std::string("Number of reads in filter file (") + std::to_string(num_reads_filter) + ") does not match the number of reads in the BCL file (" + std::to_string(num_reads) + ").";
      throw std::length_error(msg.c_str());
    }
  }

  // 5. Extend alignments 1 by 1
  //-------------------------------------------------
  uint64_t num_seeds = 0;
  for (uint64_t i = 0; i < num_reads; ++i) {

    ReadAlignment* ra = input.get_alignment();
    if (filters.size() > 0 && filters.has_next()) {
      // filter file was found -> apply filter
      if(filters.next()) {
        ra->extend_alignment(basecalls.next(), index, settings);
        num_seeds += ra->seeds.size();
      }
      else {
        basecalls.next();
        ra->disable(*settings);
      }
    }
    // filter file was not found -> treat every alignment as valid
    else {
      ra->extend_alignment(basecalls.next(), index, settings);
      num_seeds += ra->seeds.size();
    }

    output.write_alignment(ra);
    delete ra;
  }


  // 6. Close files
  //-------------------------------------------------
  if (!(input.close() && output.close())) {
    std::cerr << "Could not finish alignment!" << std::endl;
  }
  
  // 7. Delete old alignment file, if requested
  //-------------------------------------------
  if ( ! settings->keep_aln_files ) {
    std::remove(in_fname.c_str());
  }

  return num_seeds;
}

//-------------------------------------------------------------------//
//------  Streamed SAM generation -----------------------------------//
//-------------------------------------------------------------------//

uint64_t alignments_to_sam(uint16_t ln, uint16_t tl, std::string rt, CountType rl, KixRun* index, AlignmentSettings* settings) {
  // set the file names
  std::string temp;
  if (settings->temp_dir == "")
    temp = rt;
  else
    temp = settings->temp_dir;

  std::string sam_dir = settings->out_dir;
  std::string filter_fname = filter_name(rt, ln, tl);
  std::string alignment_fname = alignment_name(temp, ln, tl, rl);
  std::string sam_fname = sam_tile_name(sam_dir, ln, tl, settings->write_bam);

  // check if files exist
  if ( !file_exists(alignment_fname) )
    throw std::runtime_error(std::string("Could not create SAM file. Alignment file not found: ")+ alignment_fname);
  if ( !file_exists(filter_fname) )
    std::cerr << "Could not find .filter file: " <<  filter_fname << std::endl;

  // open the alignment file
  iAlnStream input ( settings->block_size, settings->compression_format );
  input.open(alignment_fname);
  uint64_t num_reads = input.get_num_reads();

  // open the filter file, if applicable
  FilterParser filters;
  if (file_exists(filter_fname)) {
    filters.open(filter_fname);
    // extract the number of reads from the filter file
    uint32_t num_reads_filter = filters.size();
    
    if (num_reads != num_reads_filter){
      std::string msg = std::string("Number of reads in filter file (") + std::to_string(num_reads_filter) + ") does not match the number of reads in the BCL file (" + std::to_string(num_reads) + ").";
      throw std::length_error(msg.c_str());
    }
  }


  /////////////////////////////////////////////
  // generate sam file

  // set BamFileOut object
  seqan::CharString samFileName = sam_fname;
  seqan::BamFileOut samFileOut(seqan::toCString(samFileName));
  seqan::StringSet<seqan::CharString> referenceNames;
  for (unsigned i=0; i<seqan::length(index->seq_names); i++) {
      seqan::appendValue(referenceNames, index->seq_names[i]);
  }
  seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > referenceNamesCache(referenceNames);
  seqan::BamIOContext<seqan::StringSet<seqan::CharString> > bamIOContext(referenceNames, referenceNamesCache);
  samFileOut.context = bamIOContext;


  /////////////////
  // set SAM header
  seqan::BamHeaderRecord headerRecord;
  
  // @HD header.
  seqan::clear(headerRecord);
  headerRecord.type = seqan::BAM_HEADER_FIRST;
  seqan::resize(headerRecord.tags, 2);
  headerRecord.tags[0].i1 = "VN";
  headerRecord.tags[0].i2 = "1.5";
  headerRecord.tags[1].i1 = "GO";
  headerRecord.tags[1].i2 = "query";
  seqan::writeHeader(samFileOut, headerRecord);
  
  // @SQ header.
  std::stringstream ss;
  for ( uint64_t i = 0; i < index->seq_names.size(); i++  ) {
    ss.str(std::string()); // clear string stream
    seqan::clear(headerRecord);
    headerRecord.type = seqan::BAM_HEADER_REFERENCE;
    seqan::resize(headerRecord.tags, 2);
  	headerRecord.tags[0].i1 = "SN";
  	headerRecord.tags[0].i2 = index->seq_names[i];

  	headerRecord.tags[1].i1 = "LN";
  	ss << index->seq_lengths[i];
  	headerRecord.tags[1].i2 = ss.str();

    seqan::writeHeader(samFileOut, headerRecord);
  }
  
  // @PG header.
  seqan::clear(headerRecord);
  headerRecord.type = seqan::BAM_HEADER_PROGRAM;
  seqan::resize(headerRecord.tags, 3);
  headerRecord.tags[0].i1 = "ID";
  headerRecord.tags[0].i2 = "hilive";
  headerRecord.tags[1].i1 = "PN";
  headerRecord.tags[1].i2 = "HiLive";
  headerRecord.tags[2].i1 = "VN";
  ss.str(std::string());
  ss << HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR;
  headerRecord.tags[2].i2 = ss.str();
  seqan::writeHeader(samFileOut, headerRecord);
  

  /////////////////
  // prepare sam entries
  uint64_t num_alignments = 0;

  for (uint64_t i = 0; i < num_reads; ++i) {
    ReadAlignment * ra = input.get_alignment();
	
    // if either: filter file is open and all filter flags were loaded and the filter flag is > 0
    //    or: the filter file is not available
    // then proceed else skip
    if (!((filters.size()>0 && filters.next()) || filters.size() == 0))
        continue;


    // Read name format <instrument‐name>:<run ID>:<flowcell ID>:<lane‐number>:<tile‐number>:<x‐pos>:<y‐pos>
    //readname << "<instrument>:<run-ID>:<flowcell-ID>:" << ln << ":" << tl << ":<xpos>:<ypos>:" << i;
    std::stringstream readname;
    readname << "lane." << ln << "|tile." << tl << "|read." << i;

    /////////////////
    // set sam entries
    seqan::BamAlignmentRecord record;
    bool printedFirstSeed = false;
    for (SeedVecIt it = ra->seeds.begin(); it != ra->seeds.end(); ) {

    	if ( (*ra->seeds.begin())->gid == TRIMMED ) {
    		if (  (ra->seeds.size() == 1) ) {
    			settings->trimmedReads.push_back(i);
    		}
    		++it;
    		continue;
    	}

        seqan::clear(record);

        record.qName = readname.str();

        record.rID = (*it)->gid;

        record.beginPos = ra->get_SAM_start_pos(*it, *settings)-1; // seqan expects 0-based positions, but in HiLive we use 1-based
        if (record.beginPos < 0) {
            it = ra->seeds.erase(it);
            continue;
        }

        unsigned nm_i = 0;
        record.cigar = (*it)->returnSeqanCigarString(&nm_i);

        // flag and seq
        record.flag = 0;
        record.seq = ra->getSequenceString(*settings);
        if ((*it)->start_pos < 0) { // if read matched reverse complementary
            seqan::reverseComplement(record.seq);
            record.flag |= 16;
        }
        if (printedFirstSeed) { // if current seed is secondary alignment
            record.flag |= 256;
            seqan::clear(record.seq);
            record.qual = "*";
        }


        // check if cigar string sums up to read length
        unsigned cigarElemSum = 0;
        unsigned deletionSum = 0;
        unsigned asi_score = (*it)->num_matches;
        for (seqan::Iterator<seqan::String<seqan::CigarElement<> > >::Type elem = seqan::begin(record.cigar); elem != end(record.cigar); ++elem) {
            if ((elem->operation == 'M') || (elem->operation == 'I') || (elem->operation == 'S') || (elem->operation == '=') || (elem->operation == 'X')) 
                cigarElemSum += elem->count;
            if ( (elem->operation == 'I') )
            	asi_score += elem->count - 1;
            if (elem->operation == 'D') {
                deletionSum += elem->count;
                asi_score -= 1;
            }
        }
        if (cigarElemSum != settings->seqlen) {
            std::cerr << "WARNING: Excluded an alignment of read " << record.qName << " at position " << ra->get_SAM_start_pos(*it, *settings) << " because its cigar vector had length " << cigarElemSum << std::endl;
            it = ra->seeds.erase(it);
            continue;
        }
        if (deletionSum >= settings->seqlen) {
            std::cerr << "WARNING: Excluded an alignment of read " << record.qName << " at position " << ra->get_SAM_start_pos(*it, *settings) << " because its cigar vector had " << deletionSum << " deletions" << std::endl;
            it = ra->seeds.erase(it);
            continue;
        }

        // tags
        seqan::BamTagsDict dict;
        seqan::appendTagValue(dict, "AS", ( asi_score ) );
        if (settings->seqlen < settings->rlen) // if demultiplexing is on
            seqan::appendTagValue(dict, "BC", ra->getBarcodeString(*settings));
        seqan::appendTagValue(dict, "NM", nm_i);
        record.tags = seqan::host(dict);


        // write record to disk
        seqan::writeRecord(samFileOut, record);
        ++num_alignments;
        printedFirstSeed = true;
        ++it;
    }
  }
  
  std::ofstream statsfile;
  statsfile.open(sam_fname+".stats");
  statsfile << "Number of reads\t" << num_reads << std::endl;
  statsfile << "Number of alignments\t" << num_alignments << std::endl;
  statsfile.close();

  return 0;
}
