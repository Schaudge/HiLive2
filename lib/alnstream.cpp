#include "alnstream.h"


//-------------------------------------------------------------------//
//------  The output Alignment Stream class  ------------------------//
//-------------------------------------------------------------------//


oAlnStream::oAlnStream(uint16_t ln, uint16_t tl, uint16_t cl, CountType rl, uint32_t nr, uint64_t bs, uint8_t fmt):
  lane(ln), tile(tl), cycle(cl), rlen(rl), num_reads(nr), num_written(0), buffer(bs,0), buf_size(bs), buf_pos(0), format(fmt), fstream(NULL), zfstream(Z_NULL), fname(""), flocked(false) {}


oAlnStream::~oAlnStream() {
	funlock();
}


uint64_t oAlnStream::lz4write(const char* source, uint64_t size) {
  // allocate buffer for the compressed data
  std::vector<char> buf (LZ4_COMPRESSBOUND(size),0);
  
  // compress the data
  uint32_t compressed_size = LZ4_compress (source, buf.data(), size);
  if (!compressed_size)
    throw std::runtime_error("Error compressing data with LZ4.");
  
  // write the block size
  if ( !fwrite(&compressed_size, 1, sizeof(uint32_t), fstream) )
    throw std::runtime_error("Error writing block size to file while compressing data with LZ4.");

  // write the data chunk
  if ( !fwrite(buf.data(), 1, compressed_size, fstream) )
    throw std::runtime_error("Error writing data to file while compressing with LZ4.");
  
  return size;
}


uint64_t oAlnStream::open(std::string f_name) {

	fname = f_name;
	flock();

	// open the new Alignment file
	switch (format) {
	case 0: case 2:
		fstream = fopen(fname.c_str(), "wb");
		if (!fstream) {
			funlock();
			throw file_open_error( "Error opening file " + fname + " for writing.");
			return 0;
		}
		break;
	case 1:
		zfstream = gzopen(fname.c_str(), "wb1"); //Don't compress too much, not enough bang for the buck
		if (zfstream == Z_NULL) {
			funlock();
			throw file_open_error( "Error opening file " + fname + " for writing.");
			return 0;
		}
		break;
	default:
		funlock();
		throw file_format_error("Output file format not recognized.");
	}

	// write the header:

	// calculate total size first
	unsigned long int total_size = 0;

	total_size += sizeof(uint16_t); // lane
	total_size += sizeof(uint16_t); // tile
	total_size += sizeof(CountType); // cycle

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

	// write the read length
	memcpy(d,&rlen,sizeof(CountType));
	d += sizeof(CountType);

	// write the number of reads
	memcpy(d,&num_reads,sizeof(uint32_t));
	d += sizeof(int32_t);

	// write all data
	uint64_t written = 0;
	switch (format) {
	case 0: case 2:  written = fwrite(data.data(), 1, data.size(), fstream); break;
	case 1: written = gzwrite(zfstream, data.data(), data.size()); break;
	}
  
	return written;
}


uint64_t oAlnStream::write_alignment(ReadAlignment * al) {

  if ( (!fstream && (format == 0 || format == 2)) || (zfstream == Z_NULL && format == 1) ){
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
    switch (format) {

    case 0: fwrite(buffer.data(), 1, buffer.size(), fstream); break;
    case 1: gzwrite(zfstream, buffer.data(), buffer.size()); break;
    case 2: lz4write(buffer.data(), buffer.size()); break;

    }

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
      switch (format) {

      case 0: fwrite(buffer.data(), 1, buffer.size(), fstream); break;
      case 1: gzwrite(zfstream, buffer.data(), buffer.size()); break;
      case 2: lz4write(buffer.data(), buffer.size()); break;

      }
      buf_pos = 0;
    }
  }

  num_written++;

  return num_written;  
}


bool oAlnStream::close() {
  if ( ((format == 0 || format == 2) && fstream) || (format == 1 && zfstream != Z_NULL) ) {
    // write remaining buffer content to file
    switch (format) {

    case 0: fwrite(buffer.data(), 1, buf_pos, fstream); break;
    case 1: gzwrite(zfstream, buffer.data(), buf_pos); break;
    case 2: lz4write(buffer.data(), buf_pos); break;

    }
    buf_pos = 0;
    if (num_written == num_reads) {
      switch (format) {
      case 0: case 2: fclose(fstream); break;
      case 1: gzclose(zfstream); break;
      }
      funlock();
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

void oAlnStream::flock() {
	fileLocks.lock(fname);
	flocked = true;
}


void oAlnStream::funlock() {
	if ( flocked ) {
		fileLocks.unlock(fname);
	}
}

//-------------------------------------------------------------------//
//------  The input Alignment Stream class  -------------------------//
//-------------------------------------------------------------------//

iAlnStream::iAlnStream(uint64_t bs, uint8_t fmt):lane(0), tile(0), cycle(0), rlen(0), num_reads(0), num_loaded(0), buffer(bs,0), buf_size(bs), buf_pos(bs), format(fmt), fstream(NULL), zfstream(Z_NULL), fname(""), flocked(false) {}

iAlnStream::~iAlnStream() {
	funlock();
}


// read function for lz4 decompression, reads one block of data
uint64_t iAlnStream::lz4read_block() {
  // get the size of the next block
  uint32_t compressed_size = 0;
  if ( !fread(&compressed_size,sizeof(uint32_t),1,fstream) )
    return 0;
  
  // allocate buffer for the compressed data
  std::vector<char> cbuf (compressed_size,0);

  // read the data
  if ( !fread(cbuf.data(),compressed_size,1,fstream) )
    throw std::runtime_error("Malformed input file. Could not read next block.");
  
  // decompress the data
  int64_t r_size = LZ4_decompress_safe (cbuf.data(), buffer.data(), compressed_size, buffer.size());
  if ( r_size < 0 )
    throw std::runtime_error("Error while decompressing LZ4 compressed block.");
  
  // update the current buffer size
  buf_size = r_size;

  return (uint64_t)r_size;
}


uint64_t iAlnStream::open(std::string f_name) {

	if ( !file_exists(f_name) ) {
		throw file_not_exist_error( " File " + fname + " does not exist.");
	}

	fname = f_name;
	flock();

	// open the new Alignment file
	switch (format) {
	case 0: case 2:
		fstream = fopen(fname.c_str(), "rb");
		if (!fstream) {
			funlock();
			throw file_open_error( "Error opening file " + fname + " for reading.");
			return 0;
		}
		break;
	case 1:
		zfstream = gzopen(fname.c_str(), "rb");
		if (zfstream == Z_NULL) {
			funlock();
			throw file_open_error( "Error opening file " + fname + " for reading.");
			return 0;
		}
		break;
	default:
		funlock();
		throw file_format_error("Input file format not recognized.");
	}

	// load the header:

	uint64_t bytes = 0;
	switch (format) {
	case 0: case 2:
	{
		// read the lane
		bytes += fread(&lane,sizeof(uint16_t),1,fstream);
		// read the tile
		bytes += fread(&tile,sizeof(uint16_t),1,fstream);
		// read the cycle
		bytes += fread(&cycle,sizeof(CountType),1,fstream);
		// read the read length
		bytes += fread(&rlen,sizeof(CountType),1,fstream);
		// read the number of reads
		bytes += fread(&num_reads,sizeof(uint32_t),1,fstream);
		break;
	}
	case 1:
	{
		// read the lane
		bytes += gzread(zfstream,&lane,sizeof(uint16_t));
		// read the tile
		bytes += gzread(zfstream,&tile,sizeof(uint16_t));
		// read the cycle
		bytes += gzread(zfstream,&cycle,sizeof(CountType));
		// read the read length
		bytes += gzread(zfstream,&rlen,sizeof(CountType));
		// read the number of reads
		bytes += gzread(zfstream,&num_reads,sizeof(uint32_t));
		break;
	}
	}

	return bytes;
}

ReadAlignment* iAlnStream::get_alignment() {

  if ( (format==0 && !fstream) || (format==1 && zfstream == Z_NULL) ){
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
    switch (format) {
    case 0:

      fread(buffer.data(),1,buf_size,fstream);
      break;
    case 1:
      gzread(zfstream,buffer.data(),buf_size);
      break;
    case 2:
      lz4read_block();
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
      switch (format) {
      case 0:

        fread(buffer.data(),1,buf_size,fstream);
        break;
      case 1:
        gzread(zfstream,buffer.data(),buf_size);
        break;
      case 2:
	lz4read_block();
	break;
      }
      buf_pos = 0;
    }

  }

  // finally, deserialize the alignment. Set total number of cycles to rlen and increase cycle number by 1.
  ReadAlignment* ra = new ReadAlignment(rlen, cycle+1);

  ra->deserialize(data.data());

  num_loaded++;

  return ra;
}


bool iAlnStream::close() {

  if ( ((format==0 || format==2) && fstream) || (format==1 && zfstream != Z_NULL)) {
    if (num_loaded == num_reads) {
      switch (format) {
      case 0: case 2:
    	  fclose(fstream);
    	  break;
      case 1: gzclose(zfstream); break;
      }
      funlock();
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


void iAlnStream::flock() {
	fileLocks.lock(fname);
	flocked = true;
}


void iAlnStream::funlock() {
	if ( flocked ) {
		fileLocks.unlock(fname);
	}
}




//-------------------------------------------------------------------//
//------  The StreamedAlignment class  ------------------------------//
//-------------------------------------------------------------------//

std::string StreamedAlignment::get_bcl_file(uint16_t cycle, uint16_t mate) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L00" << lane << "/C" << getSeqCycle(cycle, mate) << ".1/s_"<< lane <<"_" << tile << ".bcl";
  return path_stream.str();
}


std::string StreamedAlignment::get_alignment_file(uint16_t cycle, uint16_t mate, std::string base){
  if (base == "") {
    base = globalAlignmentSettings.get_root();
  }
  std::ostringstream path_stream;
  path_stream << base << "/L00" << lane << "/s_"<< lane << "_" << tile << "." << mate << "."<< cycle << ".align";
  return path_stream.str();
}


std::string StreamedAlignment::get_filter_file() {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L00" << lane << "/s_"<< lane << "_" << tile << ".filter";
  return path_stream.str();
}


void StreamedAlignment::create_directories() {
  std::ostringstream path_stream;
  if (globalAlignmentSettings.get_temp_dir() == "") {
    path_stream << globalAlignmentSettings.get_root();
  }
  else {
    path_stream << globalAlignmentSettings.get_temp_dir();
  }
  path_stream << "/L00" << lane;

  boost::filesystem::create_directories(path_stream.str());
  boost::filesystem::create_directories(globalAlignmentSettings.get_out_dir());
}


void StreamedAlignment::init_alignment(uint16_t mate) {
	std::string out_fname = get_alignment_file(0, mate, globalAlignmentSettings.get_temp_dir());

  // get the number of reads in this tile by looking in the first bcl file
  std::string first_cycle = get_bcl_file(1, 0);

  // extract the number of reads
  uint32_t num_reads = num_reads_from_bcl(first_cycle);
  
  // open output alignment stream
  oAlnStream output (lane, tile, 0, rlen, num_reads, globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format());
  output.open(out_fname);

  // write empty read alignments for each read
  for (uint32_t i = 0; i < num_reads; ++i) {
    ReadAlignment * ra = new ReadAlignment(rlen, 0);
    output.write_alignment(ra);
    delete ra;
  }

  if(!output.close()) {
    std::cerr << "Error: Could not create initial alignment file." << std::endl;
  }
} 


uint64_t StreamedAlignment::extend_alignment(uint16_t cycle, uint16_t read_no, uint16_t mate, KixRun* index) {

  // 1. Open the input file
  //-----------------------

  std::string in_fname = get_alignment_file(cycle-1, mate, globalAlignmentSettings.get_temp_dir());
  std::string bcl_fname = get_bcl_file(cycle, read_no);
  std::string filter_fname = get_filter_file();

  iAlnStream input ( globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format() );

  input.open(in_fname);

  assert(input.get_cycle() == cycle-1);
  assert(input.get_lane() == lane);
  assert(input.get_tile() == tile);
  assert(input.get_rlen() == rlen);

  uint32_t num_reads = input.get_num_reads();


  // 2. Open output stream
  //----------------------------------------------------------

  std::string out_fname = get_alignment_file(cycle, mate, globalAlignmentSettings.get_temp_dir());
  oAlnStream output (lane, tile, cycle, rlen, num_reads, globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format());
  output.open(out_fname);


  // 3. Read the full BCL file (this is not too much)
  //-------------------------------------------------

  BclParser basecalls;
  basecalls.open(bcl_fname);


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

//	  iAlnStream input2( globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format() );
//	  input2.open(in_fname);

	  bool testRead = false;


    ReadAlignment* ra = input.get_alignment();

    if (filters.size() > 0 && filters.has_next()) {
      // filter file was found -> apply filter
      if(filters.next()) {
        ra->extend_alignment(basecalls.next(), index, testRead);
        num_seeds += ra->seeds.size();
      }
      else {
        basecalls.next();
        ra->disable();
      }
    }
    // filter file was not found -> treat every alignment as valid
    else {
      ra->extend_alignment(basecalls.next(), index, testRead);
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

  if ( ! ( globalAlignmentSettings.is_keep_aln_files_cycle(getSeqCycle(cycle, globalAlignmentSettings.getSeqByMate(mate).id)-1) || globalAlignmentSettings.is_output_cycle(getSeqCycle(cycle, globalAlignmentSettings.getSeqByMate(mate).id)-1)) ) {
    std::remove(in_fname.c_str());
  }

  return num_seeds;
}


void StreamedAlignment::extend_barcode(uint16_t bc_cycle, uint16_t read_cycle, uint16_t read_no, uint16_t mate) {

	// 1. Open the input file
	//-----------------------

	std::string in_fname = get_alignment_file(read_cycle, mate, globalAlignmentSettings.get_temp_dir());
	std::string bcl_fname = get_bcl_file(bc_cycle, read_no);
	std::string filter_fname = get_filter_file();

	  iAlnStream input ( globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format() );
	  input.open(in_fname);
	  assert(input.get_cycle() == read_cycle);
	  assert(input.get_lane() == lane);
	  assert(input.get_tile() == tile);

	  uint32_t num_reads = input.get_num_reads();


	  // 2. Open output stream
	  //----------------------------------------------------------
	  std::string out_fname = in_fname + ".temp";
	  oAlnStream output (lane, tile, read_cycle, input.get_rlen(), num_reads, globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format());
	  output.open(out_fname);


	  // 3. Read the full BCL file (this is not too much)
	  //-------------------------------------------------
	  BclParser basecalls;
	  basecalls.open(bcl_fname);


	  // 4. Extend barcode sequence
	  //-------------------------------------------------
	  for (uint64_t i = 0; i < num_reads; ++i) {
		char bc = basecalls.next();
	    ReadAlignment* ra = input.get_alignment();
	    ra->appendNucleotideToSequenceStoreVector(bc, true);

	    // filter invalid barcodes if new barcode fragment is completed
	    // TODO: Is done for each mate. Check if it's worth to change it (runtime should not be too high?)
	    if ( !globalAlignmentSettings.get_keep_all_barcodes() && bc_cycle == globalAlignmentSettings.get_seqs()[read_no].length && ra->getBarcodeIndex() == NO_MATCH ) {
	    	ra->disable();
	    }

	    output.write_alignment(ra);
	    delete ra;
	  }

	  // 5. Close files
	  //-------------------------------------------------
	  if (!(input.close() && output.close())) {
	    std::cerr << "Could not finish alignment!" << std::endl;
	  }

	  // 6. Move temp out file to the original file.
	  //-------------------------------------------
	  atomic_rename(out_fname.c_str(), in_fname.c_str());

}


StreamedAlignment& StreamedAlignment::operator=(const StreamedAlignment& other) {
  if(&other == this)
    return *this;

  lane = other.lane;
  tile = other.tile;
  rlen = other.rlen;

  return *this;
}

//-------------------------------------------------------------------//
//------  Streamed SAM generation -----------------------------------//
//-------------------------------------------------------------------//

AlnOut::AlnOut(std::vector<CountType> lns, std::vector<CountType> tls, CountType cycl, KixRun* idx ) : cycle(cycl),index(idx) {

	// Fill list of specified barcodes
	std::vector<std::string> barcodes;
	for ( unsigned i = 0; i < globalAlignmentSettings.get_barcodeVector().size(); i++ ) {
		barcodes.push_back(globalAlignmentSettings.get_barcodeString(i));
	}

	// Get the finished cycles for each mate
	for ( CountType mate = 1; mate <= globalAlignmentSettings.get_mates(); mate++ ) {
		mateCycles.push_back( getMateCycle( mate, cycle ) );
	}

	// Add a waiting task for all lanes and tiles.
	for ( auto ln : lns ) {
		for ( auto tl : tls ) {
			add_task(Task(ln,tl,cycle), WAITING);
		}
	}

};


AlnOut::~AlnOut() {
	if ( !is_finalized() && !finalize() ) {
		std::cerr << "Could not finish output for cycle " << cycle << "." << std::endl;
	}
}


void AlnOut::init() {

	std::lock_guard<std::mutex> lock(if_lock);

	if ( initialized )
		return;

	// Init the bamIOContext (the same object can be used for all output streams)
	bfos.set_context(index->getSeqNames(), index->getSeqLengths());

	// Init the header (the same object can be used for all output streams)
	seqan::BamHeader header = getBamHeader();

	// Init output stream for each barcode (plus undetermined if keep_all_barcodes is set)
	for ( unsigned barcode=0; barcode < (barcodes.size() + 1); barcode ++) {
		if ( barcode < barcodes.size() || globalAlignmentSettings.get_keep_all_barcodes() ) {

			std::string barcode_string = ( barcode == barcodes.size() ) ? "undetermined" : barcodes[barcode];

			// Open file in Bam output stream and write the header
			bfos.emplace_back( getBamTempFileName(barcode_string, cycle).c_str() );
			bfos[barcode].writeHeader(header);

		}
	}

	initialized = true;
}


bool AlnOut::set_task_status( Task t, ItemStatus status ) {
	std::lock_guard<std::mutex> lock(tasks_lock);
	if ( tasks.find(t) == tasks.end() )
		return false;
	tasks[t] = status;
	return true;
}


bool AlnOut::set_task_status_from_to( Task t, ItemStatus oldStatus, ItemStatus newStatus ) {
	std::lock_guard<std::mutex> lock(tasks_lock);
	if ( tasks.find(t) == tasks.end() )
		return false;
	if ( tasks[t] == oldStatus ) {
		tasks[t] = newStatus;
		return true;
	}
	return false;
}


Task AlnOut::get_next( ItemStatus getStatus, ItemStatus setToStatus ) {
	std::lock_guard<std::mutex> lock(tasks_lock);
	for ( auto it = tasks.begin(); it != tasks.end(); ++it ) {
		if ( it->second == getStatus ) {
			tasks[it->first] = setToStatus;
			return it->first;
		}
	}
	return NO_TASK;
}


bool AlnOut::add_task( Task t, ItemStatus status ) {
	std::lock_guard<std::mutex> lock(tasks_lock);
	if ( tasks.find(t) != tasks.end() )
		return false;
	tasks[t] = status;
	return true;
}


bool AlnOut::sort_tile( CountType ln, CountType tl, CountType mate, CountType cycle, bool overwrite ) {


	std::string in_fname = alignment_name(ln, tl, cycle, mate);
	std::string out_fname = alignment_name(ln, tl, cycle, mate) + ".sorted";

	// Stop if sorted file already exist
	if ( file_exists( out_fname ) && !overwrite )
		return true;

	iAlnStream input ( globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format() );

	input.open(in_fname);

	assert(input.get_cycle() == cycle);
	assert(input.get_lane() == ln);
	assert(input.get_tile() == tl);

	uint32_t num_reads = input.get_num_reads();

	oAlnStream output(ln, tl, cycle, input.get_rlen(), num_reads, globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format());
	output.open(out_fname);

	for ( uint32_t i = 0; i < num_reads; i++ ) {

		try {
			ReadAlignment * ra = input.get_alignment();
			ra->sort_seeds_by_as();
			output.write_alignment(ra);
			delete ra;
		} catch (const std::exception & ex) {
			return false;
		}

	}

	if (!(input.close() && output.close()))
		return false;

	return true;

}


void AlnOut::write_tile_to_bam ( Task t ) {

	if ( !is_initialized() )
		init();

	try {
		__write_tile_to_bam__ (t);
		set_task_status( t, FINISHED );
	} catch ( const std::exception& e) {
		set_task_status( t, FAILED );
		std::cerr << "Writing of task " << t << " failed: " << e.what() << std::endl;
	}
}


void AlnOut::__write_tile_to_bam__ ( Task t) {

	std::vector<std::lock_guard<std::mutex>> locks;

	CountType lane = t.lane;
	CountType tile = t.tile;

	////////////////////////////////////////////////////
	//  Main loop //////////////////////////////////////
	////////////////////////////////////////////////////

	// set the filter file
	std::string filter_fname = filter_name(lane, tile);
	FilterParser filters;
	if (file_exists(filter_fname)) {
		filters.open(filter_fname);
	}


	// set the alignment files
	std::vector<iAlnStream*> alignmentFiles;
	unsigned numberOfAlignments = 0;

	for (unsigned mateIndex = 1; mateIndex <= mateCycles.size(); mateIndex++) {

		if ( globalAlignmentSettings.getSeqByMate(mateIndex) == NULLSEQ )
			return;

		CountType mateCycle = mateCycles[mateIndex-1];

		if ( !sort_tile( lane, tile, mateIndex, mateCycle, globalAlignmentSettings.get_force_resort()) ) {
			continue;
		}

		// Open sorted alignment file
		std::string alignment_fname = alignment_name(lane, tile, mateCycle, mateIndex) + ".sorted";
		if ( !file_exists(alignment_fname) ) {
			continue;
		}

		iAlnStream* input = new iAlnStream( globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format() );
		input->open(alignment_fname);

		// compare number of reads in alignment file with number of reads in filter file, if filter file exists
		if (file_exists(filter_fname) && input->get_num_reads() != filters.size()) {
			throw std::length_error("Unequal number of reads (.filer vs .align)");
		}

		// compare number of reads in alignment file with number of reads in previous alignment file
		if (mateIndex != 1 && input->get_num_reads() != numberOfAlignments) {
			throw std::length_error("Unequal number of reads (between mates)");
		}

		numberOfAlignments = input->get_num_reads(); // set this after last if-then construct
		alignmentFiles.push_back(input);
	}

	// for all reads in a tile
	/////////////////////////////////////////////////////////////////////////////
	for (uint64_t i = 0; i < numberOfAlignments; i++) {

		std::vector<seqan::BamAlignmentRecord> records;

		std::vector<ReadAlignment*> mateAlignments;
		for (auto e:alignmentFiles) {
			mateAlignments.push_back(e->get_alignment());
		}

		// if the filter file is available and the filter flag is 0 then skip
		if (filters.size() != 0 && filters.next() == false)
			continue;

		// compute barcode sequence as it should be written to BC tag
		std::string barcode = globalAlignmentSettings.format_barcode(mateAlignments[0]->getBarcodeString());

		// Barcode index for the read
		CountType barcodeIndex = mateAlignments[0]->getBarcodeIndex();

		// If read has undetermined barcode and keep_all_barcodes is not set, skip this read
		if ( barcodeIndex == NO_MATCH && !globalAlignmentSettings.get_keep_all_barcodes() )
			continue;
		else if ( barcodeIndex == NO_MATCH )
			barcodeIndex = barcodes.size(); // this is the index for the "undetermined" output stream

		// setup QNAME
		// Read name format <instrument‐name>:<run ID>:<flowcell ID>:<lane‐number>:<tile‐number>:<x‐pos>:<y‐pos>
		// readname << "<instrument>:<run-ID>:<flowcell-ID>:" << ln << ":" << tl << ":<xpos>:<ypos>:" << i;
		//TODO: where do we get the Illumina read coordinate from?
		std::stringstream readname;
		readname << "lane." << lane << "|tile." << tile << "|read." << i;

		// Track equivalent alignments for the same mate (similar positions from different seeds -> only keep the best one)
		// TODO: implement equivalent alignment window as user parameter
		PositionType equivalentAlignmentWindow = 10;
		std::set<PositionType> alignmentPositions;


		// for all mates
		/////////////////////////////////////////////////////////////////////////////
		for (unsigned mateAlignmentIndex=0; mateAlignmentIndex < mateAlignments.size(); ++mateAlignmentIndex) {

			// Alignment disabled
			if ( mateAlignments[mateAlignmentIndex]->is_disabled() )
				continue;

			// No Seeds
			if ( mateAlignments[mateAlignmentIndex]->seeds.size() == 0 )
				continue;

			readname << "|r" << mateAlignmentIndex+1;

			// Variables for output modes
			ScoreType first_seed_score = 0;
			ScoreType last_seed_score = 0;
			CountType num_diff_scores = 0;

			// Number of printed alignments for the current mate.
			unsigned printedMateAlignments = 0;

			// Unique mode interruption
			if ( mateAlignments[mateAlignmentIndex]->seeds.size() > 1 && globalAlignmentSettings.get_unique_hit_mode() )
				continue; // TODO: handle if a single seed has more than one position, maybe just use FM index span

			std::string seq = mateAlignments[mateAlignmentIndex]->getSequenceString();
			std::string qual = mateAlignments[mateAlignmentIndex]->getQualityString();

			std::vector<uint8_t> mapqs = mateAlignments[mateAlignmentIndex]->getMAPQs();
			auto mapqs_it = mapqs.begin();

			// for all seeds
			/////////////////////////////////////////////////////////////////////////////
			for (SeedVecIt it = mateAlignments[mateAlignmentIndex]->seeds.begin(); it != mateAlignments[mateAlignmentIndex]->seeds.end(); ++it, ++mapqs_it) {

				ScoreType curr_seed_score = (*it)->get_as();

				// If no alignment was printed before, the current one has the best "score"
				if ( printedMateAlignments == 0 )
					first_seed_score = curr_seed_score;

				// Stop in all best mode when AS:i score is lower than the first
				if( globalAlignmentSettings.get_all_best_hit_mode() && first_seed_score > curr_seed_score )
					goto nextmate;


				// Don't write this seed if the user-specified score or softclip ratio is not fulfilled
				CountType softclip_length = (*it)->get_softclip_length();
				if ( curr_seed_score < globalAlignmentSettings.get_min_as() || softclip_length > globalAlignmentSettings.get_max_softclip_ratio()*mateCycles[mateAlignmentIndex]) {
					continue;
				}

				// get CIGAR-String
				seqan::String<seqan::CigarElement<> > cigar = (*it)->returnSeqanCigarString();

				// Get NM:i value
				unsigned nm = (*it)->get_nm();


				// check if cigar string sums up to read length
				// TODO Potentially conflicts with the 'eachMateAligned' flag if done here.
				unsigned cigarElemSum = 0;
				unsigned deletionSum = 0;
				unsigned supposed_cigar_length = mateCycles[mateAlignmentIndex];

				for (seqan::Iterator<seqan::String<seqan::CigarElement<> > >::Type elem = seqan::begin(cigar); elem != end(cigar); ++elem) {
					if ((elem->operation == 'M') || (elem->operation == 'I') || (elem->operation == 'S') || (elem->operation == '=') || (elem->operation == 'X'))
						cigarElemSum += elem->count;

					if (elem->operation == 'D') {
						deletionSum += elem->count;
					}
				}
				if (cigarElemSum != supposed_cigar_length) {
					it = mateAlignments[mateAlignmentIndex]->seeds.erase(it);
					continue;
				}
				if (deletionSum >= supposed_cigar_length) {
					it = mateAlignments[mateAlignmentIndex]->seeds.erase(it);
					continue;

				}

				// Get positions for the current seed
				PositionPairListType pos_list;

				mateAlignments[mateAlignmentIndex]->getPositions(index, *it, pos_list);

				// handle all positions
				for ( auto p = pos_list.begin(); p != pos_list.end(); ++p ) {

					// Stop in any best mode when first alignment was already written
					if( globalAlignmentSettings.get_any_best_hit_mode() && printedMateAlignments >= 1 )
						goto nextmate;

					seqan::BamAlignmentRecord record;

					record.qName = readname.str();

					record.rID = CountType(p->first / 2);

					record.beginPos = mateAlignments[mateAlignmentIndex]->get_SAM_start_pos(index, *p, *it);

					record.mapQ = *mapqs_it;

					// skip invalid positions
					if (record.beginPos < 0 || PositionType(record.beginPos) == std::numeric_limits<PositionType>::max()) {
						continue;
					}

					// skip positions that were already written (equivalent alignments). This can be done because the best alignment for this position is written first.
					if ( alignmentPositions.find(record.beginPos - ( record.beginPos % equivalentAlignmentWindow ) ) != alignmentPositions.end() ||
							alignmentPositions.find(record.beginPos +  (equivalentAlignmentWindow - ( record.beginPos % equivalentAlignmentWindow ) ) ) != alignmentPositions.end()) {
						continue;
					}

					record.cigar = cigar;
					if ( index->isReverse(p->first) )
						seqan::reverse(record.cigar);

					// flag and seq
					record.flag = 0;
					record.seq = seq;
					record.qual = qual;

					if ( printedMateAlignments > 0 ) { // if current seed is secondary alignment
						record.flag |= 256;
						seqan::clear(record.seq);
						seqan::clear(record.qual);
					}

					if ( index->isReverse(p->first) ) { // if read matched reverse complementary
						seqan::reverseComplement(record.seq);
						seqan::reverse(record.qual);
						record.flag |= 16;
					}

					if (globalAlignmentSettings.get_mates() > 1) { // if there are more than two mates
						record.flag |= 1;
						if (mateAlignmentIndex == 0) {
							record.flag |= 64;
						} else if (mateAlignmentIndex == mateAlignments.size()-1) {
							record.flag |= 128;
						} else {
							record.flag |= 192; // 64 + 128
						}

						bool eachMateAligned = true;
						for (auto & e : mateAlignments)
							eachMateAligned = eachMateAligned && e->seeds.size() > 0;
						if (eachMateAligned)
							record.flag |= 2;
					}


					// tags
					seqan::BamTagsDict dict;

					seqan::appendTagValue(dict, "AS", curr_seed_score);

					if (barcode!="")
						seqan::appendTagValue(dict, "BC", barcode);

					seqan::appendTagValue(dict, "NM", nm);

					std::string mdz = (*it)->getMDZString();

					if ( index->isReverse(p->first))
						mdz = reverse_mdz(mdz);

					seqan::appendTagValue(dict, "MD", mdz);


					record.tags = seqan::host(dict);



					// fill records list
					records.push_back(record);

					// set variables for mode selection
					if ( last_seed_score != curr_seed_score || num_diff_scores == 0 )
						++num_diff_scores;
					last_seed_score = curr_seed_score;

					++printedMateAlignments;

				}
			}
			nextmate: {};
		}

		// Write all records as a group to keep suboptimal alignments and paired reads together.
		bfos[barcodeIndex].writeRecords(records);

		for (auto e:mateAlignments)
			delete e;
	}
	for (auto e:alignmentFiles)
		delete e;

	return;
}


Task AlnOut::write_next ( ) {
	Task t = get_next ( BCL_AVAILABLE, RUNNING );
	if ( t != NO_TASK ) {
		write_tile_to_bam ( t );
		return t;
	} else {
		return NO_TASK;
	}
}


CountType AlnOut::get_task_status_num ( ItemStatus getStatus ) {
	CountType num = 0;
	std::lock_guard<std::mutex> lock(tasks_lock);
	for ( auto it = tasks.begin(); it != tasks.end(); ++it ) {
		if ( it->second == getStatus ) {
			num += 1;
		}
	}
	return num;
}


bool AlnOut::finalize() {

	std::lock_guard<std::mutex> lock(if_lock);

	if ( finalized )
		return true;

	// Don't finish if there are unfinished tasks.
	if ( !is_finished() )
		return false;

	bool success = true;

	bfos.clear();

	// Move all output files to their final location.
	for ( unsigned barcode=0; barcode < barcodes.size() + 1; barcode ++) {
		if ( barcode < barcodes.size() || globalAlignmentSettings.get_keep_all_barcodes() ) {

			std::string barcode_string = ( barcode == barcodes.size() ) ? "undetermined" : barcodes[barcode];

			int rename = atomic_rename(getBamTempFileName(barcode_string, cycle).c_str(), getBamFileName(barcode_string, cycle).c_str());
			if ( rename == -1 ) {
				std::cerr << "Renaming temporary output file " << getBamTempFileName(barcode_string, cycle).c_str() << " to " << getBamFileName(barcode_string, cycle).c_str() << " failed." << std::endl;
				success = false;
			}
		}
	}

	// If it comes here, it counts as finalized independently from the success state (otherwise, contradictory error messages may occur).
	finalized = true;

	return success;
}
