#include "alnstream.h"


//-------------------------------------------------------------------//
//------  The output Alignment Stream class  ------------------------//
//-------------------------------------------------------------------//

// new output Alignment Stream class
oAlnStream::oAlnStream(uint16_t ln, uint16_t tl, uint16_t cl, CountType rl, uint32_t nr, uint64_t bs, uint8_t fmt):
  lane(ln), tile(tl), cycle(cl), rlen(rl), num_reads(nr), num_written(0), buffer(bs,0), buf_size(bs), buf_pos(0), format(fmt), ofile(NULL), ozfile(Z_NULL) {}


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
    if(written != buf_size)
      throw std::runtime_error("Could not write out buffer 1 in oAlnStream::write_alignment.");

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
      if(written != buf_size)
        throw std::runtime_error("Could not write out buffer 2 in oAlnStream::write_alignment.");
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
    if(written != buf_pos)
      throw std::runtime_error("Could not write out buffer in oAlnStream::close.");
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
  lane(0), tile(0), cycle(0), rlen(0), num_reads(0), num_loaded(0), buffer(bs,0), buf_size(bs), buf_pos(bs), format(fmt), ifile(NULL), izfile(Z_NULL) {}


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
  switch (format) {
  case 0: case 2:
    {
      // read the lane
      bytes += fread(&lane,sizeof(uint16_t),1,ifile);
      // read the tile
      bytes += fread(&tile,sizeof(uint16_t),1,ifile);
      // read the cycle
      bytes += fread(&cycle,sizeof(CountType),1,ifile);
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
    switch (format) {
    case 0:
      fread(buffer.data(),1,buf_size,ifile);
      break;
    case 1:
      gzread(izfile,buffer.data(),buf_size);    
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
        fread(buffer.data(),1,buf_size,ifile);
        break;
      case 1:
        gzread(izfile,buffer.data(),buf_size);    
        break;
      case 2:
        lz4read_block();    
        break;
      }
      buf_pos = 0;
    }
  }

  // finally, deserialize the alignment
  ReadAlignment* ra = new ReadAlignment();
  ra->set_total_cycles(rlen);
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
  rlen = other.rlen;
  
  return *this;
}

std::string StreamedAlignment::get_bcl_file(uint16_t cycle, uint16_t read_number) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L00" << lane << "/C" << getSeqCycle(cycle, read_number) << ".1/s_"<< lane <<"_" << tile << ".bcl";
  return path_stream.str();
}


// get the path to the alignment file. The alignment file is located in
// <base>/L00<lane>/s_<lane>_<tile>.<cycle>.align
// if base == "": base = globalAlignmentSettings.get_root()
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


// create directories required to store the alignment files (only if not stored in globalAlignmentSettings.get_root())
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
  boost::filesystem::create_directories(globalAlignmentSettings.get_out_dir().string());
}


// initialize empty alignment. Creates alignment files for a virtual Cycle 0
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
    ReadAlignment * ra = new ReadAlignment();
    ra->set_total_cycles(rlen);
    output.write_alignment(ra);
    delete ra;
  }

  if(!output.close()) {
    std::cerr << "Error: Could not create initial alignment file." << std::endl;
  }
} 



// extend an existing alignment from cycle <cycle-1> to <cycle>. returns the number of seeds
uint64_t StreamedAlignment::extend_alignment(uint16_t cycle, uint16_t read_no, uint16_t mate, KixRun* index) {

  // 1. Open the input file
  //-----------------------
  std::string in_fname = get_alignment_file(cycle-1, mate, globalAlignmentSettings.get_temp_dir());
  std::string bcl_fname = get_bcl_file(cycle, read_no); // TODO: correct cycle
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
  if ( ! ( globalAlignmentSettings.get_keep_aln_files() || globalAlignmentSettings.is_output_cycle(cycle-1)) ) {
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
	  std::rename(out_fname.c_str(), in_fname.c_str());

}


//-------------------------------------------------------------------//
//------  Streamed SAM generation -----------------------------------//
//-------------------------------------------------------------------//

uint64_t alignments_to_sam(std::vector<uint16_t> lns, std::vector<uint16_t> tls, KixRun* index, CountType cycle) {

	std::ofstream logfile;
	std::string logfile_fname = globalAlignmentSettings.get_out_dir().string() + "/hilive_out.log";
	logfile.open( logfile_fname, std::ofstream::app );

	logfile << "Start to write output for cycle " << std::to_string(cycle) << "." << std::endl;

	std::string file_suffix = globalAlignmentSettings.get_write_bam() ? ".bam" : ".sam";
	std::string file_cycle = cycle >= globalAlignmentSettings.get_cycles() ? "" : "_cycle" + std::to_string(cycle);

	// Fill list of specified barcodes
	std::vector<std::string> barcodes;

	// Add user-specified barcode strings
	for ( unsigned i = 0; i < globalAlignmentSettings.get_barcodeVector().size(); i++ ) {
		barcodes.push_back(globalAlignmentSettings.get_barcodeString(i));
	}

	// Mates to handle
	CountType mates = 0;
	CountType last_mate_cycle;

	auto seqs = globalAlignmentSettings.get_seqs();
	auto seqs_it = seqs.begin();

	for ( CountType cycles_left = cycle; cycles_left > 0;) {

		// add mate if element is no barcode
		if ( seqs_it->mate != 0 )
			mates += 1;

		// reduce number of remaining cycles. If not enough cycles left, set num cycles of the last mate and break.
		if ( cycles_left >= seqs_it->length )
			cycles_left -= seqs_it->length;
		else {
			last_mate_cycle = cycles_left;
			break;
		}

		// go to next Seq element. If the last one is processed, add length of the last mate and break.
		if ( ++seqs_it == seqs.end() ) {
			last_mate_cycle = globalAlignmentSettings.getSeqByMate(mates).length;
			break;
		}
	}

	// Init the bamIOContext (the same object can be used for all output streams)
	seqan::StringSet<seqan::CharString> referenceNames;
	seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > referenceNamesCache(referenceNames);
	seqan::BamIOContext<seqan::StringSet<seqan::CharString> > bamIOContext(referenceNames, referenceNamesCache);

	seqan::contigNames(bamIOContext) = index->seq_names;
	seqan::contigLengths(bamIOContext) = index->seq_lengths;

	// Init the header (the same object can be used for all output streams)
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


	// Vector that contains the output streams (USE POINTERS !)
	std::vector<std::unique_ptr<seqan::BamFileOut>> bfos;

	// Init output stream for each barcode (plus undetermined if keep_all_barcodes is set)
	for ( unsigned barcode=0; barcode < barcodes.size() + 1; barcode ++) {
		if ( barcode < barcodes.size() || globalAlignmentSettings.get_keep_all_barcodes() ) {
			std::string out_fname;
			std::string barcode_string = barcode < barcodes.size() ? barcodes[barcode] : "undetermined";
			out_fname = globalAlignmentSettings.get_out_dir().string() + "/hilive_out_" + barcode_string + file_cycle + ".temp" + file_suffix;
			std::unique_ptr<seqan::BamFileOut> bfo( new seqan::BamFileOut(out_fname.c_str()));
			bfos.push_back( std::move(bfo) );
			bfos[barcode]->context = bamIOContext;
			seqan::writeHeader(*bfos[barcode], header);
		}
	}


	////////////////////////////////////////////////////
	//  Main loop //////////////////////////////////////
	////////////////////////////////////////////////////

	uint64_t num_alignments = 0;
	unsigned totalNumberOfReads = 0;

	// for all lanes
	/////////////////////////////////////////////////////////////////////////////
	for (auto ln:lns) {

		// for all tiles
		/////////////////////////////////////////////////////////////////////////////
		for (auto tl:tls) {

			// set the filter file
			std::string filter_fname = filter_name(ln, tl);
			FilterParser filters;
			if (file_exists(filter_fname)) {
				filters.open(filter_fname);
			}
			else
				logfile << "Could not find .filter file: " <<  filter_fname  << ". Treated all reads as valid."<< std::endl;


			// set the alignment files
			std::vector<iAlnStream*> alignmentFiles;
			unsigned numberOfAlignments = 0;
			for (unsigned mateIndex = 1; mateIndex <= mates; ++mateIndex) {
				if ( globalAlignmentSettings.getSeqByMate(mateIndex) == NULLSEQ ) return 0;

				CountType mateCycle = mateIndex==mates ? last_mate_cycle : globalAlignmentSettings.getSeqByMate(mateIndex).length;

				std::string alignment_fname = alignment_name(ln, tl, mateCycle, mateIndex);
				if ( !file_exists(alignment_fname) ) {
					logfile << "Alignment file not found: " << alignment_fname << ". Ignored all related alignments." << std::endl;
					continue;
				}
				iAlnStream* input = new iAlnStream( globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format() );
				input->open(alignment_fname);

				// compare number of reads in alignment file with number of reads in filter file, if filter file exists
				if (file_exists(filter_fname) && input->get_num_reads() != filters.size()) {
					std::string msg = std::string("Number of reads in filter file (") + std::to_string(filters.size()) + ") does not match the number of reads in the alignment file (" + std::to_string(input->get_num_reads()) + ").";
					logfile << "Writing output for cycle " << std::to_string(cycle) << "not successful: " << msg.c_str() << std::endl;
					throw std::length_error(msg.c_str());
				}

				// compare number of reads in alignment file with number of reads in previous alignment file
				if (mateIndex != 1 && input->get_num_reads() != numberOfAlignments) {
					std::string msg = std::string("Number of reads in alignment file (") + std::to_string(input->get_num_reads()) + ") does not match the number of reads in previous alignment file (" + std::to_string(input->get_num_reads()) + ").";
					logfile << "Writing output for cycle " << std::to_string(cycle) << "not successful: " << msg.c_str() << std::endl;
					throw std::length_error(msg.c_str());
				}

				numberOfAlignments = input->get_num_reads(); // set this after last if-then construct
				alignmentFiles.push_back(input);
			}
			totalNumberOfReads += numberOfAlignments;

			// for all reads in a tile
			/////////////////////////////////////////////////////////////////////////////
			for (uint64_t i = 0; i < numberOfAlignments; ++i) {

				std::vector<ReadAlignment*> mateAlignments;
				for (auto e:alignmentFiles) {
					mateAlignments.push_back(e->get_alignment());
				}

				// if the filter file is available and the filter flag is 0 then skip
				if (filters.size() != 0 && filters.next() == false)
					continue;

				// compute barcode sequence as it should be written to BC tag
				std::string barcode = mateAlignments[0]->getBarcodeString(); // barcode how HiLive read it from .bcl files
				if (barcode!="") { // if demultiplexing is on
					// insert "-" as delimiter between the single barcodes
					barcode = globalAlignmentSettings.format_barcode(barcode);
				}

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
				readname << "lane." << ln << "|tile." << tl << "|read." << i;


				// for all mates
				/////////////////////////////////////////////////////////////////////////////
				seqan::BamAlignmentRecord record;
				unsigned printedMates = 0;
				for (unsigned mateAlignmentIndex=0; mateAlignmentIndex < mateAlignments.size(); ++mateAlignmentIndex) {

					// for all seeds
					/////////////////////////////////////////////////////////////////////////////
					for (SeedVecIt it = mateAlignments[mateAlignmentIndex]->seeds.begin(); it != mateAlignments[mateAlignmentIndex]->seeds.end(); ) {
						if ( (*it)->gid == TRIMMED ) {
							if (mateAlignments[mateAlignmentIndex]->seeds.size() == 1) {
								globalAlignmentSettings.add_trimmedRead(i);
							}
							++it;
							continue;
						}

						seqan::clear(record);

						record.qName = readname.str();

						record.rID = (*it)->gid;

						record.beginPos = mateAlignments[mateAlignmentIndex]->get_SAM_start_pos(*it)-1; // seqan expects 0-based positions, but in HiLive we use 1-based
						if (record.beginPos < 0) {
							it = mateAlignments[mateAlignmentIndex]->seeds.erase(it);
							continue;
						}

						unsigned nm_i = 0;
						record.cigar = (*it)->returnSeqanCigarString(&nm_i);

						// flag and seq
						record.flag = 0;
						record.seq = mateAlignments[mateAlignmentIndex]->getSequenceString();
						if ((*it)->start_pos < 0) { // if read matched reverse complementary
							seqan::reverseComplement(record.seq);
							record.flag |= 16;
						}
						if (printedMates >= mateAlignments.size()) { // if current seed is secondary alignment
							record.flag |= 256;
							seqan::clear(record.seq);
							seqan::clear(record.qual);
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
							for (auto e:mateAlignments)
								eachMateAligned = eachMateAligned && e->seeds.size() > 0;
							if (eachMateAligned)
								record.flag |= 2;
						}


						// check if cigar string sums up to read length
						// TODO
						// Jakob: I have not seen such a warning in a long time and a correct algorithm should prevent these cases anyway.
						// Furthermore, this is a filtering step and potentially conflicts with the 'eachMateAligned' flag if done here.
						// However, Tobi has done something here
						unsigned cigarElemSum = 0;
						unsigned deletionSum = 0;
						unsigned supposed_cigar_length = mateAlignmentIndex == CountType(mates-1) ? last_mate_cycle : globalAlignmentSettings.getSeqByMate(mateAlignmentIndex + 1).length;
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
						if (cigarElemSum != supposed_cigar_length) {
							logfile << "WARNING: Excluded an alignment of read " << record.qName << " at position " << mateAlignments[mateAlignmentIndex]->get_SAM_start_pos(*it) << " because its cigar vector had length " << cigarElemSum << std::endl;
							it = mateAlignments[mateAlignmentIndex]->seeds.erase(it);
							continue;
						}
						if (deletionSum >= supposed_cigar_length) {
							logfile << "WARNING: Excluded an alignment of read " << record.qName << " at position " << mateAlignments[mateAlignmentIndex]->get_SAM_start_pos(*it) << " because its cigar vector had " << deletionSum << " deletions" << std::endl;
							it = mateAlignments[mateAlignmentIndex]->seeds.erase(it);
							continue;
						}

						// tags
						seqan::BamTagsDict dict;
						seqan::appendTagValue(dict, "AS", ( asi_score ) );
						if (barcode!="") { // if demultiplexing is on
							seqan::appendTagValue(dict, "BC", barcode);
						}
						seqan::appendTagValue(dict, "NM", nm_i);
						record.tags = seqan::host(dict);


						// write record to disk
						seqan::writeRecord(*bfos[barcodeIndex], record);

						++num_alignments;
						++printedMates;
						++it;
					}
				}

				for (auto e:mateAlignments)
					delete e;
			}
			for (auto e:alignmentFiles)
				delete e;
		}

	}

	// Init output stream for each barcode (plus undetermined if keep_all_barcodes is set)
	for ( unsigned barcode=0; barcode < barcodes.size() + 1; barcode ++) {
		if ( barcode < barcodes.size() || globalAlignmentSettings.get_keep_all_barcodes() ) {
			std::string barcode_string = barcode < barcodes.size() ? barcodes[barcode] : "undetermined";
			std::string old_out_fname = globalAlignmentSettings.get_out_dir().string() + "/hilive_out_" + barcode_string + file_cycle + ".temp" + file_suffix;
			std::string new_out_fname = globalAlignmentSettings.get_out_dir().string() + "/hilive_out_" + barcode_string + file_cycle + file_suffix;
			std::rename(old_out_fname.c_str(), new_out_fname.c_str());
		}
	}

	// TODO maybe find a way to generate statsfiles when generating multiple output files.
	std::ofstream statsfile;
	std::string statsfile_fname = globalAlignmentSettings.get_out_dir().string() + "/hilive_out";
	statsfile_fname += cycle == globalAlignmentSettings.get_cycles() ? ".stats" : "_cycle" + std::to_string(cycle) + ".stats";
	statsfile.open( statsfile_fname );
	statsfile << "Number of reads\t" << totalNumberOfReads << std::endl;
	statsfile << "Number of alignments\t" << num_alignments << std::endl;
	statsfile.close();

	// TODO What todo with completely trimmed reads?
	// --> Amount could for example be written in the stats file(s)

	logfile << "Writing output for cycle " << std::to_string(cycle) << " finished." << std::endl;
	logfile.close();

	return 1;

}
