#include "tools_static.h"

///////////////////////////////////
////////// File handling //////////
///////////////////////////////////

std::vector<char> read_binary_file(const std::string &fname) {

  // get file size
  uint64_t size = get_filesize(fname);

  // open binary file
  FILE* f;
  f = fopen(fname.c_str(), "rb");

  if (!f) {
	  std::stringstream error;
	  error << "Error reading binary file " << fname << ": Could not open file.";
	  throw std::ios::failure ( error.str() );
  }

  // allocate memory
  std::vector<char> data (size);

  // read all data at once
  uint64_t read = fread(data.data(), 1, size, f);

  if (read != size){
	  std::stringstream error;
	  error << "Error reading binary file " << fname << ": Read " << read << " bytes while file has " << size << " bytes.";
	  throw std::ios::failure ( error.str() );
  }

  fclose(f);

  return data;
}

uint64_t write_binary_file(const std::string &fname, const std::vector<char> & data) {

  // open binary file
  FILE* ofile;
  ofile = fopen(fname.c_str(), "wb");

  if (!ofile) {
	  std::stringstream error;
	  error << "Error serializing object to file " << fname << ": Could not open file for writing.";
	  throw std::ios::failure ( error.str() );
  }

  // write all data
  uint64_t written = fwrite(data.data(), 1, data.size(), ofile);

  // close file
  fclose(ofile);

  if (written != data.size()){
	  std::stringstream error;
	  error << "Error serializing object to file " << fname << ": Wrote " << written << " bytes while data contains " << data.size() << " bytes.";
	  throw std::ios::failure ( error.str() );
  }

  return written;
}

std::string get_file_suffix ( OutputFormat format ) {
	switch ( format ) {
	case OutputFormat::SAM:
		return ".sam";
		break;
	case OutputFormat::BAM:
		return ".bam";
		break;
	case OutputFormat::CRAM:
		return ".cram";
		break;
	default:
		return ".txt";
		break;
	}
}


////////////////////////////////////////////////
////////// Property trees / XML files //////////
////////////////////////////////////////////////

bool read_xml(boost::property_tree::ptree & xml_in, std::string xml_fname) {

	if ( !file_exists(xml_fname) ) {
		std::cout << "XML file not found: " << xml_fname << std::endl;
		return false;
	}

	try {
		boost::property_tree::read_xml (xml_fname, xml_in);
	} catch ( const std::exception &ex) {
		std::cerr << "Error loading xml file " << xml_fname << ": " << std::endl << ex.what() << std::endl;
		return false;
	}

	return true;

}

bool write_ini(boost::property_tree::ptree & ini_out, std::string ini_fname) {

	try {
		boost::property_tree::write_ini( ini_fname, ini_out );
	} catch ( const std::exception &ex ) {
		std::cerr << "Error writing config file " << ini_fname << ": " << std::endl << ex.what() << std::endl;
		return false;
	}

	return true;

}


/////////////////////////////////
////////// Other stuff //////////
/////////////////////////////////

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
  if (!res) {
    std::cerr << "Error extracting number of reads from BCL file " << bcl << std::endl;
    return 0;
  }

  // close file
  fclose (ifile);

  return num_reads;
}

AlignmentMode to_alignmentMode ( std::string value ) {
	if ( std::toupper(value[0]) == 'B' )
		return BALANCED;
	if ( std::toupper(value[0]) == 'A' )
		return ACCURATE;
	if ( std::toupper(value[0]) == 'F' )
		return FAST;
	if ( std::toupper(value[0]) == 'V' ) {
		if ( std::toupper(value[1]) == 'F' || std::toupper(value[5]) == 'F')
			return VERYFAST;
		if ( std::toupper(value[1]) == 'A' || std::toupper(value[5]) == 'A')
			return VERYACCURATE;
	}
	throw std::runtime_error("Invalid alignment mode " + value + ".");
}
