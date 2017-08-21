#include "tools_static.h"


/////////////////////////////////
////////// Comparators //////////
/////////////////////////////////

bool gp_compare (GenomePosType i,GenomePosType j) {
	if ( i.pos == j.pos )
		return i.gid < j.gid;
	return (i.pos < j.pos);
}


/////////////////////////////////////
////////// Type convertion //////////
/////////////////////////////////////

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}


///////////////////////////////////
////////// File handling //////////
///////////////////////////////////

std::ifstream::pos_type get_filesize(const std::string &fname)
{
  std::ifstream in(fname, std::ios::binary | std::ios::ate);
  return in.tellg();
}

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

bool file_exists(const std::string &fname) {
  return boost::filesystem::exists(fname);

}

std::string absolute_path(std::string fname) {
	boost::filesystem::path input_path(fname);
	return boost::filesystem::canonical(fname).string();
}

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

uint64_t write_binary_file(const std::string &fname, const std::vector<char> & data) {

  // open binary file
  FILE* ofile;
  ofile = fopen(fname.c_str(), "wb");

  if (!ofile) {
    std::cerr << "Error serializing object to file " << fname << ": Could not open file for writing." << std::endl;
    return 1;
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

bool write_xml(boost::property_tree::ptree & xml_out, std::string xml_fname) {

	try {
		boost::property_tree::write_xml( xml_fname, xml_out );
	} catch ( const std::exception &ex ) {
		std::cerr << "Error writing xml file " << xml_fname << ": " << std::endl << ex.what() << std::endl;
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

