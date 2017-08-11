/**
 * This class provides functions that are dependent of the alignmentSettings!
 * For functions that are not dependent of any other HiLive class, please use the tools_static class!
 * Please do NOT add further includes to this file since this will lead to unwanted dependencies!
 */

#ifndef TOOLS_STATIC_H
#define TOOLS_STATIC_H

/* DONT ADD ANY INCLUDES */
#include "headers.h"
#include "definitions.h"
/* DONT ADD ANY INCLUDES */


/** Compare function to sort GenomePosType objects by position. */
bool gp_compare (GenomePosType i,GenomePosType j);

/** Extract the number of reads from a BCL file. */
uint32_t num_reads_from_bcl(std::string bcl);

/** Split a std::string to a std::vector<std::string>. */
void split(const std::string &s, char delim, std::vector<std::string> &elems);


///////////////////////////////////////////////
// General File Handling
///////////////////////////////////////////////

/** Get total size of a file (in bytes) */
std::ifstream::pos_type get_filesize(const std::string &fname);

/** Check if a given path is a directory. */
bool is_directory(const std::string &path);

/** Check if a given path is a file. */
bool file_exists(const std::string &fname);


///////////////////////////////////////////////
// Binary File Handling
///////////////////////////////////////////////

/** Read a binary file and stores its content in a char vector. */
std::vector<char> read_binary_file(const std::string &fname);

/** Write data from a char vector into a binary file. */
uint64_t write_binary_file(const std::string &fname, const std::vector<char> & data);


///////////////////////////////////////////////
// XML File Handling
///////////////////////////////////////////////

/** Read a file in XML format. Results are stored as property tree. */
bool read_xml(boost::property_tree::ptree & xml_in, std::string xml_fname);

/** Write a file in XML format. Input must be a property tree. */
bool write_xml(boost::property_tree::ptree & xml_out, std::string xml_fname);

/** Convert a variable to an XML node. T must be a data type that can be cast to a string-like output format. */
template<typename T> boost::property_tree::ptree getXMLnode (T variable) {

	boost::property_tree::ptree node;

	try {
		node.put("", variable);
	} catch ( const std::exception &ex ) {
		std::cerr << "Failed to convert variable to XML output format." << std::endl;
	}

	return node;

}

/** Convert a vector to an XML node. T must be a data type that can be cast to a string-like output format. */
template<typename T> boost::property_tree::ptree getXMLnode_vector (std::vector<T> vector) {

  	boost::property_tree::ptree node;

  	for ( auto el = vector.begin(); el != vector.end(); ++el ) {
  		node.add_child("element", getXMLnode ( *el ));
  	}

  	return node;

}

//bool node_exist(boost::property_tree::ptree & ptree, std::string node) {
//	boost::optional< const boost::property_tree::ptree& > child = ptree.get_child_optional( node );
//	return ! (! child);
//}

#endif /* TOOLS_STATIC_H */
