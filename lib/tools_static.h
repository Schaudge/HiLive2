/**
 * This class provides functions that are independent of any other HiLive class.
 * Please do NOT add further includes to this file since this will lead to unwanted dependencies!
 */

#ifndef TOOLS_STATIC_H
#define TOOLS_STATIC_H

/* DONT ADD ANY INCLUDES */
#include "headers.h"
#include "definitions.h"
/* DONT ADD ANY INCLUDES */


/////////////////////////////////
////////// Comparators //////////
/////////////////////////////////

/**
 * Compare function to sort GenomePosType objects by position.
 * If position if equal, compare by gid.
 * @param i First position to compare
 * @param j Second position to compare
 * @return true, if first position is "smaller" than second position.
 */
bool gp_compare (GenomePosType i,GenomePosType j);


/////////////////////////////////////
////////// Type convertion //////////
/////////////////////////////////////

/**
 * Split a std::string to a std::vector<std::string>.
 * @param s Reference to the input string.
 * @param delim A split delimiter.
 * @param elems The target vector.
 * @author Tobias Loka
 */
void split(const std::string &s, char delim, std::vector<std::string> &elems);


///////////////////////////////////
////////// File handling //////////
///////////////////////////////////


/**
 * Get total size of a file (in bytes)
 * @param fname Name of the file.
 * @return Size of the file.
 */
std::ifstream::pos_type get_filesize(const std::string &fname);

/**
 * Check if a given path is a directory.
 * @param Path of interest.
 * @return true, if the given path is a directory.
 */
bool is_directory(const std::string &path);

/**
 * Check if a given path is a file.
 * @param Path of interest.
 * @return true, if the given path is a file.
 */
bool file_exists(const std::string &fname);

/**
 * Convert a relative to an absolute path.
 * @param fname Input path.
 * @return Absolute path to fname.
 * @author Tobias Loka
 * TODO: Not tested and used yet.
 */
std::string absolute_path(std::string fname);


/**
 * Read a binary file and stores its content in a char vector.
 * @param fname Path to the file.
 * @return All data from the file as char vector.
 */
std::vector<char> read_binary_file(const std::string &fname);

/**
 * Write data from a char vector into a binary file.
 * @param fname Path to the file.
 * @param data Data to be saved in the file.
 * @return Number of written bytes.
 */
uint64_t write_binary_file(const std::string &fname, const std::vector<char> & data);


////////////////////////////////////////////////
////////// Property trees / XML files //////////
////////////////////////////////////////////////

/**
 * Read a file in XML format. Results are stored as property tree.
 * @param xml_in Reference to the property tree to store the XML data.
 * @param xml_fname Name of the input file.
 * @return true on success
 * @author Tobias Loka
 */
bool read_xml(boost::property_tree::ptree & xml_in, std::string xml_fname);

/**
 * Write a property tree to an XML file.
 * @param xml_out Property tree that contains the data.
 * @param xml_fname Name of the output file.
 * @return true on success
 * @author Tobias Loka
 */
bool write_xml(boost::property_tree::ptree & xml_out, std::string xml_fname);

/**
 * Convert a variable of a non-vector type to a property tree.
 * @param variable The variable to convert.
 * @return The property tree for the input variable
 * @author Tobias Loka
 * TODO: check if the exception handling makes sense.
 */
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

/**
 * Convert a variable of a vector type to a property tree.
 * The subnodes have key "el".
 * @param vector The vector to convert.
 * @return The property tree for the input variable
 * @author Tobias Loka
 */
/** Convert a vector to an XML node. T must be a data type that can be cast to a string-like output format. */
template<typename T> boost::property_tree::ptree getXMLnode_vector (std::vector<T> vector) {

  	boost::property_tree::ptree node;

  	for ( auto el = vector.begin(); el != vector.end(); ++el ) {
  		node.add_child("el", getXMLnode ( *el ));
  	}

  	return node;

}


/////////////////////////////////
////////// Other stuff //////////
/////////////////////////////////

/**
 * Extract the number of reads from a BCL file.
 * @param bcl Path to the bcl file.
 * @return Number of reads in the bcl file.
 */
uint32_t num_reads_from_bcl(std::string bcl);

/**
<<<<<<< HEAD
 * Trim from start (in place).
 * @param s String to be trimmed.
 * @author Tobias Loka
 */
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
}

/**
 * Trim from end (in place).
 * @param s String to be trimmed.
 * @author Tobias Loka
 */
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

 /**
  * Trim from both ends (in place).
  * @param s String to be trimmed.
  * @author Tobias Loka
  */
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

/**
 * Convert the flowcell layout data to a plain vector of tile numbers.
 * @param surfaceCount The surface count.
 * @param swatchCount The swath count.
 * @param tileCount The tile count.
 * @return A vector of tile numbers.
 */
std::vector<CountType> flowcell_layout_to_tile_numbers( CountType surfaceCount, CountType swathCount, CountType tileCount );

CountType prob2mapq(float prob, float max_prob = 0.99993f);

#endif /* TOOLS_STATIC_H */

