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


bool compare_records_by_pos(const seqan::BamAlignmentRecord & l, const seqan::BamAlignmentRecord & r);

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
 * @param path of interest.
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
 * Write a property tree to an ini file.
 * @param ini_out Property tree that contains the data.
 * @param ini_fname Name of the output file.
 * @return true on success
 * @author Tobias Loka
 */
bool write_ini(boost::property_tree::ptree & ini_out, std::string ini_fname);

/**
 * Add a value of type T to a given property tree.
 * @param ptree Reference to a given property tree.
 * @param value The value that will be added to the property tree.
 * @return True, if the value was successfully added to the property tree. False otherwise.
 */
template<typename T> bool putConfigNode (boost::property_tree::ptree & ptree, std::string key, T value) {

	bool success = true;

	try {
		ptree.put(key, value);
	} catch ( const std::exception &ex ) {
		std::cerr << "WARN: Failed to convert value to config output format." << std::endl;
		success = false;
	}

	return success;

}


/////////////////////////////////
////////// Other stuff //////////
/////////////////////////////////

/**
 * Convert a string of delimited values to a vector of a desired data type.
 * The string is split at every occurence of one of the delimiter defined in delim_list.
 * T will be created by streaming from a stringstream, thus >> must be defined for T.
 * @param values The values as string delimited by characters included in delim_list
 * @param delim_list List of delimiting characters [default: definitions.h --> split_chars]
 * @return Vector of the desired type for which >> must be defined (e.g. this holds for strings and numeric values)
 */
template<typename T> std::vector<T> to_vector ( std::string values, std::string delim_list=split_chars ) {

	std::vector<std::string> str_vector;
	std::vector<T> T_vector;

	boost::split(str_vector, values, [&delim_list](char c){return delim_list.find(c)!=delim_list.npos;});
	for ( auto & el : str_vector ) {
		if ( el.size() == 0 )
			continue;
		T next_el;
		std::stringstream ss(el);
		ss >> next_el;
		T_vector.push_back(next_el);
	}

	return T_vector;
}

/**
 * Convert a vector of a desired data type to a string of delimited values.
 * The string will be created by streaming objects of type T to a stringstream, thus << must be defined for T.
 * @param vector The vector of values of type T
 * @param delim Character which will be used as delimiter in the resulting string [default: ',']
 * @return String with delimited values of original type T.
 */
template<typename T> std::string to_string ( std::vector<T> vector, char delim = ',' ) {
	std::stringstream ss;
	for ( auto el=vector.begin(); el!=vector.end(); ++el ) {
		ss << (*el);
		if ( el != --vector.end() )
			ss << delim;
	}
	return ss.str();
}

/**
 * Get a number as a std::string with N digits (e.g., 6 -> "006" for N=3).
 * @param value The number to be formatted.
 * @param N The number of digits.
 * @return The formatted number as std::string.
 */
template<typename T, typename=typename std::enable_if<std::is_arithmetic<T>::value, T>::type> std::string to_N_digits ( T value, CountType N ) {
	std::stringstream ss;
	ss << std::setw(N) << std::setfill('0') << value;
	return ss.str();
}

/**
 * Check if a number contains a certain SAM flag.
 * @param value The value to check
 * @param flag The flag to check for
 * @return true, if value contains flag
 */
template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type> bool hasSAMFlag( T value, SAMFlag flag ) {
	return ( ( value & flag) == flag );
}

/**
 * Add a SAM flag to the total flag value.
 * @param value The total flag value
 * @param flag The flag to add
 * @return The new total flag value
 */
template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type> T addSAMFlag( T value, SAMFlag flag ) {
	return ( value | flag );
}

/**
 * Remove a SAM flag from the total flag value.
 * @param value The total flag value
 * @param flag The flag to remove
 * @return The new total flag value
 */
template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type> T removeSAMFlag( T value, SAMFlag flag ) {
	return ( addSAMFlag(value, flag) ^ flag );
}

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

