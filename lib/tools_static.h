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
inline bool gp_compare (GenomePosType i,GenomePosType j) {
	if ( i.pos == j.pos )
		return i.gid < j.gid;
	return (i.pos < j.pos);
}

/**
 * Compare BamAlignmentRecords (from SeqAn library) by their position.
 * Primary field is the reference ID (rID), secondary field is the position (beginPos).
 * @param l First record.
 * @param r Second record.
 * @return true, if l has a lower positions than r.
 */
inline bool compare_records_by_pos(const seqan::BamAlignmentRecord & l, const seqan::BamAlignmentRecord & r) {
	if ( l.rID == r.rID )
		return l.beginPos < r.beginPos;
	return l.rID < r.rID;
}


/////////////////////////////////////
////////// Type convertion //////////
/////////////////////////////////////

/**
 * Convert a vector of a desired data type to a string of delimited values.
 * The string will be created by streaming objects of type T to a stringstream, thus << must be defined for T.
 * @param vector The vector of values of type T
 * @param delim Character which will be used as delimiter in the resulting string [default: ',']
 * @return String with delimited values of original type T.
 */
template<typename T> std::string join ( std::vector<T> vector, char delim = ',' ) {
	std::stringstream ss;
	for ( auto el=vector.begin(); el!=vector.end(); ++el ) {
		ss << (*el);
		if ( el != --vector.end() )
			ss << delim;
	}
	return ss.str();
}

/**
 * Split a std::string to a std::vector<std::string>.
 * This template is only enabled for arithmetic (numeric) data types and std::string.
 * @param target Reference to the target vector to store the split values.
 * @param s The input string.
 * @param delim_list A list of split delimiters.
 * @author Tobias Loka
 */
template <
	typename T,
	typename = typename std::enable_if<std::is_arithmetic<T>::value || std::is_same<T, std::string>::value, T>::type
> void split ( std::vector<T> &target, const std::string &s, std::string delim_list = split_chars ) {
	std::size_t prev = 0, pos;
	std::string next_value;
	while ( prev < std::string::npos ) {
		pos = s.find_first_of(delim_list, prev);
		if ( pos > prev ) {
			next_value = s.substr(prev, pos-prev);
			if ( ! next_value.empty() ) {
				try {
					target.push_back( boost::lexical_cast<T>( s.substr(prev, pos-prev)));
				}
				catch (boost::bad_lexical_cast & ex) {
					std::cerr << "WARN: Ignored invalid value during string splitting: " << s.substr(prev, pos-prev) << "(" << ex.what() << ")" << std::endl;
				}
			}
		}
		prev = pos == std::string::npos ? std::string::npos : pos+1;
	}
}

template<
	typename T,
	typename=typename std::enable_if<std::is_arithmetic<T>::value, T>::type
> Operations get_operation( T operation ) {
	switch ( operation ) {
	case MATCH:
		return MATCH;
		break;
	case NO_MATCH:
		return NO_MATCH;
		break;
	case DELETION:
		return DELETION;
		break;
	case INSERTION:
		return INSERTION;
		break;
	default:
		return MATCH;
	}
}


///////////////////////////////////
////////// File handling //////////
///////////////////////////////////

/**
 * Get total size of a file (in bytes)
 * @param fname Name of the file.
 * @return Size of the file.
 */
inline std::ifstream::pos_type get_filesize(const std::string &fname)
{
  std::ifstream in(fname, std::ios::binary | std::ios::ate);
  return in.tellg();
}

/**
 * Check if a given path is a directory.
 * @param path of interest.
 * @return true, if the given path is a directory.
 */
inline bool is_directory(const std::string &path) {
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

/**
 * Check if a given path is a file.
 * @param Path of interest.
 * @return true, if the given path is a file.
 */
inline bool file_exists(const std::string &fname) {
	return boost::filesystem::exists(fname);
}

/**
 * Convert a relative to an absolute path.
 * @param fname Input path.
 * @return Absolute path to fname.
 * @author Tobias Loka
 */
inline std::string absolute_path(std::string fname) {
	boost::filesystem::path input_path(fname);
	// TODO: Change to boost::filesystem::weakly_canonical for boost version >1_60
	return boost::filesystem::absolute(fname).string();
}

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

/**
 * Get the suffix for a specified file format.
 * @param format The desired file format.
 * @return File suffix as string (e.g., ".bam" for BAM format)
 */
std::string get_file_suffix ( OutputFormat format );


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
 * Get a number as a std::string with N digits (e.g., 6 -> "006" for N=3).
 * @param value The number to be formatted.
 * @param N The number of digits.
 * @param fill_char Character to fill with (e.g., 6 -> ",,1" for fill_char=',') [default='0']
 * @return The formatted number as std::string.
 */
template<typename T, typename=typename std::enable_if<std::is_arithmetic<T>::value, T>::type> std::string to_N_digits ( T value, CountType N, char fill_char = '0' ) {
	std::stringstream ss;
	ss << std::setw(N) << std::setfill(fill_char) << value;
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
inline std::vector<CountType> flowcell_layout_to_tile_numbers( CountType surfaceCount, CountType swathCount, CountType tileCount ) {
	std::vector<uint16_t> tiles_vec;
	for (uint16_t surf = 1; surf <= surfaceCount; surf++)
		for (uint16_t swath = 1; swath <= swathCount; swath++)
			for (uint16_t tile = 1; tile <= tileCount; tile++)
				tiles_vec.push_back(surf*1000 + swath*100 + tile);
	return tiles_vec;
}

/**
 * Convert the highest tile number to a vector of tiles under the assumption that all combinations of surface, swath and tile are used.
 * @param max_tile Maximum tile number.
 * @return Vector of all tiles that are included in the surface count, swath count and tile count defined by the maximum of the given input.
 */
inline std::vector<CountType> maxTile_to_tiles ( CountType max_tile ) {
	return flowcell_layout_to_tile_numbers( max_tile/1000, (max_tile%1000)/100, max_tile % 100 );
}

/**
 * Compute a MAPQ value from a given probability.
 * @param prob Probability to convert.
 * @param max_prob Maximal probability. This is necessary to have an upper boundary for the MAPQ. The default is 0.99993f, resulting in a MAPQ of 42.
 * @return The MAPQ value for the given probability.
 */
inline CountType prob2mapq(float prob, float max_prob = 0.99993f) {

	// Catch negative values and save computation time for value <0.1 that always have a MAPQ of 0
	if ( prob <= 0.1f)
		return 0;
	// Save computation time for several values up to 0.5
	else if ( prob <= 0.29f)
		return 1;
	else if ( prob <= 0.435f )
		return 2;
	else if ( prob <= 0.55f )
		return 3;

	// Otherwise calculate the correct value
	return ( float( (-10.0f) * std::log10( 1.0f - std::min(max_prob, prob ))) + 0.5f);
}

/**
 * Convert the output format to std::string.
 * @param format The output format to convert.
 * @return String of the output format.
 */
inline std::string to_string ( OutputFormat format ) {
	switch ( format ) {
	case BAM:
		return "BAM";
	case SAM:
		return "SAM";
	case CRAM:
		return "CRAM";
	default:
		return "BAM";
	}
}

/**
 * Convert the alignment mode to std::string.
 * @param mode The alignment mode to convert.
 * @param bestn The bestn value [default: 0]
 * @return String of the alignment mode.
 */
inline std::string to_string ( AlignmentMode mode, CountType bestn = 0 ) {
	switch ( mode ) {
	case ANYBEST:
		return "ANYBEST";
	case ALLBEST:
		return "ALLBEST";
	case ALL:
		return "ALL";
	case UNIQUE:
		return "UNIQUE";
	case UNKNOWN:
		return "UNKNOWN";
	case BESTN:
		return "BESTN" + std::to_string(bestn);
	default:
		return "ANYBEST";
	}
}

/**
 * Set the value of an immutable variable to a value without throwing an exception.
 * @param immutable The immutable to change the value.
 * @param value The new value for the immutable.
 * @return true, if the value could be set. False otherwise (e.g., if the immutable value was already set before).
 */
template<typename T>
bool set_immutable(Immutable<T> & immutable, T value) {
	  try {
		  immutable.set(value);
	  }
	  catch (immutable_error& e) {
		  std::cerr << "WARN: " << e.what() << std::endl;
		  return false ;
	  }
	  return true;
}

/**
 * Get the value of an immutable variable without throwing an exception.
 * If the value was not set, the default initialization value of type T will be returned.
 * @param immutable The immutable variable.
 * @return Value of the immutable variable. Default initialization value of T if the value of the immutable variable was not set yet.
 */
template<typename T>
T get_immutable(const Immutable<T> & immutable) {
	  try {
		  return immutable.get();
	  }
	  catch (immutable_error& e) {
		  std::cerr << "WARN: " << e.what() << std::endl;
		  return T();
	  }
}

/**
 * Convert a base call quality value to the respective char in PHRED syntax.
 * This function considers the settings of full quality or 2-bit quality in the globalAlignmentSettings.
 * @param bc_qual The base call quality as stored in HiLive.
 * @return PHRED char ( "!" - "I" )
 */
inline char to_phred_quality ( uint8_t bc_qual ) {
	char phred_score = '!';
	phred_score += bc_qual;
	return phred_score;
}

/**
 * Create a vector containing all lanes for Illumina HiSeq (1-8)
 * @return Vector containing all lanes for Illumina HiSeq.
 */
inline std::vector<uint16_t> all_lanes() {
  std::vector<uint16_t> ln;
  for (uint16_t l=0; l < 8; l++)
    ln.push_back(l+1);
  return ln;
}

/**
 * Create a vector containing all tiles for Illumina HiSeq (1101-2316)
 * @return Vector containing all tiles for Illumina HiSeq.
 */
inline std::vector<uint16_t> all_tiles() {
  return maxTile_to_tiles(2316);
}

#endif /* TOOLS_STATIC_H */

