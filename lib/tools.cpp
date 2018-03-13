#include "tools.h"

////////////////////////////////////
////////// SAM/BAM output //////////
////////////////////////////////////

seqan::BamHeader getBamHeader() {
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

	return header;
}

std::string reverse_mdz(std::string mdz) {

	std::string reverse = "";
	std::string next_match = "";

	for ( auto string_it = mdz.rbegin(); string_it!=mdz.rend(); ++string_it) {
		if ( seq_chars.find(*string_it) != std::string::npos ) {
			std::reverse(next_match.begin(), next_match.end());
			reverse += next_match;
			next_match = "";
			reverse += comp(*string_it);
		}
		else {
			next_match += *string_it;
		}
	}
	std::reverse(next_match.begin(), next_match.end());
	reverse += next_match;
	return reverse;
}


/////////////////////////////////
////////// Other stuff //////////
/////////////////////////////////

int atomic_rename( const char *oldname, const char *newname ) {

	if ( !file_exists(oldname) )
		throw file_not_exist_error("Can't rename file: " + std::string(oldname) + ". File does not exist.");

	std::lock_guard<std::mutex> old_lock(fileLocks.at(std::string(oldname)));
	std::lock_guard<std::mutex> new_lock(fileLocks.at(std::string(newname)));

	return std::rename(oldname, newname);

}
