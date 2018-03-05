#include "tools.h"


//////////////////////////////////////
////////// Sequence Hashing //////////
//////////////////////////////////////

std::string unhash(HashIntoType myHash, unsigned hashLen)
{
	std::string seq = "";

	unsigned mask = 15;
	for (unsigned i = 1; i<pow(2,2*hashLen); i *= 4) {
//		kmer.push_back(revtwobit_repr(myHash & mask));
		seq.push_back(revfourbit_repr(myHash & mask));
		myHash = myHash >> 4;
	}
	std::reverse(seq.begin(), seq.end());
	return seq;
}



////////////////////////////////////////////
////////// File name construction //////////
////////////////////////////////////////////

std::string bcl_name(uint16_t ln, uint16_t tl, uint16_t cl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L" << to_N_digits(ln,3) << "/C" << cl << ".1/s_"<< ln <<"_" << tl << ".bcl";
  return path_stream.str();
}

std::string alignment_name(uint16_t ln, uint16_t tl, uint16_t cl, uint16_t mt){
  std::ostringstream path_stream;
  std::string base = globalAlignmentSettings.get_temp_dir() != "" ? globalAlignmentSettings.get_temp_dir() : globalAlignmentSettings.get_root();
  path_stream << base << "/L" << to_N_digits(ln,3) << "/s_"<< ln << "_" << tl << "." << mt << "."<< cl << ".align";
  return path_stream.str();
}

uint16_t getSeqCycle(uint16_t cycle, uint16_t seq_id) {
	uint16_t seq_cycle = cycle;
	for ( int i = 0; i < seq_id; i++ )
		seq_cycle += globalAlignmentSettings.getSeqById(i).length;
	return seq_cycle;
}

uint16_t getMateCycle( uint16_t mate_number, uint16_t seq_cycle ) {

	// Invalid mate
	if ( mate_number == 0 || mate_number > globalAlignmentSettings.get_mates() )
		return 0;

	// Iterate through all sequence elements (including barcodes)
	for ( CountType id = 0; id < globalAlignmentSettings.get_seqs().size(); id++ ) {

		// Current sequence element
		SequenceElement seq = globalAlignmentSettings.getSeqById(id);

		// Seq is mate of interest
		if ( seq.mate == mate_number )
			return ( seq.length > seq_cycle ? seq_cycle : seq.length );

		// Not enough cycles left to reach mate of interest
		else if ( seq.length >= seq_cycle )
			return 0;

		// Reduce number of cycles by the Seq length
		else
			seq_cycle -= seq.length;

	}

	// Should not be reached
	return 0;
}

std::string filter_name(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "/L" << to_N_digits(ln,3) << "/s_"<< ln << "_" << tl << ".filter";
  return path_stream.str();
}

std::string position_name(uint16_t ln, uint16_t tl) {
  std::ostringstream path_stream;
  path_stream << globalAlignmentSettings.get_root() << "../L" << to_N_digits(ln,3) << "/s_"<< ln << "_" << tl << ".clocs";
  return path_stream.str();
}

std::string get_settings_name() {
	std::ostringstream path_stream;
	std::string base = globalAlignmentSettings.get_temp_dir() != "" ? globalAlignmentSettings.get_temp_dir() : globalAlignmentSettings.get_root();
	path_stream << base << "/hilive_config.ini";
	return path_stream.str();
}

std::string get_out_log_name() {
	return ( globalAlignmentSettings.get_out_dir() + "/hilive_out.log" );
}


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

std::string getTileBamTempFileName(CountType ln, CountType tl, std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_output_format() == OutputFormat::BAM ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_temp_dir() << "/L" << to_N_digits(ln,3) << "/s_" << std::to_string(ln) << "_" << std::to_string(tl) << "." << std::to_string(cycle) << "." << barcode << ".temp" << file_suffix;
	return fname.str();
}

std::string getTileBamFileName(CountType ln, CountType tl, std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_output_format() == OutputFormat::BAM ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_temp_dir() << "/L" << to_N_digits(ln,3) << "/s_" << std::to_string(ln) << "_" << std::to_string(tl) << "." << std::to_string(cycle) << "." << barcode << file_suffix;
	return fname.str();
}

std::string getBamTempFileName(std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_output_format() == OutputFormat::BAM ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_out_dir() << "/hilive_out_" << "cycle" << std::to_string(cycle) << "_" << barcode << ".temp" << file_suffix;
	return fname.str();
}

std::string getBamFileName(std::string barcode, CountType cycle) {
	std::ostringstream fname;
	std::string file_suffix = globalAlignmentSettings.get_output_format() == OutputFormat::BAM ? ".bam" : ".sam";
	fname << globalAlignmentSettings.get_out_dir() << "/hilive_out_" << "cycle" << std::to_string(cycle) << "_" << barcode << file_suffix;
	return fname.str();
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


/////////////////////////////
////////// Scoring //////////
/////////////////////////////

uint16_t getMinSingleErrorPenalty() {

	// Mismatch: +1 mismatch, -1 match
	uint16_t mismatch_penalty = globalAlignmentSettings.get_mismatch_penalty() + globalAlignmentSettings.get_match_score();

	// Deletion: +1 deletion, maximum number of matches can still be reached
	uint16_t deletion_penalty = globalAlignmentSettings.get_deletion_opening_penalty() + globalAlignmentSettings.get_deletion_extension_penalty();

	// Insertion: +1 insertion, -1 match
	uint16_t insertion_penalty = globalAlignmentSettings.get_insertion_opening_penalty() + globalAlignmentSettings.get_insertion_extension_penalty() + globalAlignmentSettings.get_match_score();

	return std::min(mismatch_penalty, std::min(insertion_penalty, deletion_penalty));
}

uint16_t getMaxSingleErrorPenalty() {

	// Mismatch: +1 mismatch, -1 match
	uint16_t mismatch_penalty = globalAlignmentSettings.get_mismatch_penalty() + globalAlignmentSettings.get_match_score();

	// Deletion: +1 deletion, maximum number of matches can still be reached
	uint16_t deletion_penalty = globalAlignmentSettings.get_deletion_opening_penalty() + globalAlignmentSettings.get_deletion_extension_penalty();

	// Insertion: +1 insertion, -1 match
	uint16_t insertion_penalty = globalAlignmentSettings.get_insertion_opening_penalty() + globalAlignmentSettings.get_insertion_extension_penalty() + globalAlignmentSettings.get_match_score();

	return std::max(mismatch_penalty, std::max(insertion_penalty, deletion_penalty));
}

ScoreType getMaxPossibleScore( CountType cycles ) {
	return cycles * globalAlignmentSettings.get_match_score();
}

CountType getMinSoftclipPenalty( CountType softclip_length ) {
	return ceil( float(softclip_length) / globalAlignmentSettings.get_anchor_length() ) * getMinSingleErrorPenalty();
}

ScoreType getMinCycleScore( CountType cycle, CountType read_length ) {

	if ( cycle < globalAlignmentSettings.get_anchor_length() )
		return globalAlignmentSettings.get_min_as();

	ScoreType maxScore = getMaxPossibleScore(read_length);
	ScoreType minCycleScore = maxScore - ( ceil((cycle - globalAlignmentSettings.get_anchor_length()) / float(globalAlignmentSettings.get_error_rate())) * getMinSingleErrorPenalty() );
	return std::max(minCycleScore, globalAlignmentSettings.get_min_as());
}

int atomic_rename( const char *oldname, const char *newname ) {

	if ( !file_exists(oldname) )
		throw file_not_exist_error("Can't rename file: " + std::string(oldname) + ". File does not exist.");

	std::lock_guard<std::mutex> old_lock(fileLocks.at(std::string(oldname)));
	std::lock_guard<std::mutex> new_lock(fileLocks.at(std::string(newname)));

	return std::rename(oldname, newname);

}



/////////////////////////////////
////////// Other stuff //////////
/////////////////////////////////

char to_phred_quality ( uint8_t bc_qual ) {
	char phred_score = '!';
	phred_score += bc_qual;
	return phred_score;
}

bool isSeedingCycle(CountType cycle) {

	// Don't seed cycles smaller than the anchor length
	if ( cycle < globalAlignmentSettings.get_anchor_length() )
		return false;

	// Create seeds when reaching the anchor length for the first time
	if ( cycle == globalAlignmentSettings.get_anchor_length() )
		return true;

	// Create seeds every seeding_interval cycle after the first anchor
	if ( ( cycle - globalAlignmentSettings.get_anchor_length() ) % globalAlignmentSettings.get_seeding_interval() == 0 )
		return true;

	return false;

}
