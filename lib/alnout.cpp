#include "alnout.h"

//-------------------------------------------------------------------//
//------  Streamed SAM generation -----------------------------------//
//-------------------------------------------------------------------//

AlnOut::AlnOut(std::vector<CountType> lns, std::vector<CountType> tls, CountType cycl) : cycle(cycl) {

	// Fill list of specified barcodes
	std::vector<std::string> barcodes;
	for ( unsigned i = 0; i < globalAlignmentSettings.get_barcode_vector().size(); i++ ) {
		barcodes.push_back(globalAlignmentSettings.get_barcode_string(i));
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

	barcodes = globalAlignmentSettings.get_barcode_string_vector();

	// Init the bamIOContext (the same object can be used for all output streams)
	bfos.set_context(idx->getSeqNames(), idx->getSeqLengths());

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


	std::string in_fname = get_align_fname(ln, tl, cycle, mate);
	std::string out_fname = get_align_fname(ln, tl, cycle, mate) + ".sorted";

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

CountType AlnOut::openiAlnStream( CountType lane, CountType tile, CountType mateCycle, CountType mate, iAlnStream* istream) {

	// Check if mate is valid.
	if ( globalAlignmentSettings.get_seq_by_mate(mate) == NULLSEQ )
		return 1;

	// Check if mate cycle is >0
	if ( mateCycle == 0 )
		return 2;

	// Sort file if necessary
	if ( ! sort_tile( lane, tile, mate, mateCycle, globalAlignmentSettings.get_force_resort() ) )
		return 3;

	// Open sorted alignment file
	std::string alignment_fname = get_align_fname(lane, tile, mateCycle, mate) + ".sorted";

	// Check if sorted file exist
	if ( !file_exists(alignment_fname) ) {
		return 3;
	}

	// Open the stream
	istream->open(alignment_fname);

	return 0;
}

std::vector<iAlnStream*> AlnOut::openiAlnStreams( CountType lane, CountType tile, bool filter_exist, unsigned filter_size) {

	std::vector<iAlnStream*> alignmentFiles;
	unsigned numberOfAlignments = 0;

	for (unsigned mateIndex = 1; mateIndex <= mateCycles.size(); mateIndex++) {

		// Init the stream object
		iAlnStream* input = new iAlnStream( globalAlignmentSettings.get_block_size(), globalAlignmentSettings.get_compression_format());

		// Try to open the stream
		CountType ret = openiAlnStream(lane, tile, mateCycles[mateIndex-1], mateIndex, input);

		// On success
		if ( ret == 0 ) {

			// compare number of reads in alignment file with number of reads in filter file, if filter file exists
			if ( filter_exist && input->get_num_reads() != filter_size ) {
				std::stringstream msg;
				msg << "Unequal number of reads in filter file (" << filter_size << ") and alignment file (" << input->get_num_reads() << ")";
				throw std::length_error(msg.str());
			}

			// compare number of reads in alignment file with number of reads in previous alignment file
			if (mateIndex != 1 && input->get_num_reads() != numberOfAlignments) {
				throw std::length_error("Unequal number of reads (between mates)");
			}

			alignmentFiles.push_back(input);
			numberOfAlignments = input->get_num_reads();
			continue;
		}

		// Throw an exception if a mate is not valid (should not happen)
		else if ( ret == 1 ) {
			throw std::runtime_error("Mate number is not valid: " + mateIndex);
		}

		// Mate not sequenced yet (mateCycle==0)
		else if ( ret == 2 ) {
			continue;
		}

		// Sorted file is not available
		else if ( ret == 3 ) {
			throw std::runtime_error("Sorted align file not found.");
		}

		// Undefined return value
		else {
			throw std::logic_error("Undefined return value of function AlnOut::openiAlnStream: " + ret);
		}

	}

	return alignmentFiles;

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

void AlnOut::__write_tile_to_bam__ ( Task t ) {

	//TODO: Write real paired-end output (?)

	CountType lane = t.lane;
	CountType tile = t.tile;

	////////////////////////////////////////////////////
	//  Main loop //////////////////////////////////////
	////////////////////////////////////////////////////

	// set the filter file
	std::string filter_fname = get_filter_fname(lane, tile);
	FilterParser filters;
	if (file_exists(filter_fname)) {
		filters.open(filter_fname);
	}

	// Open the input streams for the sorted alignment files.
	std::vector<iAlnStream*> alignmentFiles = openiAlnStreams( lane, tile, file_exists(filter_fname), filters.size());
	unsigned numberOfAlignments = alignmentFiles[0]->get_num_reads();

	// for all reads in a tile
	/////////////////////////////////////////////////////////////////////////////
	for (uint64_t i = 0; i < numberOfAlignments; i++) {

		std::vector<ReadAlignment*> mateAlignments;
		for (auto e:alignmentFiles) {
			mateAlignments.push_back(e->get_alignment());
		}

		std::vector<std::vector<seqan::BamAlignmentRecord>> mateRecords(mateAlignments.size(), std::vector<seqan::BamAlignmentRecord>(0));

		// if the filter file is available and the filter flag is 0 then skip
		if (filters.size() != 0 && filters.next() == false)
			continue;

		// compute barcode sequence as it should be written to BC tag
		std::string barcode = globalAlignmentSettings.format_barcode(mateAlignments[0]->getBarcodeString());

		// Barcode index for the read
		CountType barcodeIndex = mateAlignments[0]->getBarcodeIndex();

		// If read has undetermined barcode and keep_all_barcodes is not set, skip this read
		if ( barcodeIndex == UNDETERMINED && !globalAlignmentSettings.get_keep_all_barcodes() )
			continue;
		else if ( barcodeIndex == UNDETERMINED )
			barcodeIndex = barcodes.size(); // this is the index for the "undetermined" output stream

		// setup QNAME
		// Read name format <instrument‐name>:<run ID>:<flowcell ID>:<lane‐number>:<tile‐number>:<x‐pos>:<y‐pos>
		// readname << "<instrument>:<run-ID>:<flowcell-ID>:" << ln << ":" << tl << ":<xpos>:<ypos>:" << i;
		//TODO: Get coordinates from clocs file, run-ID from runInfo.xml, flowcell ID from runInfo.xml, instrument from runInfo.xml
		std::stringstream readname;
		readname << "lane." << lane << "|tile." << tile << "|read." << i;

		// Track equivalent alignments for the same mate (similar positions from different seeds -> only keep the best one)
		// TODO: implement equivalent alignment window as user parameter
		PositionType equivalentAlignmentWindow = 10;

		// TODO: Improve the prevention of reporting similar alignments.
//		std::set<PositionType> alignmentPositions;
		std::set<GenomePosType> alignmentPositions;

		// for all mates
		/////////////////////////////////////////////////////////////////////////////
		for (unsigned mateAlignmentIndex=0; mateAlignmentIndex < mateAlignments.size(); ++mateAlignmentIndex) {

			// Decrease the ReadAlignment cycle since it is automatically increased when loading the file TODO: this should be changed.
			mateAlignments[mateAlignmentIndex]->cycle -= 1;

			// Init record with information about the read
			seqan::BamAlignmentRecord mate_record;
			mate_record.qName = readname.str();
			mate_record.flag = 0;

			// Set correct segment (paired) flags for the record
			if ( mateAlignments.size() > 1) { // if there are at least two mates already sequenced

				mate_record.flag = addSAMFlag(mate_record.flag, SAMFlag::MULT_SEG);
				if (mateAlignmentIndex == 0) {
					mate_record.flag |= addSAMFlag(mate_record.flag, SAMFlag::FIRST_SEG);
				} else if (mateAlignmentIndex == mateAlignments.size()-1) {
					mate_record.flag = addSAMFlag(mate_record.flag, SAMFlag::LAST_SEG);
				} else {
					mate_record.flag = addSAMFlag(mate_record.flag, SAMFlag::FIRST_SEG);
					mate_record.flag = addSAMFlag(mate_record.flag, SAMFlag::LAST_SEG);
				}

			}

			// Add read specific information to the BAM record
			mateAlignments[mateAlignmentIndex]->addReadInfoToRecord(mate_record);

			// Alignment disabled or no seeds
			if ( mateAlignments[mateAlignmentIndex]->is_disabled() || mateAlignments[mateAlignmentIndex]->seeds.size() == 0 ) {

				// Don't report disabled reads if their sequences are not kept.
				if ( mateAlignments[mateAlignmentIndex]->is_disabled() && !globalAlignmentSettings.get_keep_all_sequences() )
					continue;

				// Report unmapped reads if activated
				if ( globalAlignmentSettings.get_report_unmapped() ) {
					mate_record.flag = addSAMFlag(mate_record.flag, SAMFlag::SEG_UNMAPPED);
					mateRecords[mateAlignmentIndex].push_back(mate_record);
				}
				continue;
			}

			// Variables for output modes
			ScoreType first_seed_score = 0;
			ScoreType last_seed_score = 0;
			CountType num_diff_scores = 0;

			// Number of printed alignments for the current mate.
			unsigned printedMateAlignments = 0;

			// Unique mode interruption
			// TODO: Think about reporting of unmapped and non-unique reads if report-unmapped is activated
			// TODO: this filtering approach is problematic if similar seeds exist ...
			if ( globalAlignmentSettings.is_mode(UNIQUE) && (
					mateAlignments[mateAlignmentIndex]->seeds.size() > 1 ||
					(mateAlignments[mateAlignmentIndex]->seeds.size() == 1 && mateAlignments[mateAlignmentIndex]->seeds.front()->getNumPositions() > 1))) {
				continue;
			}

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
				if( globalAlignmentSettings.is_mode(ALLBEST) && first_seed_score > curr_seed_score ) {
					goto nextmate;
				}


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
				std::vector<GenomePosType> pos_list;

				if ( globalAlignmentSettings.is_mode(ANYBEST) ) {
					pos_list = (*it)->getPositions(0, 1); // retrieve only one position from the index
				} else
					pos_list = (*it)->getPositions(); // retrieve all positions from the index

				// handle all positions
				for ( auto p = pos_list.begin(); p != pos_list.end(); ++p ) {

					// Stop in any best mode when first alignment was already written
					if( globalAlignmentSettings.is_mode(ANYBEST) && printedMateAlignments >= 1 ) {
						goto nextmate;
					}

					// Stop in bestn mode if first n alignments were already written
					if ( globalAlignmentSettings.is_mode(BESTN) && printedMateAlignments >= globalAlignmentSettings.get_best_n() ) {
						goto nextmate;
					}

					seqan::BamAlignmentRecord record = mate_record;

					record.rID = CountType(p->gid / 2);

					record.beginPos = mateAlignments[mateAlignmentIndex]->get_SAM_start_pos(*p, *it);

					record.mapQ = *mapqs_it;

					// skip invalid positions
					if (record.beginPos < 0 || PositionType(record.beginPos) == std::numeric_limits<PositionType>::max()) {
						continue;
					}

					// skip positions that were already written (equivalent alignments). This can be done because the best alignment for this position is written first.
					if ( alignmentPositions.find(GenomePosType(p->gid, record.beginPos - ( record.beginPos % equivalentAlignmentWindow ) )) != alignmentPositions.end() ||
							alignmentPositions.find(GenomePosType(p->gid, record.beginPos +  (equivalentAlignmentWindow - ( record.beginPos % equivalentAlignmentWindow ) ) ) ) != alignmentPositions.end()) {
						continue;
					}

					record.cigar = cigar;
					if ( idx->isReverse(p->gid) )
						seqan::reverse(record.cigar);

					if ( printedMateAlignments > 0 ) { // if current seed is secondary alignment
						record.flag = addSAMFlag(record.flag, SAMFlag::SEC_ALIGNMENT);
						seqan::clear(record.seq);
						seqan::clear(record.qual);
					}

					if ( idx->isReverse(p->gid) ) { // if read matched reverse complementary
						seqan::reverseComplement(record.seq);
						seqan::reverse(record.qual);
						record.flag = addSAMFlag(record.flag, SAMFlag::SEQ_RC);
					}

					// Dictionary for additional SAM tags
					seqan::BamTagsDict dict;

					// Alignment Score
					seqan::appendTagValue(dict, "AS", curr_seed_score);

					// Barcode sequence
					if (barcode!="")
						seqan::appendTagValue(dict, "BC", barcode);

					// Number of mismatches
					seqan::appendTagValue(dict, "NM", nm);

					// MD:Z string
					std::string mdz = (*it)->getMDZString();
					if ( idx->isReverse(p->gid))
						mdz = reverse_mdz(mdz);
					seqan::appendTagValue(dict, "MD", mdz);

					record.tags = seqan::host(dict);

					// fill records list
					mateRecords[mateAlignmentIndex].push_back(record);

					// set variables for mode selection
					if ( last_seed_score != curr_seed_score || num_diff_scores == 0 )
						++num_diff_scores;
					last_seed_score = curr_seed_score;

					++printedMateAlignments;
					alignmentPositions.insert(GenomePosType(p->gid, record.beginPos - ( record.beginPos % equivalentAlignmentWindow )));

				}
			}
			nextmate: {
				// Report unmapped reads if activated
				if ( mateRecords[mateAlignmentIndex].size()==0 && globalAlignmentSettings.get_report_unmapped() ) {
					mate_record.flag = addSAMFlag(mate_record.flag, SAMFlag::SEG_UNMAPPED);
					mateRecords[mateAlignmentIndex].push_back(mate_record);
				}
			};
		}

		// Set flags related to the next mate.
		setMateSAMFlags(mateRecords);

		// Write all records as a group to keep suboptimal alignments and paired reads together.
		bfos[barcodeIndex].writeRecords(mateRecords);

		for (auto e:mateAlignments)
			delete e;

	}
	for (auto e:alignmentFiles)
		delete e;

	return;
}

void AlnOut::setMateSAMFlags( std::vector<std::vector<seqan::BamAlignmentRecord>> & mateRecords ) {
	for ( CountType i=0; i<mateRecords.size(); i++ ) {

		// Stop if not paired
		if ( mateRecords.size() <= 1 )
			break;

		// Set related mate index
		CountType next = i==mateRecords.size()-1 ? 0 : i+1;

		for ( auto & record_pointer : mateRecords[i] ) {

			if ( mateRecords[next].size() > 0 ) {

				CountType mateFlag = mateRecords[next][0].flag;

				// set next mate rID
				record_pointer.rNextId = mateRecords[next][0].rID;

				// set next mate pos
				record_pointer.pNext = mateRecords[next][0].beginPos;

				// other mate entry is reverse complemented
				if ( hasSAMFlag( mateFlag, SAMFlag::SEQ_RC ) ) {
					record_pointer.flag = addSAMFlag(record_pointer.flag, SAMFlag::NEXT_SEQ_RC);
				}

				// other mate entry is flagged as unmapped
				if ( hasSAMFlag( mateFlag, SAMFlag::SEG_UNMAPPED ) ) {
					record_pointer.flag = addSAMFlag(record_pointer.flag, SAMFlag::NEXT_SEG_UNMAPPED);

					record_pointer.pNext = record_pointer.beginPos;
					record_pointer.rNextId = record_pointer.rID;

				// the current read is flagged as unmapped but the mate is not
				} else if ( hasSAMFlag(record_pointer.flag, SAMFlag::SEG_UNMAPPED) ) {
					record_pointer.beginPos = mateRecords[next][0].beginPos;
					record_pointer.rID = mateRecords[next][0].rID;
				}


			// No record for the mate available (this implies, that it is unmapped)
			} else {
				record_pointer.flag = addSAMFlag(record_pointer.flag, SAMFlag::NEXT_SEG_UNMAPPED);
			}
		}

	}
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
