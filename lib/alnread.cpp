#include "alnread.h"

seqan::String<seqan::CigarElement<> > Seed::returnSeqanCigarString() {

	bool extended_cigar = globalAlignmentSettings.get_extended_cigar();

	typedef seqan::String<seqan::CigarElement<> > TSeqanCigarString;
	TSeqanCigarString seqanCigarString;
	seqan::CigarElement<> cigarElem;

	for (CigarVector::const_iterator it = cigar_data.begin(); it != cigar_data.end(); ++it) {

		// Alignment begins with NO_MATCH => Softclip
		if (it == cigar_data.begin() && (*it).offset==NO_MATCH) {
			cigarElem.operation='S';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}

		// Mismatch => Alignment match
		if ((*it).offset==NO_MATCH) {
			cigarElem.operation= extended_cigar ? 'X' : 'M';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}

		// Deletion
		else if((*it).offset==DELETION) {
			cigarElem.operation='D';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}

		// Insertion
		else if((*it).offset==INSERTION) {
			cigarElem.operation='I';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}

		// Match
		else {
			cigarElem.operation= extended_cigar ? '=' : 'M';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}

	}

    // collapse Neighboring match regions
	for (unsigned k = 1; k<length(seqanCigarString); k++)
		if ((seqanCigarString[k-1].operation == 'M') && (seqanCigarString[k].operation == 'M')) {
			unsigned temp = seqanCigarString[k].count;
			erase(seqanCigarString, k);
			seqanCigarString[k-1].count += temp;
			k--;
		}
	
    return seqanCigarString;
}

void Seed::cout(){
	std::cout << "----- SEED START -----" << std::endl;
//	std::cout << "gid: " << this->gid << std::endl;
//	std::cout << "start_pos: " << this->start_pos << std::endl;
	std::cout << "num_errors: " << this->num_errors << std::endl;
	std::cout << "CIGAR: ";
	for ( auto el : this->cigar_data ) {
		std::cout << el.length;
		if ( el.offset == NO_MATCH ) {
			std::cout << "X ";
		}else {
			std::cout << "M(" << el.offset << ") ";
		}
	}
	std::cout << std::endl << "------ SEED END ------" << std::endl;
}

uint16_t Seed::serialize_size() {

  // Calculate total size
  uint16_t total_size = 0;

  // Size of FM-Index vertex descriptor
  total_size += sizeof(FMVertexDescriptor);

  // Number of errors
  total_size += sizeof(CountType);
  
  if (cigar_data.size() >= 256)
    throw std::overflow_error("CIGAR information contains more than 255 elements!");

  uint8_t cigar_len = cigar_data.size();

  // Number of CIGAR elements
  total_size += sizeof(uint8_t);

  // Length and offset for each CIGAR element
  total_size += cigar_len*(sizeof(CountType));

  return total_size;
}

std::vector<char> Seed::serialize() {

  // Total number of bytes after serialization
  uint16_t total_size = serialize_size();

  // Number of CIGAR elements
  uint8_t cigar_len = (uint8_t) cigar_data.size();

  // Char vector to store the data
  std::vector<char> data (total_size);
  char* d = data.data();

  // Write vertex descriptor
  memcpy(d,&vDesc,sizeof(FMVertexDescriptor));
  d += sizeof(FMVertexDescriptor);

  // Write number of errors
  memcpy(d,&num_errors,sizeof(CountType));
  d += sizeof(CountType);

  // Write number of CIGAR elements
  memcpy(d,&cigar_len,sizeof(uint8_t));
  d += sizeof(uint8_t);
  
  // Write the CIGAR elements themselves
  for (auto it = cigar_data.begin(); it != cigar_data.end(); ++it) {

	  // Serialize the CIGAR elements information
	  CountType serialized_value = it->length;
	  serialized_value = (serialized_value << 2);

	  // MATCH = 0; does not have to be handled
	  if (it->offset == NO_MATCH)
		  serialized_value |= CountType(3);
	  else if (it->offset == DELETION)
		  serialized_value |= CountType(2);
	  else if (it->offset == INSERTION)
		  serialized_value |= CountType(1);

	  memcpy(d,&(serialized_value),sizeof(CountType));
	  d += sizeof(CountType);

  }
  
  return data;
}

uint16_t Seed::deserialize(char* d) {
  
  // Total number of bytes read
  uint16_t bytes = 0; 

  // FM-index Vertex Descriptor
  memcpy(&vDesc,d+bytes,sizeof(FMVertexDescriptor));
  bytes += sizeof(FMVertexDescriptor);

  // Number of errors
  memcpy(&num_errors,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // Number of CIGAR elements
  uint8_t cigar_len = 0;
  memcpy(&cigar_len,d+bytes,sizeof(uint8_t));
  bytes += sizeof(uint8_t);

  // The CIGAR elements themselves
  cigar_data.clear();
  for (uint8_t i = 0; i < cigar_len; ++i) {

    CigarElement cig;

    // Deserialize CIGAR element
    CountType serialized_value;
    memcpy(&(serialized_value),d+bytes,sizeof(CountType));
    bytes += sizeof(CountType);

    DiffType offset_value = (serialized_value & CountType(3));

    // Store offset
    if (offset_value == 3)
    	cig.offset = NO_MATCH;
    else if (offset_value == 2)
    	cig.offset = DELETION;
    else if (offset_value == 1)
    	cig.offset = INSERTION;
    else
    	cig.offset = 0;

    // Store length
    cig.length = (serialized_value >> 2);

    // Add to CIGAR vector
    cigar_data.emplace_back(cig);
  }

  return bytes;  
}

uint64_t ReadAlignment::serialize_size() {

  // Total size
  uint64_t total_size = 0;
  
  // Flag
  total_size += 1;

  // Cycle number
  total_size += sizeof(CountType);

  // Last invalid cycle
  total_size += sizeof(CountType);

  // Sequence length
  total_size += sizeof(CountType);

  // The sequence information itself
  total_size += sequenceStoreVector.size()*(sizeof(uint8_t));

  // Barcode length
  total_size += sizeof(CountType);

  // The barcode sequence itself
  total_size += barcodeStoreVector.size()*(sizeof(uint8_t));

  // Number of seeds
  total_size += sizeof(uint32_t);

  // Size of the single seeds for each seed
  for (auto & s : seeds) {
    total_size += sizeof(uint16_t) + s->serialize_size();
  }
  
  return total_size;
}

std::vector<char> ReadAlignment::serialize() {

  // Total size
  uint64_t total_size = serialize_size();

  // Number of seeds
  uint32_t num_seeds = (uint32_t) seeds.size();

  // Char vector to store the data
  std::vector<char> data (total_size);
  char* d = data.data();
  
  // The flag
  memcpy(d,&flags,1);
  d++;

  // Sequence length
  memcpy(d,&sequenceLen,sizeof(CountType));
  d += sizeof(CountType);

  // The sequence itself
  for (auto it = sequenceStoreVector.begin(); it != sequenceStoreVector.end(); ++it) {
    memcpy(d,&(*it),sizeof(uint8_t));
    d += sizeof(uint8_t);
  }

  // Barcode length
  memcpy(d,&barcodeLen,sizeof(CountType));
  d += sizeof(CountType);

  // Barcode sequence itself
  for (auto it = barcodeStoreVector.begin(); it != barcodeStoreVector.end(); ++it) {
    memcpy(d,&(*it),sizeof(uint8_t));
    d += sizeof(uint8_t);
  }

  // Number of seeds
  memcpy(d,&num_seeds,sizeof(uint32_t));
  d += sizeof(uint32_t);
  
  // The seeds themselves
  for (auto it = seeds.begin(); it != seeds.end(); ++it) {

	  std::vector<char> seed_data = (*it)->serialize();
	  uint16_t seed_size = seed_data.size();

	  memcpy(d,&seed_size,sizeof(uint16_t));
	  d += sizeof(uint16_t);

	  memcpy(d,seed_data.data(),seed_size);
	  d += seed_size;
  }
  
  return data;
}

uint64_t ReadAlignment::deserialize(char* d) {

  // Total number of bytes
  uint64_t bytes = 0; 
  
  // The flag
  memcpy(&flags,d,1);
  bytes++;

  // Sequence length
  sequenceLen = 0;
  memcpy(&sequenceLen,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // The sequence itself
  unsigned seqVec_size = (unsigned) std::ceil((float) sequenceLen / 2.0);
  sequenceStoreVector.clear();
  sequenceStoreVector.reserve(seqVec_size);
  for (unsigned i = 0; i <seqVec_size; ++i) {
    uint8_t elem;
    memcpy(&(elem),d+bytes,sizeof(uint8_t));
    bytes += sizeof(uint8_t);
    sequenceStoreVector.push_back(elem);
  }

  // Barcode length
  barcodeLen = 0;
  memcpy(&barcodeLen,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // The barcode itself
  unsigned barVec_size = (unsigned) std::ceil((float) barcodeLen / 2.0);
  barcodeStoreVector.clear();
  barcodeStoreVector.reserve(barVec_size);
  for (unsigned i = 0; i <barVec_size; ++i) {
    uint8_t elem;
    memcpy(&(elem),d+bytes,sizeof(uint8_t));
    bytes += sizeof(uint8_t);
    barcodeStoreVector.push_back(elem);
  }

  // Number of seeds
  uint32_t num_seeds = 0;
  memcpy(&num_seeds,d+bytes,sizeof(uint32_t));
  bytes += sizeof(uint32_t);

  // The seeds themselves
  seeds.clear();
  for (uint32_t i = 0; i < num_seeds; ++i) {

	  uint16_t seed_size = 0;
	  memcpy(&seed_size,d+bytes,sizeof(uint16_t));
	  bytes += sizeof(uint16_t);
    
	  std::vector<char> seed_data (seed_size,0);
	  memcpy(seed_data.data(),d+bytes,seed_size);
	  bytes += seed_size;
    
	  USeed s (new Seed);
	  s->deserialize(seed_data.data());

	  // Must be sorted after extension anyways, so just push back
	  seeds.push_back(std::move(s));
  }

  return bytes;  
}

// TODO: handle Ns (use quality - maybe implement a unhash5() function or so considering quality)
std::string ReadAlignment::getSequenceString() {

    std::string seq = "";

    // Append 2 bases at a time (fitting into 1 byte)
    for (unsigned i = 0; i<sequenceStoreVector.size(); i++)
        seq.append(unhash(sequenceStoreVector[i], 2));

    // Delete overhang
    for (unsigned i = sequenceLen; i<2*sequenceStoreVector.size(); ++i)
        seq.pop_back();

    return seq;

}

std::string ReadAlignment::getBarcodeString() {


    std::string seq = "";

    // Append 2 bases at a time (fitting into 1 byte)
    for (unsigned i = 0; i<barcodeStoreVector.size(); i++)
        seq.append(unhash(barcodeStoreVector[i], 2));

    // Delete overhang
    for (unsigned i = barcodeLen; i<2*barcodeStoreVector.size(); ++i)
        seq.pop_back();

    return seq;

}

void ReadAlignment::appendNucleotideToSequenceStoreVector(char bc, bool appendToBarCode) {

	uint8_t nucl = bc & 3;
	uint8_t qual = bc >> 2;

	// Convert nucl and qual to four-bit value (first two bits describing the quality (0=N; 1=invalid; 2=valid); second two bits describing the nucleotide)
	uint8_t four_bit_repr = ( ((qual != 0) + (qual > globalAlignmentSettings.get_min_qual())) << 2 ) | nucl;

	CountType & len = appendToBarCode ? barcodeLen : sequenceLen;
	std::vector<uint8_t> & seqVector = appendToBarCode ? barcodeStoreVector : sequenceStoreVector;

    // check if all bits from sequenceStoreVector are used
    if (len % 2 == 0) { // yes, all used => new 8 Bit block needs to be created
        uint8_t newBlock = four_bit_repr << 4; // 'empty' bits are on the right side
        seqVector.push_back(newBlock);
        ++len;
    }

    else { // not all bits are used. There is enough space for the new basecall
        seqVector.back() = seqVector.back() | four_bit_repr;
        ++len;
    }
}

CountType ReadAlignment::getMaxNumErrors(USeed s) {

	// Determine length of the softclip
	CountType softclip_length = s->cigar_data.front().offset == NO_MATCH ? s->cigar_data.front().length : 0;

	// Determine the cycle without considering the softclip
	CountType seedCycle = cycle - softclip_length;

	// Determine the number of errors the seed startet with when created due to softclip length
	CountType softClipErrors = softclip_length == 0 ? 0 : ( ( softclip_length - 1 ) / globalAlignmentSettings.get_anchor_length() ) + 1;

	if (seedCycle <= globalAlignmentSettings.get_anchor_length() )
		return 0;

	return std::min(globalAlignmentSettings.get_min_errors(), CountType(( (seedCycle - globalAlignmentSettings.get_anchor_length() - 1) / globalAlignmentSettings.get_error_rate() ) + 1 + softClipErrors) );
}

void ReadAlignment::extendSeed(char base, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds){

	// The new base / nucleotide
	CountType tbr = twobit_repr(base);

	// Extend alignment match (can be match or NO_MATCH)
	getMatchSeeds(tbr, origin, allowedErrors, index, newSeeds);

	// Extend insertion
	getInsertionSeeds(tbr, origin, allowedErrors, index, newSeeds);

	// Extend deletion
	getDeletionSeeds(tbr, origin, allowedErrors, index, newSeeds);

}

void ReadAlignment::getMatchSeeds(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds) {

	// iterate through possible bases
	for (int b=0; b<4; b++) {

		FMTopDownIterator it(index->idx, origin->vDesc); 	// Create index iterator

		// Only handle when the path exist in the genome
		if (goDown(it,revtwobit_repr(b))) {

			// handle a MATCH
			if ( b == base_repr ) {

				// copy data from origin seed
				USeed s (new Seed);
				s->vDesc = it.vDesc;
				s->cigar_data = origin->cigar_data;
				s->num_errors = origin->num_errors;

				// Insert MATCH region if necessary
				if ( s->cigar_data.back().offset >= DELETION )
					s->cigar_data.emplace_back(0, 0); // TODO: add offset of previous match regions (if not deprecated)

				s->cigar_data.back().length += 1; // Increase match region length

				newSeeds.emplace_back(s); // Push the new seed to the vector.
			}

			// handle a NO_MATCH
			else {

				// Stop when maximal number of errors is reached
				if ( origin->num_errors >= allowedErrors )
					continue;

				// DON'T FINISH DELETION AND INSERTION REGIONS BY NO_MATCH!!!
				if ( origin->cigar_data.back().offset == DELETION )
					continue;

				if ( origin->cigar_data.back().offset == INSERTION )
					continue;

				// copy data from origin seed
				USeed s (new Seed);
				s->vDesc = it.vDesc;
				s->cigar_data = origin->cigar_data;
				s->num_errors = origin->num_errors + 1;

				// Insert MATCH region if necessary
				if ( s->cigar_data.back().offset != NO_MATCH )
					s->cigar_data.emplace_back(0, NO_MATCH);

				s->cigar_data.back().length += 1; // Increase match region length

				newSeeds.emplace_back(s); // Push the new seed to the vector.
			}
		}
	}
}

void ReadAlignment::getInsertionSeeds(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds) {

	// Don't handle insertions after deletions.
	if ( origin->cigar_data.back().offset == DELETION )
		return;

	// If all occurences of the seed match the current base, insertion is not required.
	FMTopDownIterator it(index->idx, origin->vDesc);
	uint64_t range = origin->vDesc.range.i2 - origin->vDesc.range.i1;
	if ( seqan::goDown(it, revtwobit_repr(base_repr)) ) {
		if ( it.vDesc.range.i2 - it.vDesc.range.i1 == range )
			return;
	}

	// handle all non-deletion regions
	if ( origin->cigar_data.back().offset != DELETION ) {

		// Extend insertion region as long as number of errors is allowed.
		if ( origin->num_errors < allowedErrors ) {

			// Don't go down the tree since an insertion means to stay at the same index position!

			// copy data from origin seed
			USeed s (new Seed);
			s->vDesc = origin->vDesc;
			s->cigar_data = origin->cigar_data;
			s->num_errors = origin->num_errors + 1; // Increase number of errors

			// Insert INSERTION region if necessary
			if ( s->cigar_data.back().offset != INSERTION )
				s->cigar_data.emplace_back(0, INSERTION);

			s->cigar_data.back().length += 1; // Increase insertion length

			newSeeds.emplace_back(s); // Push the new seed to the vector.
		}
	}
}

void ReadAlignment::getDeletionSeeds(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds) {

	// Don't handle insertion regions
	if ( origin->cigar_data.back().offset == INSERTION )
		return;

	// iterate through possible bases
	for (int b=0; b<4; b++) {

		// Don't start iteration with a match base (no deletion region required there!) or when having already all allowed errors.
		if ( b == base_repr || origin->num_errors >= allowedErrors )
			continue;

		// Create index iterator
		FMTopDownIterator it(index->idx, origin->vDesc);

		// Only consider paths existing in the index
		if ( goDown(it, revtwobit_repr(b))) {

			// copy data from origin seed
			USeed s (new Seed);
			s->vDesc = it.vDesc;
			s->cigar_data = origin->cigar_data;
			s->num_errors = origin->num_errors + 1; // Increase number of errors

			s->cigar_data.emplace_back(1, DELETION); // init deletion region
			recursive_goDown(base_repr, s, allowedErrors, index, newSeeds);
		}
	}
}

void ReadAlignment::recursive_goDown(CountType base_repr, USeed origin, CountType allowedErrors, KixRun* index, SeedVec & newSeeds) {

	if ( origin->num_errors > allowedErrors ) // too many errors -> no seeds
		return;

	// Only handle deletion regions
	if ( origin->cigar_data.back().offset != DELETION )
		return;

	// iterate through possible bases
		for (int b=0; b<4; b++) {

			// Create index iterator
			FMTopDownIterator it(index->idx, origin->vDesc);

			// Only consider paths existing in the index
			if ( goDown(it, revtwobit_repr(b))) {

				// Finish deletion region when a match occurs
				if ( b == base_repr ) {

					// copy data from origin seed
					USeed s (new Seed);
					s->vDesc = it.vDesc;
					s->cigar_data = origin->cigar_data;
					s->num_errors = origin->num_errors ;

					s->cigar_data.emplace_back(1, 0); // init MATCH region
					newSeeds.emplace_back(s);

				}

				// Continue deletion region otherwise
				else {

					// stop if maximum number of errors reached
					if ( origin->num_errors >= allowedErrors)
						continue;

					// copy data from origin seed
					USeed s (new Seed);
					s->vDesc = it.vDesc;
					s->cigar_data = origin->cigar_data;
					s->num_errors = origin->num_errors + 1; // Increase number of errors

					s->cigar_data.back().length += 1; // extend deletion region
					recursive_goDown(base_repr, s, allowedErrors, index, newSeeds);

				}
			}
		}
}

void ReadAlignment::createSeeds(KixRun* index, SeedVec & newSeeds) {

	if ( cycle < globalAlignmentSettings.get_anchor_length() )
		return;

	// use as exact mode? Becomes quite complex when using a high number of errors.
	uint32_t min_mismatches = cycle == globalAlignmentSettings.get_anchor_length() ? 0 : ( ( cycle - globalAlignmentSettings.get_anchor_length() - 1 ) / globalAlignmentSettings.get_anchor_length() ) + 1;

	// use as fast mode?
//	uint32_t min_mismatches = this->getMaxNumErrors(settings);

	if ( min_mismatches > globalAlignmentSettings.get_min_errors() )
		return;

	std::string anchor_seq = getSequenceString().substr( sequenceLen - globalAlignmentSettings.get_anchor_length(), globalAlignmentSettings.get_anchor_length());

	FMTopDownIterator it(index->idx);
	if ( seqan::goDown(it, anchor_seq) ) {
		USeed seed ( new Seed() );
		seed->num_errors = min_mismatches;
		if ( cycle != globalAlignmentSettings.get_anchor_length() )
			seed->cigar_data.emplace_back( ( cycle - globalAlignmentSettings.get_anchor_length()) , NO_MATCH ); // add softclip
		seed->cigar_data.emplace_back(globalAlignmentSettings.get_anchor_length(), 0);
		seed->vDesc = it.vDesc;
		newSeeds.emplace_back(seed);
	}

}

void ReadAlignment::extend_alignment(char bc, KixRun* index, bool testPrint) {

	if ( testPrint )
		std::cout << "ExtendAlignment 1" << std::endl;

    // cycle is not allowed to be > total_cycles
    assert( total_cycles >= cycle );

    appendNucleotideToSequenceStoreVector(bc);

    // do not update the alignments when reading the first kmer_span-1 cycles
    if ( cycle < globalAlignmentSettings.get_anchor_length() )
        return;

	SeedVec newSeeds;
	newSeeds.reserve(9*seeds.size()); // worst case estimation for only 1 error

	// extend existing seeds
    char base = revtwobit_repr(bc & 3);
    for ( auto seed = seeds.begin(); seed != seeds.end(); ++seed ) {
    	extendSeed(base, (*seed), getMaxNumErrors(*seed), index, newSeeds);
    }

    // create new seeds in defined intervals
//    if ( cycle % globalAlignmentSettings.get_anchor_length() == 0 ) {
    if ( cycle >= globalAlignmentSettings.get_anchor_length() ) {
    	createSeeds(index, newSeeds);
    }

	// sort the newSeeds vector
//    newSeeds.sort(PComp<USeed>);
    std::sort(newSeeds.begin(), newSeeds.end(), PComp<USeed>);

	// put only unique seeds into the seeds set.
	seeds.clear();

	if(newSeeds.size() == 0)
		return;

	if ( testPrint ) {
		std::cout << "---NEWSEEDS---" << std::endl;
		for(auto sd = newSeeds.begin(); sd != newSeeds.end(); ++sd) {
			std::cout << "CIGAR: " ;
			for (auto el = (*sd)->cigar_data.begin(); el != (*sd)->cigar_data.end(); ++el ) {
				std::cout << el->length;
				if(el->offset == NO_MATCH)
					std::cout << "X";
				else if (el->offset == DELETION)
					std::cout << "D";
				else if (el->offset == INSERTION)
					std::cout << "I";
				else
					std::cout << "M";
		    }
			std::cout << " NUM_ERRORS: " << (*sd)->num_errors;
			std::cout << " RANGE: " << (*sd)->vDesc.range;
			std::cout << std::endl;
		}
	}

	auto svit = newSeeds.begin();
	USeed lastUniqueSeed = (*svit);
	seeds.reserve(newSeeds.size());

	++svit;
	while( svit != newSeeds.end() ) {

		if ( lastUniqueSeed->vDesc != (*svit)->vDesc ) {
			seeds.emplace_back(lastUniqueSeed);
			lastUniqueSeed = (*svit);
		}
		++svit;
	}
	seeds.emplace_back(lastUniqueSeed); // pushback the last unique seed that was not pushed in the while loop.

	if ( testPrint ) {
		std::cout << "---SEEDS---" << std::endl;
		for(auto sd = seeds.begin(); sd != seeds.end(); ++sd) {
			std::cout << "CIGAR: " ;
			for (auto el = (*sd)->cigar_data.begin(); el != (*sd)->cigar_data.end(); ++el ) {
				std::cout << el->length;
				if(el->offset == NO_MATCH)
					std::cout << "X";
				else if (el->offset == DELETION)
					std::cout << "D";
				else if (el->offset == INSERTION)
					std::cout << "I";
				else
					std::cout << "M";
		    }
			std::cout << " NUM_ERRORS: " << (*sd)->num_errors;
			std::cout << " RANGE: " << (*sd)->vDesc.range;
			std::cout << std::endl;
		}
	}

	return;
}

CountType ReadAlignment::getBarcodeIndex() {

	// Get the barcodes of the read
	std::string read_bc = getBarcodeString();

	uint16_t fragment_errors = 0;
	uint16_t fragment_pos = 0;
	uint16_t fragment_num = 0;
	uint16_t matching_bc = NO_MATCH;

	// Iterate through all user-defined (multi-)barcodes
	// That's quite complicated since the read barcodes are consecutive and the user barcodes are divided in vectors. // TODO: change that?
	for ( uint16_t barcodeIndex = 0; barcodeIndex < globalAlignmentSettings.get_barcodeVector().size(); barcodeIndex++ ) {

		// reset values for the barcode
		fragment_errors = 0;
		fragment_pos = 0;
		fragment_num = 0;
		matching_bc = barcodeIndex;

		// for each base of the read barcode
		for ( uint16_t nucl = 0; nucl < read_bc.length(); nucl++ ) {

			// reset values for each barcode fragment
			if ( fragment_pos >= (globalAlignmentSettings.get_barcodeVector()[barcodeIndex])[fragment_num].length() ) {
				fragment_pos = 0;
				fragment_num += 1;
				fragment_errors = 0;
				assert( fragment_num < (globalAlignmentSettings.get_barcodeVector()[barcodeIndex]).size() );
			}

			// compare nucleotides and increase the number of fragment errors if not equal
			if ( read_bc.at(nucl) != (globalAlignmentSettings.get_barcodeVector()[barcodeIndex])[fragment_num].at(fragment_pos) ) {
				fragment_errors++;
			}

			// if too many errors in a fragment, break the loop for the barcode
			if ( fragment_errors > globalAlignmentSettings.get_barcode_errors()[fragment_num] ) {
				matching_bc = NO_MATCH;
				break;
			}

			fragment_pos += 1; // increment the fragment position

		}

		// if one barcode fulfilled the criteria, we can stop.
		if ( matching_bc != NO_MATCH )
			break;
	}

	return matching_bc;
}

void ReadAlignment::disable() {
  seeds.clear();
  flags = 0;
  sequenceLen=0;
  sequenceStoreVector.clear();
}

PositionType ReadAlignment::get_SAM_start_pos(KixRun* index, PositionPairType p, USeed & sd) {

	// Only valid if CIGAR string exist (should always be the case)
	if ( sd->cigar_data.size() == 0 )
		return std::numeric_limits<PositionType>::max();

	// Retrieve position from the index
	PositionType sam_pos = p.second;

	// REVERSE POSITIONS: position is already correct
	if ( index->isReverse(p.first) ) {
		return sam_pos;
	}

	// FORWARD POSITIONS: Compute start position
	else {
		sam_pos = index->getSequenceLength(p.first) - p.second;

		CountType cigLen = 0;
		for ( auto el = sd->cigar_data.begin(); el != sd->cigar_data.end(); ++el ) {
			if ( el->offset != INSERTION ) {
				cigLen += el->length;
			}
		}
		sam_pos -= cigLen;
	}

    return sam_pos;
}

void ReadAlignment::getPositions(KixRun* index, USeed sd, PositionPairListType & position_list) {

	seqan::Pair<unsigned> hitInterval = sd->vDesc.range;
	for (; hitInterval.i1 < hitInterval.i2; ++hitInterval.i1) {
	    	std::pair<GenomeIdType, PositionType> el ( seqan::getFibre(index->idx, seqan::FibreSA())[hitInterval.i1].i1, seqan::getFibre(index->idx, seqan::FibreSA())[hitInterval.i1].i2);
	    	position_list.push_back(el);
	}

}

void ReadAlignment::getSeeds_errorsorted(SeedVec & seeds_sorted) {
//	SeedVec seeds_sorted;
	seeds_sorted.insert(seeds_sorted.end(), seeds.begin(), seeds.end());
//	seeds_sorted.sort(seed_comparison_by_error);
	std::sort(seeds_sorted.begin(), seeds_sorted.end(), seed_comparison_by_error);
}

// Calculate the mapping quality for all alignments of the read based on the other alignments and the number of matching positions.
int16_t MAPQ(const SeedVec &sv){
  return sv.size();
}
