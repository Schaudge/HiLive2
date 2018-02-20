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
	std::cout << "Max. alignment score: " << this->max_as << std::endl;
	std::cout << "CIGAR: ";
	for ( auto el : this->cigar_data ) {
		std::cout << el.length;
		if ( el.offset == NO_MATCH ) {
			std::cout << "X ";
		}
		else if ( el.offset == INSERTION ) {
			std::cout << "I ";
		}
		else if ( el.offset == DELETION ) {
			std::cout << "D ";
		}
		else {
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

  // Maximal AS
  total_size += sizeof(ScoreType);
  
  if (cigar_data.size() >= 256)
    throw std::overflow_error("CIGAR information contains more than 255 elements!");

  uint8_t cigar_len = cigar_data.size();

  // Number of CIGAR elements
  total_size += sizeof(uint8_t);

  // Length and offset for each CIGAR element
  total_size += cigar_len*(sizeof(CountType));

  // Size of mdz_nucleotides vector
  total_size += sizeof(uint8_t);

  // MD:Z nucleotides
  total_size += mdz_nucleotides.size() * sizeof(uint8_t);

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
  memcpy(d,&max_as,sizeof(ScoreType));
  d += sizeof(ScoreType);

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
  
  // Write size of mdz_nucleotides vector
  memcpy(d,&mdz_length,sizeof(uint8_t));
  d += sizeof(uint8_t);

  // MD:Z nucleotides
  for ( auto it = mdz_nucleotides.begin(); it != mdz_nucleotides.end(); ++it) {
	  uint8_t next_byte = *it;

	  memcpy(d,&next_byte,sizeof(uint8_t));
	  d += sizeof(uint8_t);
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
  memcpy(&max_as,d+bytes,sizeof(ScoreType));
  bytes += sizeof(ScoreType);

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

  // Read size of mdz_nucleotides vector
  memcpy(&mdz_length, d+bytes, sizeof(uint8_t));
  bytes += sizeof(uint8_t);

  // MD:Z nucleotides
  for ( CountType i = 0; i < (mdz_length+3)/4; i++ ) {
	  uint8_t next_byte;

	  memcpy(&next_byte, d+bytes, sizeof(uint8_t));
	  bytes += sizeof(uint8_t);

	  mdz_nucleotides.push_back(next_byte);

  }

  return bytes;  
}

ScoreType Seed::get_as() {

	ScoreType as = 0;

	auto cigar_it = cigar_data.begin();


	for ( ;cigar_it != cigar_data.end(); ++cigar_it ) {

		if ( cigar_it->length == 0 )
			continue;

		switch ( cigar_it->offset ) {

		case NO_MATCH:

			// Softclip
			if ( cigar_it == cigar_data.begin() || std::next(cigar_it) == cigar_data.end() )
				as -= ( globalAlignmentSettings.get_softclip_opening_penalty() + ( ( cigar_it->length ) * globalAlignmentSettings.get_softclip_extension_penalty() ) );

			// Regular mismatch
			else
				as -= cigar_it->length * globalAlignmentSettings.get_mismatch_penalty();

			break;

		case INSERTION:

			// Opening
			as -= globalAlignmentSettings.get_insertion_opening_penalty();

			// Extension
			as -= ( ( cigar_it->length ) * globalAlignmentSettings.get_insertion_extension_penalty() );

			break;

		case DELETION:

			// Opening
			as -= globalAlignmentSettings.get_deletion_opening_penalty();

			// Extension
			as -= ( ( cigar_it->length ) * globalAlignmentSettings.get_deletion_extension_penalty() );

			break;

		default:

			as += ( cigar_it->length * globalAlignmentSettings.get_match_score() );

			break;
		}

	}

	return as;
}

CountType Seed::get_nm() {

	CountType nm = 0;

	// Don't count mismatches at front or end of the CIGAR string
	for ( auto el = ++(cigar_data.begin()); el != cigar_data.end(); ++el )
		nm += el->offset >= DELETION ? el->length : 0;

	return nm;
}

CountType Seed::get_softclip_length() {
	if ( cigar_data.size() == 0 )
		return 0;

	CountType sc_length = cigar_data.front().offset == NO_MATCH ? cigar_data.front().length : 0;
	return sc_length;
}
void Seed::add_mdz_nucleotide(char nucl) {

	uint8_t n = twobit_repr(nucl) << (6 - ( 2 * (mdz_length % 4) ) );
	if ( mdz_length % 4 == 0 )
		mdz_nucleotides.push_back(n);
	else {
		mdz_nucleotides.back() = mdz_nucleotides.back() | n;
	}
	mdz_length += 1;

}

std::string Seed::getMDZString() {

	std::string mdz_string = "";
	uint8_t nucleotide_pos = 0;
	auto mdz_it = mdz_nucleotides.begin();
	CountType match_counter = 0;

	for ( auto el = cigar_data.begin(); el != cigar_data.end(); ++el ) {

		// Softclip --> not included in MD:Z
		if ( el == cigar_data.begin() && el->offset == NO_MATCH )
			continue;

		switch ( el->offset ) {

		// DELETION or NO_MATCH: Add previous match length and nucleotides for the current region.
		case DELETION:
		case NO_MATCH:

			if ( match_counter != 0 )
				mdz_string += std::to_string(match_counter);

			for ( CountType i=0; i<el->length; i++ ) {
				mdz_string += revtwobit_repr( ( (*mdz_it) >> ( 6 - (2*nucleotide_pos) ) ) & 3);
				if ( nucleotide_pos >= 3 ) {
					nucleotide_pos = 0;
					++mdz_it;
				}
				else {
					nucleotide_pos += 1;
				}
			}

			match_counter = 0;
			break;

		// INSERTION: Do nothing.
		case INSERTION:
			break;

		// Match: Count the region length.
		default:
			match_counter += el->length;
			break;
		}

	}

	// Add last match sequence if exist
	if ( match_counter != 0 )
		mdz_string += std::to_string(match_counter);

	return mdz_string;
}

std::vector<GenomePosType> Seed::getPositions( CountType firstPosition, CountType lastPosition ) {

	std::vector<GenomePosType> position_list;

	seqan::Pair<unsigned> hitInterval = vDesc.range;

	// Return empty list if there is nothing to return.
	if ( lastPosition <= firstPosition || hitInterval.i2 - hitInterval.i1 <= firstPosition )
		return position_list;

	// set span values
	unsigned min_i1 = hitInterval.i1 + firstPosition;
	unsigned max_i2 = lastPosition == MAX_NUM_POSITIONS ? hitInterval.i2 : std::min(hitInterval.i2, hitInterval.i1 + lastPosition);

	for (; min_i1 < max_i2; ++min_i1) {
		GenomePosType el ( seqan::getFibre(idx->idx, seqan::FibreSA())[min_i1].i1, seqan::getFibre(idx->idx, seqan::FibreSA())[min_i1].i2);
		position_list.push_back(el);
	}

	return position_list;

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
  unsigned seqVec_size = sequenceLen;
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
  unsigned barVec_size = barcodeLen;
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

std::string ReadAlignment::getSequenceString() {

	std::string seq = "";
	uint8_t two_bit_mask = 3;

	// iterate through all sequence bytes
	for (unsigned i = 0; i<sequenceLen; i++ ) {

		// Next basecall (6 bits quality; 2 bits nucleotide)
		uint8_t next = sequenceStoreVector[i];

		if ( next < 4 ) { // qual == 0 --> N-call
			seq.append("N");
		} else {		  // qual > 0  --> write nucleotide
			seq += revtwobit_repr(next & two_bit_mask);
		}
	}

	// return barcode sequence
	return seq;

}

std::string ReadAlignment::getBarcodeString() {

	std::string seq = "";
	uint8_t two_bit_mask = 3;

	// iterate through all sequence bytes
	for (unsigned i = 0; i<barcodeLen; i++ ) {

		// Next basecall (6 bits quality; 2 bits nucleotide)
		uint8_t next = barcodeStoreVector[i];

		if ( next < 4 ) { // qual == 0 --> N-call
			seq.append("N");
		} else {		  // qual > 0  --> write nucleotide
			seq += revtwobit_repr(next & two_bit_mask);
		}
	}

    // return barcode sequence
    return seq;
}

std::string ReadAlignment::getQualityString() {

	std::string qual = "";

	// iterate through all sequence bytes
	for (unsigned i = 0; i<sequenceLen; i++ ) {

		// Next Quality (shift by 2 bit nucleotide information)
		uint8_t next_qual = sequenceStoreVector[i] >> 2;

		qual += (to_phred_quality(next_qual));
	}

	// return PHRED quality sequence
	return qual;
}


void ReadAlignment::appendNucleotideToSequenceStoreVector(char bc, bool appendToBarCode) {

	// Store byte
	CountType & len = appendToBarCode ? barcodeLen : sequenceLen;
	std::vector<uint8_t> & seqVector = appendToBarCode ? barcodeStoreVector : sequenceStoreVector;
	seqVector.push_back(bc);
	++len;
	return;

}

void ReadAlignment::extendSeed(char base, USeed origin, SeedVec & newSeeds){

	// The new base / nucleotide
	CountType tbr = twobit_repr(base);

	// Extend alignment match (can be match or NO_MATCH)
	getMatchSeeds(tbr, origin, newSeeds);

	// Extend insertion (only if not in the last cycle)
	if ( cycle != total_cycles )
		getInsertionSeeds(tbr, origin, newSeeds);

	// Extend deletion (may occur in the last cycle since it must always be closed by a match when created)
	getDeletionSeeds(tbr, origin, newSeeds);

}

void ReadAlignment::getMatchSeeds(CountType base_repr, USeed origin, SeedVec & newSeeds) {

	// iterate through possible bases
	for (int b=0; b<4; b++) {

		FMTopDownIterator it(idx->idx, origin->vDesc); 	// Create index iterator

		// Only handle when the path exist in the genome
		if (goDown(it,seqan::DnaString(revtwobit_repr(b)))) {

			// handle a MATCH
			if ( b == base_repr ) {

				// copy data from origin seed
				USeed s (new Seed);
				s->vDesc = it.vDesc;
				s->cigar_data = origin->cigar_data;
				s->mdz_nucleotides = origin->mdz_nucleotides;
				s->max_as = origin->max_as;
				s->mdz_length = origin->mdz_length;

				// Insert MATCH region if necessary
				if ( s->cigar_data.back().offset >= DELETION )
					s->cigar_data.emplace_back(0, 0); // TODO: add offset of previous match regions (if not deprecated)

				s->cigar_data.back().length += 1; // Increase match region length

				newSeeds.emplace_back(s); // Push the new seed to the vector.
			}

			// handle a NO_MATCH
			else {

				// Compute new max_as
				ScoreType new_max_as = origin->max_as - globalAlignmentSettings.get_mismatch_penalty() - globalAlignmentSettings.get_match_score();

				// If max_as is no longer valid
				if ( new_max_as < getMinCycleScore(cycle, total_cycles) )
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
				s->mdz_nucleotides = origin->mdz_nucleotides;
				s->max_as = new_max_as;
				s->mdz_length = origin->mdz_length;

				// Insert MATCH region if necessary
				if ( s->cigar_data.back().offset != NO_MATCH )
					s->cigar_data.emplace_back(0, NO_MATCH);

				s->cigar_data.back().length += 1; // Increase match region length

				// Add nucleotide for MDZ tag
				s->add_mdz_nucleotide(revtwobit_repr(b));


				newSeeds.emplace_back(s); // Push the new seed to the vector.
			}
		}
	}
}

void ReadAlignment::getInsertionSeeds(CountType base_repr, USeed origin, SeedVec & newSeeds) {

	// Don't handle insertions after deletions.
	if ( origin->cigar_data.back().offset == DELETION )
		return;

	// Don't handle if insertion region is getting too long
	if ( origin->cigar_data.back().offset == INSERTION && origin->cigar_data.back().length > globalAlignmentSettings.get_max_gap_length())
		return;

	// If all occurences of the seed match the current base, insertion is not required.
	FMTopDownIterator it(idx->idx, origin->vDesc);
	uint64_t range = origin->vDesc.range.i2 - origin->vDesc.range.i1;
	if ( seqan::goDown(it, seqan::DnaString(revtwobit_repr(base_repr))) ) {
		if ( it.vDesc.range.i2 - it.vDesc.range.i1 == range )
			return;
	}

	// handle all non-deletion regions
	if ( origin->cigar_data.back().offset != DELETION ) {

		// Compute new maximal alignment score when having an insertion
		ScoreType new_max_as = origin->max_as - globalAlignmentSettings.get_insertion_extension_penalty() - globalAlignmentSettings.get_match_score();
		if ( origin->cigar_data.back().offset != INSERTION)
			new_max_as -= globalAlignmentSettings.get_insertion_opening_penalty();

		// Extend insertion region as long as number of errors is allowed.
		if ( new_max_as >= getMinCycleScore(cycle, total_cycles) ) {

			// Don't go down the tree since an insertion means to stay at the same index position!

			// copy data from origin seed
			USeed s (new Seed);
			s->vDesc = origin->vDesc;
			s->cigar_data = origin->cigar_data;
			s->mdz_length = origin->mdz_length;
			s->mdz_nucleotides = origin->mdz_nucleotides;
			s->max_as = new_max_as; // Increase number of errors

			// Insert INSERTION region if necessary
			if ( s->cigar_data.back().offset != INSERTION )
				s->cigar_data.emplace_back(0, INSERTION);

			s->cigar_data.back().length += 1; // Increase insertion length

			newSeeds.emplace_back(s); // Push the new seed to the vector.
		}
	}
}

void ReadAlignment::getDeletionSeeds(CountType base_repr, USeed origin, SeedVec & newSeeds) {

	// Don't handle insertion regions
	if ( origin->cigar_data.back().offset == INSERTION )
		return;

	// iterate through possible bases
	for (int b=0; b<4; b++) {

		// Can only be opening (extension is performed in recursive_goDown() function)
		ScoreType new_max_as = origin->max_as - globalAlignmentSettings.get_deletion_opening_penalty() - globalAlignmentSettings.get_deletion_extension_penalty();

		// Don't start iteration with a match base (no deletion region required there!) or when having already all allowed errors.
		if ( b == base_repr || new_max_as < getMinCycleScore(cycle, total_cycles) )
			continue;

		// Create index iterator
		FMTopDownIterator it(idx->idx, origin->vDesc);

		// Only consider paths existing in the index
		if ( goDown(it, seqan::DnaString(revtwobit_repr(b)))) {

			// copy data from origin seed
			USeed s (new Seed);
			s->vDesc = it.vDesc;
			s->cigar_data = origin->cigar_data;
			s->mdz_length = origin->mdz_length;
			s->mdz_nucleotides = origin->mdz_nucleotides;
			s->max_as = new_max_as; // Increase number of errors

			s->cigar_data.emplace_back(1, DELETION); // init deletion region

			// Add MDZ nucleotide
			s->add_mdz_nucleotide(revtwobit_repr(b));

			recursive_goDown(base_repr, s, newSeeds);
		}
	}
}

void ReadAlignment::recursive_goDown(CountType base_repr, USeed origin, SeedVec & newSeeds) {

	if ( origin->max_as < globalAlignmentSettings.get_min_as() ) // too many errors -> no seeds
		return;

	// Only handle deletion regions
	if ( origin->cigar_data.back().offset != DELETION )
		return;

	// iterate through possible bases
		for (int b=0; b<4; b++) {

			// Create index iterator
			FMTopDownIterator it(idx->idx, origin->vDesc);

			// Only consider paths existing in the index
			if ( goDown(it, seqan::DnaString(revtwobit_repr(b)))) {

				// Finish deletion region when a match occurs
				if ( b == base_repr ) {

					// copy data from origin seed
					USeed s (new Seed);
					s->vDesc = it.vDesc;
					s->cigar_data = origin->cigar_data;
					s->mdz_length = origin->mdz_length;
					s->mdz_nucleotides = origin->mdz_nucleotides;
					s->max_as = origin->max_as ;

					s->cigar_data.emplace_back(1, 0); // init MATCH region
					newSeeds.emplace_back(s);

				}

				// Continue deletion region otherwise (only as long as it is shorter than maximal InDel length)
				else if ( origin->cigar_data.back().length < globalAlignmentSettings.get_max_gap_length() ){

					ScoreType new_max_as = origin->max_as - globalAlignmentSettings.get_deletion_extension_penalty();

					// stop if maximum number of errors reached
					if ( new_max_as < getMinCycleScore(cycle, total_cycles) )
						continue;

					// copy data from origin seed
					USeed s (new Seed);
					s->vDesc = it.vDesc;
					s->cigar_data = origin->cigar_data;
					s->mdz_length = origin->mdz_length;
					s->mdz_nucleotides = origin->mdz_nucleotides;
					s->max_as = new_max_as; // Increase number of errors

					s->cigar_data.back().length += 1; // extend deletion region
					s->add_mdz_nucleotide(revtwobit_repr(b));
					recursive_goDown(base_repr, s, newSeeds);

				}
			}
		}
}

void ReadAlignment::createSeeds(SeedVec & newSeeds) {

	// Stop if not enough cycles to create anchor
	if ( cycle < globalAlignmentSettings.get_anchor_length() )
		return;

	// Stop if not a seeding cycle
	if ( cycle != globalAlignmentSettings.get_anchor_length() && ( cycle - globalAlignmentSettings.get_anchor_length() ) % globalAlignmentSettings.get_seeding_interval() != 0 )
		return;

	CountType softclip_cycles = cycle - globalAlignmentSettings.get_anchor_length();
	CountType max_match_cycles = total_cycles - softclip_cycles;
	ScoreType max_as = ( max_match_cycles * globalAlignmentSettings.get_match_score() ) - getMinSoftclipPenalty(softclip_cycles);

	if ( max_as < getMinCycleScore(cycle, total_cycles) )
		return;

	std::string anchor_seq = getSequenceString().substr( sequenceLen - globalAlignmentSettings.get_anchor_length(), globalAlignmentSettings.get_anchor_length());

	FMTopDownIterator it(idx->idx);
	if ( seqan::goDown(it, seqan::DnaString(anchor_seq)) ) {
		USeed seed ( new Seed() );
		seed->max_as = max_as;
		if ( cycle != globalAlignmentSettings.get_anchor_length() )
			seed->cigar_data.emplace_back( ( cycle - globalAlignmentSettings.get_anchor_length()) , NO_MATCH ); // add softclip
		seed->cigar_data.emplace_back(globalAlignmentSettings.get_anchor_length(), 0);
		seed->vDesc = it.vDesc;
		newSeeds.emplace_back(seed);
	}

}

void ReadAlignment::extend_alignment(char bc) {

    // cycle is not allowed to be > total_cycles
    assert( total_cycles >= cycle );

    appendNucleotideToSequenceStoreVector(bc);

    // do not update the alignments when reading the first kmer_span-1 cycles
    if ( cycle < globalAlignmentSettings.get_anchor_length() )
        return;

	SeedVec newSeeds;

	// Reserve space for the worst case scenario: 4*M + 4*D + 1*I for each existing seed plus 1 new seed
	newSeeds.reserve(9*seeds.size() + 1);

	// extend existing seeds
    char base = revtwobit_repr(bc & 3);
    for ( auto seed = seeds.begin(); seed != seeds.end(); ++seed ) {
    	extendSeed(base, (*seed), newSeeds);
    }

    // create new seeds in defined intervals (intervals are handled inside createSeeds() function)
    if ( cycle >= globalAlignmentSettings.get_anchor_length() ) {
    	createSeeds(newSeeds);
    }

	// sort the newSeeds vector
    std::sort(newSeeds.begin(), newSeeds.end(), PComp<USeed>);

	// put only the "best" of similar seeds into the vector
	seeds.clear();

	if(newSeeds.size() == 0)
		return;

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

	return;
}

CountType ReadAlignment::getBarcodeIndex() {

	// Get the barcodes of the read
	std::string read_bc = getBarcodeString();

	if ( read_bc.length() == 0 )
		return NO_MATCH;

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
  if ( ! globalAlignmentSettings.get_keep_all_sequences() ) {
	  sequenceLen=0;
	  sequenceStoreVector.clear();
  }
}

bool ReadAlignment::is_disabled() {
	return flags == 0;
}

PositionType ReadAlignment::get_SAM_start_pos(GenomePosType p, USeed & sd) {

	// Only valid if CIGAR string exist (should always be the case)
	if ( sd->cigar_data.size() == 0 )
		return std::numeric_limits<PositionType>::max();

	// Retrieve position from the index
	PositionType sam_pos = p.pos;

	// REVERSE POSITIONS: position is already correct
	if ( idx->isReverse(p.gid) ) {
		return sam_pos;
	}

	// FORWARD POSITIONS: Compute start position
	else {
		sam_pos = idx->getSequenceLength(p.gid) - p.pos;

		CountType cigLen = 0;
		for ( auto el = sd->cigar_data.begin(); el != sd->cigar_data.end(); ++el ) {
			if ( el->offset != INSERTION ) {
				cigLen += el->length;
			}
		}
		sam_pos -= cigLen;

		// Begin pos in SAM specification ignores the softclip.
		sam_pos += sd->get_softclip_length();
	}

    return sam_pos;
}

void ReadAlignment::sort_seeds_by_as() {
	std::sort(seeds.begin(), seeds.end(), seed_comparison_by_as);
}

std::vector<uint8_t> ReadAlignment::getMAPQs(){

	// Vector that contains the final MAPQ values
	std::vector<uint8_t> mapqs;

	// True, if there is exactly one alignment for the read.
	bool unique = true;

	// Stop if not seeds.
	if ( seeds.size() == 0 )
		return mapqs;
	// Check if the best alignment can still be unique.
	else if ( seeds.size() != 1 )
		unique = false;

	// Vector that contains the MAPQ factor for each seed
	// The MAPQ factor is the weight of all alignments of a read
	std::vector<float> mapqFactors;

	uint16_t minSingleErrorPenalty = getMinSingleErrorPenalty();
	ScoreType maxPossibleScore = getMaxPossibleScore(cycle);

	// Maximal number of single errors for the minimal score of the current cycle (this does not consider affine gap penalties)
	float maxErrorsWithMinPenalty = float(maxPossibleScore - getMinCycleScore(cycle, total_cycles)) / float(minSingleErrorPenalty);

	// Maximal error percentage for the minimal score of the current cycle ( this does not consider affine gap penalties)
	float max_error_percent = float(100.0f * float(maxErrorsWithMinPenalty)) / float(total_cycles);

	mapqs.reserve(seeds.size());
	mapqFactors.reserve(seeds.size());

	// Sum of all MAPQ factors
	float mapqFactorsSum = 0;

	// Minimal number of errors for the best seed. This is required to compute the MAPQ factor
	float minReadErrors = float( getMaxPossibleScore(cycle) - seeds[0]->get_as() ) / getMinSingleErrorPenalty();

	// Compute the sum of all MAPQ factors
	for ( auto & s : seeds ) {

		// Weight the minimal number of errors
		float mapqFactor = std::pow ( 0.01f, (float( getMaxPossibleScore(cycle) - s->get_as() ) / getMinSingleErrorPenalty()) - minReadErrors );

		mapqFactors.push_back(mapqFactor);

		// For the sum, multiply the factor by the number of positions for this seed.
		mapqFactorsSum += ( mapqFactor * ( s->vDesc.range.i2 - s->vDesc.range.i1 ) );

		// Check if the best alignment can still be unique
		if ( s->vDesc.range.i2 - s->vDesc.range.i1 > 1 )
			unique = false;
	}

	// Compute the MAPQ for all seeds
	for ( CountType i=0; i<seeds.size(); i++ ) {

		// Always give 42 for unique, perfectly matching alignments
		if ( unique && seeds[i]->get_as() == maxPossibleScore ) {
			mapqs.push_back(42);
		}

		// Otherwise, apply MAPQ heuristics. Base Call Quality is not considered.
		else {

			float error_percent = float(100.0f * float(maxPossibleScore - seeds[i]->get_as()) / float(minSingleErrorPenalty)) / float(cycle);

			float prob = 1.0f;
			if ( seeds[i]->get_as() != maxPossibleScore ) {
				prob = 0.49995f + ( 0.5f * ( 1.0f - std::pow( 10.0f, -5.0f - max_error_percent + (2.0f * error_percent))));
			}

			mapqs.push_back(prob2mapq(prob * mapqFactors[i] / mapqFactorsSum));
		}
	}


	return mapqs;
}


