#include "alnread.h"


seqan::String<seqan::CigarElement<> > Seed::returnSeqanCigarString() {
	typedef seqan::String<seqan::CigarElement<> > TSeqanCigarString;
	TSeqanCigarString seqanCigarString;
	seqan::CigarElement<> cigarElem;
	int last_offset = 0;
	for (CigarVector::const_iterator it = cigar_data.begin(); it != cigar_data.end(); ++it) { // for every element in cigar_data
		if (it == cigar_data.begin() && (*it).offset==NO_MATCH) { // Alignment begins with NO_MATCH region => Softclipped start 
			cigarElem.operation='S';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}
		if (it == --cigar_data.end() && (*it).offset==NO_MATCH) { // Alignment ends with NO_MATCH region => Softclipped end 
			cigarElem.operation='S';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}
		if ((*it).offset==NO_MATCH) { // ??? -> Mismatch
			cigarElem.operation='M';
			cigarElem.count=(*it).length;
			seqan::appendValue(seqanCigarString, cigarElem);
			continue;
		}
		if ((*it).offset!=NO_MATCH) { // ??? -> Match
            if (last_offset == (*it).offset) { // no offset change
                cigarElem.operation='M';
                cigarElem.count=(*it).length;
                seqan::appendValue(seqanCigarString, cigarElem);
                last_offset = (*it).offset;
                continue;
			}
			if (last_offset < (*it).offset) { // offset gets bigger => reference in alignment is longer than read => deletion in read (and thereby cigar string)
				cigarElem.operation='D';
				cigarElem.count=(*it).offset - last_offset;
				seqan::appendValue(seqanCigarString, cigarElem);
				cigarElem.operation='M';
				cigarElem.count=(*it).length;
				seqan::appendValue(seqanCigarString, cigarElem);
                last_offset = (*it).offset;
				continue;
			}
			if (last_offset > (*it).offset) { // offset gets smaller => reference in alignment is smaller than read => insertion in read (and thereby cigar string)
				cigarElem.operation='I';
				cigarElem.count=last_offset - (*it).offset;
				seqan::appendValue(seqanCigarString, cigarElem);
				cigarElem.operation='M';
				cigarElem.count=(*it).length;
				seqan::appendValue(seqanCigarString, cigarElem);
                last_offset = (*it).offset;
				continue;
			}
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
	
	// reverse Cigar String if seed maps reverse
    if (start_pos < 0)
        seqan::reverse(seqanCigarString);
    return seqanCigarString;
}

uint16_t Seed::serialize_size() {
  // calculate total size
  uint16_t total_size = 0;

  total_size += sizeof(GenomeIdType); // the target genome ID
  total_size += sizeof(PositionType); // the start position
  total_size += sizeof(CountType); // the number of matching positions
  
  if (cigar_data.size() >= 256)
    throw std::overflow_error("CIGAR information contains more than 255 elements!");
  uint8_t cigar_len = cigar_data.size();
  total_size += sizeof(uint8_t); // the size of the cigar information 
  total_size += cigar_len*(sizeof(CountType) + sizeof(DiffType)); // the cigar information itself
  return total_size;
}

std::vector<char> Seed::serialize() {
  // get the total size of the serialization
  uint16_t total_size = serialize_size();
  uint8_t cigar_len = (uint8_t) cigar_data.size();

  // create the vector to store the data
  std::vector<char> data (total_size);
  char* d = data.data();
  
  // write the target Genome ID
  memcpy(d,&gid,sizeof(GenomeIdType));
  d += sizeof(GenomeIdType);

  // write the start position
  memcpy(d,&start_pos,sizeof(PositionType));
  d += sizeof(PositionType);

  // write the number of matches
  memcpy(d,&num_matches,sizeof(CountType));
  d += sizeof(CountType);

  // write the number of cigar elements
  memcpy(d,&cigar_len,sizeof(uint8_t));
  d += sizeof(uint8_t);
  
  // write the seeds
  for (auto it = cigar_data.begin(); it != cigar_data.end(); ++it) {
    memcpy(d,&(it->length),sizeof(CountType));
    d += sizeof(CountType);

    memcpy(d,&(it->offset),sizeof(DiffType));
    d += sizeof(DiffType);
  }
  
  return data;
}



uint16_t Seed::deserialize(char* d) {
  
  // the total number of bytes read
  uint16_t bytes = 0; 
  
  // read the target Genome ID
  memcpy(&gid,d,sizeof(GenomeIdType));
  bytes += sizeof(GenomeIdType);

  // read the start position
  memcpy(&start_pos,d+bytes,sizeof(PositionType));
  bytes += sizeof(PositionType);

  // read the number of matches
  memcpy(&num_matches,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the number of cigar elements
  uint8_t cigar_len = 0;
  memcpy(&cigar_len,d+bytes,sizeof(uint8_t));
  bytes += sizeof(uint8_t);

  // read the cigar elements
  cigar_data.clear();
  for (uint8_t i = 0; i < cigar_len; ++i) {
    CigarElement cig;
    memcpy(&(cig.length),d+bytes,sizeof(CountType));
    bytes += sizeof(CountType);

    memcpy(&(cig.offset),d+bytes,sizeof(DiffType));
    bytes += sizeof(DiffType);

    cigar_data.emplace_back(cig);
  }

  return bytes;  
}


void ReadAlignment::set_total_cycles(CountType c) {
    total_cycles = c;
}


uint64_t ReadAlignment::serialize_size() {
  // calculate total size first
  uint64_t total_size = 0;
  total_size += 1; // the flag
  total_size += sizeof(CountType); // the cycle number
  total_size += sizeof(CountType); // the last_invalid cycle
  
  total_size += sizeof(CountType); // the sequence length
  total_size += sequenceStoreVector.size()*(sizeof(uint8_t)); // the sequence information itself
  total_size += sizeof(CountType); // The barcode length
  total_size += barcodeStoreVector.size()*(sizeof(uint8_t)); // the barcode sequence information

  // total number of seeds
  total_size += sizeof(uint32_t);

  // size of the single seeds
  for (auto & s : seeds) {
    total_size += sizeof(uint16_t) + s->serialize_size();
  }
  
  return total_size;
}


std::vector<char> ReadAlignment::serialize() {
  // get the total size of the serialization
  uint64_t total_size = serialize_size();
  uint32_t num_seeds = (uint32_t) seeds.size();

  // create the vector to store the data
  std::vector<char> data (total_size);
  char* d = data.data();
  
  // write the flag
  memcpy(d,&flags,1);
  d++;

  // write the cycle
  memcpy(d,&cycle,sizeof(CountType));
  d += sizeof(CountType);

  // write the last invalid cycle
  memcpy(d,&last_invalid,sizeof(CountType));
  d += sizeof(CountType);

  // write the sequence length
  memcpy(d,&sequenceLen,sizeof(CountType));
  d += sizeof(CountType);

  // write the sequenceStoreVector
  for (auto it = sequenceStoreVector.begin(); it != sequenceStoreVector.end(); ++it) {
    memcpy(d,&(*it),sizeof(uint8_t));
    d += sizeof(uint8_t);
  }

  // write the barcode length
  memcpy(d,&barcodeLen,sizeof(CountType));
  d += sizeof(CountType);

  // write the barcodeStoreVector
  for (auto it = barcodeStoreVector.begin(); it != barcodeStoreVector.end(); ++it) {
    memcpy(d,&(*it),sizeof(uint8_t));
    d += sizeof(uint8_t);
  }

  // write the number of seeds
  memcpy(d,&num_seeds,sizeof(uint32_t));
  d += sizeof(uint32_t);
  
  // write the seeds
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
  // the total number of bytes read
  uint64_t bytes = 0; 
  
  // read the flag
  memcpy(&flags,d,1);
  bytes++;

  // read the cycle
  memcpy(&cycle,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the last invalid cycle
  memcpy(&last_invalid,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the sequence length
  sequenceLen = 0;
  memcpy(&sequenceLen,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the sequence
  unsigned seqVec_size = (unsigned) std::ceil((float) sequenceLen / 4.0);
  sequenceStoreVector.clear();
  sequenceStoreVector.reserve(seqVec_size);
  for (unsigned i = 0; i <seqVec_size; ++i) {
    uint8_t elem;
    memcpy(&(elem),d+bytes,sizeof(uint8_t));
    bytes += sizeof(uint8_t);
    sequenceStoreVector.push_back(elem);
  }

  // read the barcode length
  barcodeLen = 0;
  memcpy(&barcodeLen,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the barcode
  unsigned barVec_size = (unsigned) std::ceil((float) barcodeLen / 4.0);
  barcodeStoreVector.clear();
  barcodeStoreVector.reserve(barVec_size);
  for (unsigned i = 0; i <barVec_size; ++i) {
    uint8_t elem;
    memcpy(&(elem),d+bytes,sizeof(uint8_t));
    bytes += sizeof(uint8_t);
    barcodeStoreVector.push_back(elem);
  }

  // read the number of seeds
  uint32_t num_seeds = 0;
  memcpy(&num_seeds,d+bytes,sizeof(uint32_t));
  bytes += sizeof(uint32_t);

  // read the seeds
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

    // insert into sorted list. The data read in is already sorted!!!
    // therefore I only push back
    seeds.push_back(std::move(s));
  }

  return bytes;  
}


// convert and return sequence of the seed as string
std::string ReadAlignment::getSequenceString() {
    std::string seq = "";
    // append one 4 base block at a time
    for (unsigned i = 0; i<sequenceStoreVector.size(); i++)
        seq.append(unhash(sequenceStoreVector[i], 4)); // 4 because 4 nucleotides fit in one element of sequenceStoreVector, namely a uint8_t

    // delete last few bases, because sequence is only len long
    for (unsigned i = sequenceLen; i<4*sequenceStoreVector.size(); ++i)
        seq.pop_back();

    // return sequence without barcode
    return seq;
}


std::string ReadAlignment::getBarcodeString() {
    std::string seq = "";
    // append one 4 base block at a time
    for (unsigned i = 0; i<barcodeStoreVector.size(); i++)
        seq.append(unhash(barcodeStoreVector[i], 4)); // 4 because 4 nucleotides fit in one element of sequenceStoreVector, namely a uint8_t

    // delete last few bases, because sequence is only len long
    for (unsigned i = barcodeLen; i<4*barcodeStoreVector.size(); ++i)
        seq.pop_back();

    // return barcode sequence
    return seq;
}


// append one nucleotide to sequenceStoreVector
void ReadAlignment::appendNucleotideToSequenceStoreVector(char nuc, bool appendToBarCode) {

	CountType & len = appendToBarCode ? barcodeLen : sequenceLen;
	std::vector<uint8_t> & seqVector = appendToBarCode ? barcodeStoreVector : sequenceStoreVector;

    // check if all bits from sequenceStoreVector are used
    if (len % 4 == 0) { // yes, all used => new 8 Bit block needs to be created
        uint8_t newBlock = twobit_repr(nuc);
        newBlock = newBlock << 6; // 'empty' bits are on the right side
        seqVector.push_back(newBlock);
        ++len;
    }
    else { // not all bits are used. At least the last two are unoccupied
        // append new 2 bits to the right of the old bits
        if (len % 4 == 1) // 6 bits empty
            seqVector.back() = seqVector.back() >> 4;
        if (len % 4 == 2) // 4 bits empty
            seqVector.back() = seqVector.back() >> 2;
        seqVector.back() = seqVector.back() | twobit_repr(nuc);
        ++len;
        // shift used bits until they are left aligned
        if (len % 4 == 2) // 4 bits empty
            seqVector.back() = seqVector.back() << 4;
        if (len % 4 == 3) // 2 bits empty
            seqVector.back() = seqVector.back() << 2;
    }
}


// helper function for add_new_seeds
bool seed_compare_pos (const USeed & i, const USeed & j) { 
  return (i->start_pos < j->start_pos); 
}


// Create new seeds from a list of kmer positions and add to current seeds
void ReadAlignment::add_new_seeds(GenomePosListType& pos, std::vector<bool> & posWasUsedForExtension, AlignmentSettings & settings) {
  SeedVecIt sit = seeds.begin();

  for(GenomePosListIt it = pos.begin(); it != pos.end(); ++it) {
    if (posWasUsedForExtension[it - pos.begin()]) // if current reference hit was used at least once for seed extension
        continue;
    USeed s (new Seed);

    s->gid = it->gid;
    s->start_pos = it->pos - (cycle-settings.kmer_span);
    s->num_matches = K_HiLive;
    s->cigar_data.clear();
    if (cycle-settings.kmer_span > 0)
      s->cigar_data.emplace_back(cycle-settings.kmer_span,NO_MATCH);

    // set correct matches and mismatches depending on kmer mask
    std::vector<unsigned> gapVec = settings.kmer_gaps;
    gapVec.push_back(settings.kmer_span+1);
    unsigned lastProcessedGapPosition = 0;
    for (unsigned gapIndex = 0, nextGap; gapIndex < gapVec.size(); ++gapIndex) {
        nextGap = gapVec[gapIndex];
        if (nextGap - lastProcessedGapPosition - 1 > 0)
            s->cigar_data.emplace_back(nextGap - lastProcessedGapPosition - 1,0);
        lastProcessedGapPosition = nextGap;
        s->cigar_data.emplace_back(1,NO_MATCH);
    }
    s->cigar_data.pop_back(); // remove tailing NO_MATCH from the for-loop

    // insert seed into sorted list of seeds
    while (sit != seeds.end() && seed_compare_pos(*sit, s)) // if seed exists and elem has larger starting position than (*sit)
        ++sit;
    sit = seeds.insert(sit, std::move(s));
  }
}


// Extend or create a placeholder seed for read with only trimmed matches
void ReadAlignment::create_placeholder_seed(AlignmentSettings & settings) {
	// if no seed exist, create one (only in the correct cycle and if parameter is true).
	if(seeds.size()==0){
		USeed s (new Seed);
		s->gid = TRIMMED;
		s->num_matches = 1;
		s->cigar_data.clear();
		s->cigar_data.emplace_back(settings.kmer_span,TRIMMED_MATCH);

        // this is always sorted, because there is only one seed now
        seeds.push_back(std::move(s));
	}
}


// convert a placeholder seed to a set of normal seeds
void ReadAlignment::convertPlaceholder(GenomePosListType& pos, AlignmentSettings & settings){
	if(seeds.size()!=1)
		return;
	auto s = seeds.begin();
	auto matches = (*s)->num_matches;
	if(cycle - settings.kmer_span > 0 && (*s)->gid==TRIMMED){
		seeds.clear();
		for(auto p = pos.begin(); p != pos.end(); ++p){
			USeed newSeed (new Seed());
		    newSeed->gid = p->gid;
		    newSeed->start_pos = p->pos - (cycle - settings.kmer_span);
		    newSeed->num_matches = matches;
		    newSeed->cigar_data.clear();
            newSeed->cigar_data.emplace_back(cycle - settings.kmer_span,NO_MATCH);

            // set correct matches and mismatches depending on kmer mask
            std::vector<unsigned> gapVec = settings.kmer_gaps;
            gapVec.push_back(settings.kmer_span); // without +1 because the last match will be inserted in extend_alignment
            unsigned lastProcessedGapPosition = 0;
            for (unsigned gapIndex = 0, nextGap; gapIndex < gapVec.size(); ++gapIndex) {
                nextGap = gapVec[gapIndex];
                if (nextGap - lastProcessedGapPosition - 1 > 0)
                    newSeed->cigar_data.emplace_back(nextGap - lastProcessedGapPosition - 1,0);
                lastProcessedGapPosition = nextGap;
                newSeed->cigar_data.emplace_back(1,NO_MATCH);
            }
            newSeed->cigar_data.pop_back(); // remove tailing NO_MATCH from the for-loop

            // insert into sorted list. The list is empty and pos vector is already sorted.
            // therefore I only push back
            seeds.push_back(std::move(newSeed));
		}
	}
}


// filter seeds based on filtering mode and q gram lemma. Also calls add_new_seeds.
void ReadAlignment::filterAndCreateNewSeeds(AlignmentSettings & settings, GenomePosListType & pos, std::vector<bool> & posWasUsedForExtension) {
    // TODO this function misbehaves when using demultiplexing. Reads might be preferred, kicking others only to get discarded due to barcode in the end.
    // compute possible remaining matches
    int possibleRemainingMatches = total_cycles - cycle + settings.kmer_span - 1;
    if (cycle == total_cycles)
        possibleRemainingMatches = 0;

    // compute threshold for number of matching bases a seed must already have to have a chance of getting printed in the end
    // adjust threshold via adapted q-gram-lemma
    int num_matches_threshold = total_cycles - settings.min_errors*(settings.kmer_span) - possibleRemainingMatches;

    if ((settings.all_best_hit_mode) || (settings.any_best_hit_mode)) {
        int max_num_matches = 0;
        for(SeedVecIt sd = seeds.begin() ; sd !=seeds.end(); ++sd)
            max_num_matches = std::max(max_num_matches, (int) (*sd)->num_matches);
        num_matches_threshold = std::max(num_matches_threshold, (int) max_num_matches - possibleRemainingMatches);
    }
    if (settings.all_best_n_scores_mode && settings.best_n > 0) {
        std::vector<CountType> num_matches_vector;
        for (SeedVecIt it=seeds.begin(); it!=seeds.end();++it)
            num_matches_vector.push_back((*it)->num_matches);
        std::sort(num_matches_vector.begin(), num_matches_vector.end(), std::greater<CountType>()); // sort in decreasing order
        unsigned numberOfUniques = std::unique(num_matches_vector.begin(), num_matches_vector.end()) - num_matches_vector.begin(); // uniques are in front of vector now
        if (numberOfUniques > settings.best_n) {
            CountType nth_best_num_matches = num_matches_vector[settings.best_n - 1];
            num_matches_threshold = std::max(num_matches_threshold, (int) nth_best_num_matches - possibleRemainingMatches);
        }
    }
    // else nothing has to be done for all-hit-mode
    

    // delete all seeds which do not reach threshold
    SeedVecIt it=seeds.begin();
    bool foundHit = false;
    while ( it!=seeds.end()) {
        if ((*it)->num_matches < num_matches_threshold)
            it = seeds.erase(it);
        else if (settings.discard_ohw && (cycle>settings.start_ohw) && ((*it)->num_matches <= K_HiLive)) // remove one-hit-wonders
            it = seeds.erase(it);
        else if (cycle == total_cycles && settings.any_best_hit_mode && foundHit)
            it = seeds.erase(it);
        else
            ++it;
        foundHit = true;
    }

    // if a new seed would have a chance then create it
    if (num_matches_threshold <= 1) // new seed would have 1 more match than expected by possibleRemainingMatches
        add_new_seeds(pos, posWasUsedForExtension, settings);
}


// updates cigar_data accordingly to a new matching kmer
void ReadAlignment::addMatchingKmer(USeed & s, DiffType offset, AlignmentSettings & settings) {

    s->cigar_data.emplace_back(1,offset);
    s->num_matches += 1;

    ////////////////////////////////////////////////////////////
    //// determine last occurred offset ////////////////////////
    int last_offset = 0;
    if ((*prev(prev(s->cigar_data.end()))).offset != NO_MATCH) // if last CigarElement is Match
        last_offset = (*prev(prev(s->cigar_data.end()))).offset;
    else // then the one before has to be a match
        last_offset = (*prev(prev(prev(s->cigar_data.end())))).offset;
    assert(last_offset != NO_MATCH);

    ////////////////////////////////////////////////////////////
    //// split last kmer-span bases in single CigarElements ////
    CigarVector::iterator it = --(s->cigar_data.end());
    unsigned summedLength = 1;
    while (summedLength < settings.kmer_span) {
        ++summedLength;
        --it;
        if ((*it).length > 1) {
            CigarElement elem1(1, (*it).offset);
            CigarElement elem2((*it).length-1, (*it).offset);
            s->cigar_data.insert(it, elem2);
            s->cigar_data.insert(it, elem1);
            it = s->cigar_data.erase(it); // points now at element after former it
            --it; // points now at elem1
        }
    }
    // it points now at kmer-spans-th last element of cigar_data list

    ////////////////////////////////////////////////////////////
    //// remove possibly inserted bases ////////////////////////
    if (offset < last_offset) {
        --it;
        for (int i=0; i<last_offset-offset; ++i) {
            assert((*it).offset == NO_MATCH);
            if ((*it).length > 1)
                --((*it).length);
            else
                it = --(s->cigar_data.erase(it));
        }
        ++it; // now again points at kmer-spans-th last element of cigar_data list
    }

    ////////////////////////////////////////////////////////////
    //// set matched bases as match ////////////////////////////

    CigarVector::iterator it_save = it; // for joining
    unsigned positionInKmer = 1;
    while (it != s->cigar_data.end()) {
        // if positionInKmer is not in kmer_gaps and cigar element was match
        if ((*it).offset == NO_MATCH && std::find(settings.kmer_gaps.begin(), settings.kmer_gaps.end(), positionInKmer) == settings.kmer_gaps.end()) {
            (*it).offset = offset;
            s->num_matches += 1;
        }
        ++it;
        ++positionInKmer;
    }

    ////////////////////////////////////////////////////////////
    //// join last kmer-span+1 CigarElements ///////////////////
    it = it_save;
    if (it == s->cigar_data.begin())
        ++it;
    while (it != s->cigar_data.end()) {
        if ((*prev(it)).offset == (*it).offset) {
            (*prev(it)).length += 1;
            it = s->cigar_data.erase(it);
        }
        else
            ++it;
    }
}


// Extend an existing seed (to be precise, extend the CIGAR vector / data).
bool ReadAlignment::extendSeed(USeed & s, DiffType offset, AlignmentSettings & settings){
	// Extend CIGAR for TRIMMED k-mers
	if ((offset == TRIMMED_MATCH) || ((offset == NO_MATCH) && (s->cigar_data.back().offset == TRIMMED_MATCH))){
		s->cigar_data.back().length += 1;
		if (s->cigar_data.back().offset != NO_MATCH){
			s->num_matches += 1;
			return true;
		}
		return false;
	}
    
    // if last region was match region
	if(s->cigar_data.back().offset != NO_MATCH){ 
        // MATCH -> NO_MATCH
        if(offset == NO_MATCH){
            // new mismatch region
            s->cigar_data.emplace_back(1,NO_MATCH);
            return false;
        }
        // MATCH -> MATCH  Extend match region with new match.
        else {
            int offset_change = offset - s->cigar_data.rbegin()->offset;
            // without any mismatch in between there cannot be a valid offset_change other than 0
            if (offset_change!=0) {
                // new mismatch region
                s->cigar_data.emplace_back(1,NO_MATCH);
                return false;
            }
            // else: extend current match region
            addMatchingKmer(s, offset, settings);
            return true;
        }
    }
    // so last region was mismatch region
    else {
        // NO_MATCH -> NO_MATCH  Extend mismatch region with new mismatch.
        if(offset == NO_MATCH){
            s->cigar_data.back().length += 1;
            return false;
        }
        // NO_MATCH -> MATCH
        else {
            assert((++(s->cigar_data.rbegin()))->offset != TRIMMED_MATCH && (++(s->cigar_data.rbegin()))->offset != NO_MATCH);
            int offset_change = offset - (++(s->cigar_data.rbegin()))->offset;
            // If there is an offset change, I need to have seen the appropriate mismatches before.
            if ( offset_change == 0
            || ((offset_change < 0) && (s->cigar_data.back().length >= -offset_change + settings.kmer_span - 1)) // Insertion in read
            || ((offset_change > 0) && (s->cigar_data.back().length >= settings.kmer_span - 1 )) ) { // Deletion in read
                addMatchingKmer(s, offset, settings);
                return true;
            } else {
                // criteria not fulfilled: extend existing mismatch area.
                s->cigar_data.back().length += 1;
                return false;
            }
        }
    }
}



void ReadAlignment::extend_alignment(char bc, KixRun* index, AlignmentSettings* settings) {

	// move to the next cycle
	cycle += 1;

    // cycle is not allowed to be > total_cycles
    assert( total_cycles >= cycle );

	// update the last k-mer
	uint8_t qual = ((bc >> 2) & 63); // get bits 3-8
	if ( (bc == 0) || (qual < settings->min_qual) ){ // no call if all 0 bits or quality below threshold
		last_invalid = last_invalid > cycle ? last_invalid : cycle; // TODO append an N as basecall? Could be a bad idea
	}
    unsigned mask = 3;
    if (flags != 0) // if read is valid
        appendNucleotideToSequenceStoreVector(revtwobit_repr(bc & mask)); // get the nucleotide as an actual character, disregarding the quality

    // do not update the alignments when reading the first kmer_span-1 cycles
    if (cycle < settings->kmer_span)
        return;

	// update the alignments
	GenomePosListType pos;
    std::vector<bool> posWasUsedForExtension;
	// if last kmer of read is not valid
	if (!( last_invalid+settings->kmer_span-1 < cycle )) {
		// write a NO_MATCH
        for (auto sit = seeds.begin(); sit != seeds.end(); ++sit)
            extendSeed(*sit, NO_MATCH, *settings);
	}
	else {

		// get all occurrences of last_kmer (fwd & rc) from index
        const std::string sequence = getSequenceString();
        std::string::const_iterator it_lastKmer = sequence.end() - settings->kmer_span;
        HashIntoType last_kmer = 0;
        hash_fw(it_lastKmer, sequence.end(), last_kmer, *settings);
		pos = index->retrieve_positions(sequence.substr(sequence.length()-settings->kmer_span), *settings);
        posWasUsedForExtension.resize(pos.size(), false);

		// check if the current k-mer was trimmed in the index
		if ( (pos.size() == 1) && ((*pos.begin()).gid == TRIMMED) ) {
			if(seeds.size()==0 && cycle==settings->kmer_span ){	// Create placeholder seed if first k-mer is trimmed
				create_placeholder_seed(*settings);
			} else
				// pretend that all existing seeds could be extended
				for(auto sd = seeds.begin() ; sd !=seeds.end(); ++sd)
					extendSeed(*sd, TRIMMED_MATCH, *settings);
			// clear the pos list so nothing bad happens in the next steps
			pos.clear();
		}

		// not trimmed in the index --> try to extend existing seeds
		else {
			// Convert placeholder seeds (if any) to "real" seeds
			if(seeds.size()==1 && (*(seeds.begin()))->gid==TRIMMED){
				convertPlaceholder(pos, *settings);
			}
			// find support for each candidate: iterate over seed candidates and positions simultaneously
			auto cPos1 = pos.begin(), cPos2 = pos.begin(); // sliding window [cPos1, cPos2)
      
			for (auto cSeed = seeds.begin(); cSeed!=seeds.end(); ++cSeed ) {
				PositionType seed_pos = (*cSeed)->start_pos + cycle -settings->kmer_span;

                // adjust the window in the position list
                while( (cPos1!=pos.end()) && (cPos1->pos < seed_pos - settings->window) )
                    ++cPos1;
                while( (cPos2!=pos.end()) && (cPos2->pos <= seed_pos + settings->window) )
                    ++cPos2;
        
				// search all positions in the window for the best matching extension of the seed
				DiffType best_offset = settings->window+1;  // set larger than search window
				GenomePosListIt best_match = cPos2; // set behind the last element of the window
				for(GenomePosListIt kmerHitIt = cPos1; kmerHitIt!=cPos2; ++kmerHitIt)
					if (kmerHitIt->gid == (*cSeed)->gid){
					    // if offset gets bigger => Deletion in read
						int offset = kmerHitIt->pos - seed_pos;
						if ((best_match==cPos2)||(abs(offset) < abs(best_offset))) {
							best_match = kmerHitIt;
							best_offset = offset;
						}
					}
        
				// check if a best match was found for this seed
				if (best_match != cPos2) {
						if(extendSeed(*cSeed, best_offset, *settings))
		            		// if pos was used as match, mark it so that later it does not get converted into a new seed
		            		posWasUsedForExtension[best_match-pos.begin()] = true;
                }
				else{
					// no position found to extend the current seed
					extendSeed(*cSeed, NO_MATCH, *settings);
				}
			} // END: for(seeds...)
		} // END: not trimmed
	} // END: if last kmer is valid

    filterAndCreateNewSeeds(*settings, pos, posWasUsedForExtension);

	return;
}

CountType ReadAlignment::getBarcodeIndex(AlignmentSettings* settings) {

	// Get the barcodes of the read
	std::string read_bc = getBarcodeString();

	uint16_t fragment_errors = 0;
	uint16_t fragment_pos = 0;
	uint16_t fragment_num = 0;
	uint16_t matching_bc = NO_MATCH;

	// Iterate through all user-defined (multi-)barcodes
	// That's quite complicated since the read barcodes are consecutive and the user barcodes are divided in vectors. // TODO: change that?
	for ( uint16_t barcodeIndex = 0; barcodeIndex < settings->barcodeVector.size(); barcodeIndex++ ) {

		// reset values for the barcode
		fragment_errors = 0;
		fragment_pos = 0;
		fragment_num = 0;
		matching_bc = barcodeIndex;

		// for each base of the read barcode
		for ( uint16_t nucl = 0; nucl < read_bc.length(); nucl++ ) {

			// reset values for each barcode fragment
			if ( fragment_pos >= (settings->barcodeVector[barcodeIndex])[fragment_num].length() ) {
				fragment_pos = 0;
				fragment_num += 1;
				fragment_errors = 0;
				assert( fragment_num < (settings->barcodeVector[barcodeIndex]).size() );
			}

			// compare nucleotides and increase the number of fragment errors if not equal
			if ( read_bc.at(nucl) != (settings->barcodeVector[barcodeIndex])[fragment_num].at(fragment_pos) ) {
				fragment_errors++;
			}

			// if too many errors in a fragment, break the loop for the barcode
			if ( fragment_errors > settings->barcode_errors[fragment_num] ) {
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



// disable this alignment, i.e. delete all seeds and set the last_invalid indicator to the
// end of the read. --> This read will not be aligned and consumes almost no space.
void ReadAlignment::disable() {
  last_invalid = total_cycles;
  seeds.clear();
  flags = 0;
  sequenceLen=0;
  sequenceStoreVector.clear();
}


// obtain start position of a seed according to SAM (leftmost) 
PositionType ReadAlignment::get_SAM_start_pos(USeed & sd, AlignmentSettings & settings) {
  PositionType pos = sd->start_pos;
  if (pos < 0) {
    if (sd->cigar_data.back().offset == NO_MATCH)
        pos = -pos - total_cycles + settings.kmer_span - (++sd->cigar_data.rbegin())->offset;
    else
        pos = -pos - total_cycles + settings.kmer_span - (sd->cigar_data.rbegin())->offset;
  }
  return pos;
}


// Calculate the mapping quality for all alignments of the read based on the other alignments and the number of matching positions.
int16_t MAPQ(const SeedVec &sv){
  return sv.size();
  }
