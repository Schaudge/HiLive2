#include <time.h>
#include <chrono>

#include "alnread.h"


bool seed_compare (Seed i,Seed j) { 
  return (i.start_pos < j.start_pos); 
}

bool seed_compare_pos (const USeed & i, const USeed & j) { 
  return (i->start_pos < j->start_pos); 
}

bool seed_compare_num_matches (const USeed & i, const USeed & j) { 
  return (i->num_matches < j->num_matches); 
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
  cigar_data.reserve(cigar_len);
  for (uint8_t i = 0; i < cigar_len; ++i) {
    CigarElement cig;
    memcpy(&(cig.length),d+bytes,sizeof(CountType));
    bytes += sizeof(CountType);

    memcpy(&(cig.offset),d+bytes,sizeof(DiffType));
    bytes += sizeof(DiffType);

    cigar_data.push_back(cig);
  }

  return bytes;  
}


/*
ReadAlignment& ReadAlignment::operator=(const ReadAlignment& other) {
    if(&other == this)
        return *this;
    
    rlen = other.rlen;
    last_kmer = other.last_kmer;
    last_invalid = other.last_invalid;
    cycle = other.cycle;
    
    // deep copy of seeds
    seeds.clear();
    seeds.reserve(other.seeds.size());
    for (auto sit = other.seeds.begin(); sit != other.seeds.end(); ++sit) {
        USeed s (new Seed);
        *s = **sit;
        seeds.push_back(std::move(s));
    }
    return *this;
}*/


void ReadAlignment::set_rlen(CountType r) {
    rlen = r;
}


uint64_t ReadAlignment::serialize_size() {
  // calculate total size first
  uint64_t total_size = 0;
  total_size += 1; // the flag
  total_size += sizeof(HashIntoType); // the last k-mer
  total_size += sizeof(CountType); // the cycle number
  total_size += sizeof(CountType); // the last_invalid cycle

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

  // write the last k-mer
  memcpy(d,&last_kmer,sizeof(HashIntoType));
  d += sizeof(HashIntoType);

  // write the cycle
  memcpy(d,&cycle,sizeof(CountType));
  d += sizeof(CountType);

  // write the last invalid cycle
  memcpy(d,&last_invalid,sizeof(CountType));
  d += sizeof(CountType);

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

  // read the last k-mer
  memcpy(&last_kmer,d+bytes,sizeof(HashIntoType));
  bytes += sizeof(HashIntoType);

  // read the cycle
  memcpy(&cycle,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the last invalid cycle
  memcpy(&last_invalid,d+bytes,sizeof(CountType));
  bytes += sizeof(CountType);

  // read the number of seeds
  uint32_t num_seeds = 0;
  memcpy(&num_seeds,d+bytes,sizeof(uint32_t));
  bytes += sizeof(uint32_t);

  // read the seeds
  seeds.clear();
  seeds.reserve(num_seeds);
  for (uint32_t i = 0; i < num_seeds; ++i) {
    uint16_t seed_size = 0;
    memcpy(&seed_size,d+bytes,sizeof(uint16_t));
    bytes += sizeof(uint16_t);
    
    std::vector<char> seed_data (seed_size,0);
    memcpy(seed_data.data(),d+bytes,seed_size);
    bytes += seed_size;
    
    USeed s (new Seed);
    s->deserialize(seed_data.data());

    seeds.push_back(std::move(s));
  }

  return bytes;  
}



// Create new seeds from a list of kmer positions and add to current seeds
void ReadAlignment::add_new_seeds(GenomePosListType& pos) {

  seeds.reserve(seeds.size() + pos.size());
  
  for(GenomePosListIt it = pos.begin(); it != pos.end(); ++it) {

    USeed s (new Seed);

    s->gid = it->gid;
    s->start_pos = it->pos - (cycle-K);
    s->num_matches = 1;
    s->cigar_data.clear();
    if (cycle-K > 0)
      s->cigar_data.emplace_back(cycle-K,NO_MATCH);
    s->cigar_data.emplace_back(K,0);

    seeds.push_back(std::move(s));

  }

}



void ReadAlignment::extend_alignment(char bc, KixRun* index, AlignmentSettings* settings) {

  // q-gram-lemma: seeds added after last_new_seed have at least more than min_errors errors
  CountType last_new_seed = K*(settings->min_errors+1); 

  // move to the next cycle
  cycle += 1;

  // update the last k-mer
  uint8_t qual = ((bc >> 2) & 63); // get bits 3-8
  if ( (bc == 0) || (qual < settings->min_qual) ){ // no call if all 0 bits or quality below threshold 
    last_kmer = 0;
    last_invalid = cycle;
  }
  else{
    // update the current k-mer of the read using the basecall
    update_kmer(last_kmer, bc);
  }

  // update the alignments
  if ( last_invalid+K-1 < cycle ) {

    // get all occurrences of last_kmer (fwd & rc) from index
    GenomePosListType pos = index->retrieve_positions(last_kmer);

    // pos MUST be sorted. However, pos is sorted as long as the index is sorted (should be by default)
    // DEPRECATED
    if (settings->sort_positions) {
      std::sort(pos.begin(),pos.end(),gp_compare);
    }
    
    // maximum number of k-mer matches in seeds
    CountType max_num_matches = 0;
    
    // check if the current k-mer was trimmed in the index
    if ( (pos.size() == 1) && (pos[0].gid == TRIMMED) ) {
      // clear the pos list so nothing bad happens in the next steps
      pos.clear();
      
      // pretend that all existing seeds could be extended
      for (auto sd = seeds.begin() ; sd!=seeds.end(); ++sd ) {
        (*sd)->cigar_data.back().length += 1;
        if ((*sd)->cigar_data.back().offset != NO_MATCH) {
          (*sd)->num_matches += 1;
          max_num_matches = std::max(max_num_matches, (*sd)->num_matches);
        }
      }
    }
    // not trimmed in the index --> try to extend existing seeds
    else {
      // find support for each candidate: iterate over seed candidates and positions simultaneously
      auto cPos1 = pos.begin(); // beginning of the sliding window [cPos1, cPos2)
      auto cPos2 = pos.begin(); // end of the sliding window
      auto wSeed1 = seeds.begin(); // beginning of the sliding window in the seeds [wSeed1, wSeed2)
      auto wSeed2 = seeds.begin(); // end of the sliding window
      
      for (auto cSeed = seeds.begin(); cSeed!=seeds.end(); ++cSeed ) {
        PositionType seed_pos = (*cSeed)->start_pos + cycle -K;
        
        // adjust the window in the position list
        while( (cPos1!=pos.end()) && (cPos1->pos < seed_pos - settings->window) ){
          ++cPos1;
        }
        while( (cPos2!=pos.end()) && (cPos2->pos < seed_pos + settings->window) ){
          ++cPos2;
        }
        
        // adjust the neighboring seeds window
        while( (wSeed1!=seeds.end()) && ((*wSeed1)->start_pos < (*cSeed)->start_pos - 2*settings->window) ){
          ++wSeed1;
        }
        while( (wSeed2!=seeds.end()) && ((*wSeed2)->start_pos < (*cSeed)->start_pos + 2*settings->window) ){
          ++wSeed2;
        }
        
        // search all positions in the window for the best matching extension of the seed
        DiffType best_distance = settings->window+1;  // set larger than search window
        GenomePosListIt best_match = cPos2; // set behind the last element of the window
        for(GenomePosListIt win=cPos1; win!=cPos2; ++win){
          if (win->gid == (*cSeed)->gid){
            int dist = seed_pos - win->pos; 
            if ((best_match==cPos2)||(abs(dist) < abs(best_distance))) {
              best_match = win;
              best_distance = dist;
            }
          }
        }
        
        // check if a best match was found for this seed
        if (best_match != cPos2) {
          // find the best seed from the perspective of best_match
          DiffType best_sdist = 2*settings->window+1;  // set larger than search window
          auto best_seed = wSeed2;   // set behind the last element of the window
          for(auto win=wSeed1; win!=wSeed2; ++win){
            if ((*win)->gid == best_match->gid){
              int dist = best_match->pos - ((*win)->start_pos+cycle-K); 
              if ((best_seed==wSeed2)||(abs(dist) < abs(best_sdist))) {
                best_seed = win;
                best_sdist = dist;
              }
            }
          }
        
          if (best_seed == cSeed) {
            // the current seed is the best extension
            
            if ( (*cSeed)->cigar_data.back().offset == NO_MATCH ) {
              // Before starting a new matching area check if the match is valid
              // A match after a NO_MATCH area is valid if:
              //  - offset change == 0 and NO_MATCH length >= K (mismatch)
              //  - offset change > 0 and NO_MATCH length >= Offset change + K - 1 (insertion in read)
              //  - offset change < 0 and NO_MATCH length >= K - 1 (deletion in read)
              if ((*cSeed)->cigar_data.size() > 1) { 
                int offset_change = best_distance - ((*cSeed)->cigar_data.rbegin() + 1)->offset;
                if ( ((offset_change == 0) && ((*cSeed)->cigar_data.back().length >= K))
                    || ((offset_change > 0) && ((*cSeed)->cigar_data.back().length >= offset_change + K - 1))
                    || ((offset_change < 0) && ((*cSeed)->cigar_data.back().length > K - 1 )) ) {
                  // criteria fulfilled. Start a new match area. 1 matching k-mer = K matches
                  (*cSeed)->cigar_data.emplace_back(K,best_distance);  
                  (*cSeed)->num_matches += 1;
                  
                  // remove assigned position from the list of possible matches
                  if(best_match == cPos1){
                    cPos1 = pos.erase(best_match);
                  }
                  else {
                    pos.erase(best_match);
                  }
                  --cPos2; 
                }
                else {
                  // continue existing mismatch area
                  (*cSeed)->cigar_data.back().length += 1;
                }
              }
              else {
                // start a new match area. 1 matching k-mer = K matches
                (*cSeed)->cigar_data.emplace_back(K,best_distance);
                (*cSeed)->num_matches += 1;
                
                // remove assigned position from the list of possible matches
                if(best_match == cPos1){
                  cPos1 = pos.erase(best_match);
                }
                else {
                  pos.erase(best_match);
                }
                --cPos2; 
              }
            }
            else {
              // continue existing match area
              (*cSeed)->cigar_data.back().length += 1;
              (*cSeed)->num_matches += 1;
              
              // remove assigned position from the list of possible matches
              if(best_match == cPos1){
                cPos1 = pos.erase(best_match);
              }
              else {
                pos.erase(best_match);
              }
              --cPos2; 
            }
            max_num_matches = std::max(max_num_matches, (*cSeed)->num_matches);

          }
          else{
            // best match has another favourite, don't extend this seed
            if ( (*cSeed)->cigar_data.back().offset == NO_MATCH ) {
              // continue existing mismatch area
              (*cSeed)->cigar_data.back().length += 1;
            }
            else {
              // start new mismatch area
              (*cSeed)->cigar_data.emplace_back(1,NO_MATCH);
            }
          }
        }
        else{
          // no position found to extend the current seed
          if ( (*cSeed)->cigar_data.back().offset == NO_MATCH ) {
            // continue existing mismatch area
            (*cSeed)->cigar_data.back().length += 1;
          }
          else {
            // start new mismatch area
            (*cSeed)->cigar_data.emplace_back(1,NO_MATCH);
          }
        }
      } // END: for(seeds...)
    } // END: not trimmed

    
    // get the num_matches of the N'th best seed
    // seed list gets partially sorted (and sorted back later). Noone ever said it will be fast...
    CountType nth_best_matches = 0;
    if ( settings->best_n_mode && (settings->best_n > 0) && (seeds.size() >= settings->best_n) ) {
      std::nth_element(seeds.begin(), seeds.begin()+seeds.size()-settings->best_n , seeds.end(), seed_compare_num_matches);
      nth_best_matches = (*(seeds.begin()+seeds.size()-settings->best_n))->num_matches;
    }

    // set the last_new_seed cycle according to the mapping mode
    if ( settings->best_hit_mode ) {
      last_new_seed = std::min(CountType(rlen-max_num_matches+1),last_new_seed);
    }
    else if ( settings->best_n_mode ) {
      last_new_seed = std::min(CountType(rlen-nth_best_matches+1),last_new_seed);
    }

    // create new seed candidates for each k-mer match that was not used to extend a seed
    if ( cycle <= last_new_seed ) {

      add_new_seeds(pos);

      if (pos.begin() != pos.end()) {
        max_num_matches = std::max(max_num_matches, (CountType)1);
      }

    }

    // define a lambda function implementing all discard criteria.
    // all criteria must be fulfilled to keep the seed. Go from stronger to weaker criteria.
    auto crit = [&] (USeed & s) {
      // don't filter seeds that were extended in this cycle
      if (s->cigar_data.back().offset != NO_MATCH) {
        return false;
      }

      // 1. remove one-hit-wonders
      if ( settings->discard_ohw && (cycle>settings->start_ohw)&&(s->num_matches<=1) ) {
        return true;
      }

      // 2. remove according to q-gram lemma
      if ( cycle > (K*(settings->min_errors+1) + s->num_matches) ) {
        return true;
      }

      // 3. remove according Best-Hit-criteria
      if ( settings->best_hit_mode ) {
        if (cycle > (rlen - max_num_matches + s->num_matches)) {
          return true;
        }
      }
      // 4. remove according Best-N-criteria
      else if ( settings->best_n_mode ) {
        if ( cycle > (rlen - nth_best_matches + s->num_matches) ) {
          return true;
        }
      }

      // get the number of mismatches since the last match
      int since_last_match = s->cigar_data.back().length;

      // 5. heuristic criterium
      if ((since_last_match >= K+10)&&(s->num_matches < (int)(std::sqrt(cycle-K+1)))){
        return true;
      }

      return false;
    };

    seeds.erase(std::remove_if(seeds.begin(),seeds.end(),crit) , seeds.end());

    std::sort(seeds.begin(), seeds.end(), seed_compare_pos);

  } // END: if ( last_invalid+K-1 < cycle ) ...
  else {
    
    // write a NO_MATCH if cycle > K-1
    if ( cycle > K-1 ) {
      for (auto sit = seeds.begin(); sit != seeds.end(); ++sit){
        if ( (*sit)->cigar_data.back().offset == NO_MATCH ) {
          // continue existing mismatch area
          (*sit)->cigar_data.back().length += 1;
        }
        else {
          // start new mismatch area
          (*sit)->cigar_data.emplace_back(1,NO_MATCH);
        }
      }
    }
  }

  return ;
}



// disable this alignment, i.e. delete all seeds and set the last_invalid indicator to the
// end of the read. --> This read will not be aligned and consumes almost no space.
void ReadAlignment::disable() {
  last_invalid = rlen;
  seeds.clear();
  flags = 0;
}


// generate the SAM flags for a seed
uint32_t ReadAlignment::get_SAM_flags(uint32_t sd) {
  if ( sd < seeds.size() ) {
    uint32_t flags = 0;
    // flag for the reverse strand alignment
    if (seeds[sd]->start_pos < 0)
      flags += 16;

    // Primary/Secondary alignment
    // the "Primary" alignment is the alignment with the highest number of matches
    // in case of a tie, the left-most alignment wins (including negative positions)
    bool primary = true;
    for(uint32_t i = 0; i < seeds.size(); i++) {
      if(seeds[i]->num_matches > seeds[sd]->num_matches || (seeds[i]->num_matches == seeds[sd]->num_matches && i < sd)) {
	primary = false;
      }
    }
    if (!primary)
      flags += 256;
    return flags;
  }
  else {
    throw std::length_error(std::string("Error generating SAM flags: requested alignment ID (")+std::to_string(sd)+std::string(") exceeds number of alignments (")+std::to_string(seeds.size())+std::string(")!"));
  }
}


// obtain start position of a seed according to SAM (leftmost) 
PositionType ReadAlignment::get_SAM_start_pos(uint32_t sd) {
  PositionType pos = seeds[sd]->start_pos;
  if (pos < 0) {
    pos = -pos - rlen + K;
  }
  return pos;
}


// calculate a quality score according to SAM
uint16_t ReadAlignment::get_SAM_quality(uint32_t sd) {
  CountType best_match = 0;
  for (uint32_t s = 0; s < seeds.size(); s++) {
    best_match = std::max(seeds[s]->num_matches, best_match);
  }
  double total_weighted_matches = 0.;
  for (uint32_t s = 0; s < seeds.size(); s++) {
    total_weighted_matches = std::pow(4., double(seeds[s]->num_matches - best_match));
  }

  double prob = std::pow(4., double(seeds[sd]->num_matches - best_match)) / total_weighted_matches;
  int score = -10 * std::log10(prob);
  if (score > 255 || score < 0)
    score = 255;
  return score;
}










/* Calculate the mapping quality for all alignments of the read based on the other alignments
   and the number of matching positions.
 */

int16_t MAPQ(const SeedVec &sv){
  return sv.size();
  }


/* Helper function; checks if the position at 'sit' is sane. only applicable for reads matching 
   to the backwards genome sequence (forward sanity check is MUCH easier) */
bool is_sane(std::vector<DiffType>::iterator sit,
	     const std::vector<DiffType> &matches){
  DiffType dist = 0;
  DiffType this_offset = *sit;
  // find the OFFSET of the next matching region (after the next NOMATCH region)
  bool passed_NOMATCH = false;
  while(!((*sit!=NO_MATCH)&&passed_NOMATCH) && (sit!=matches.end())){
    if(*sit == NO_MATCH) {passed_NOMATCH = true;}
    ++sit;
    ++dist;
  }
  
  if (sit == matches.end()) {
    return true;
  }
  else {
    int offset_change = *sit - this_offset;
    return ((dist-K+1 - offset_change*(offset_change>0)) > 0);
  }
}

/* Construct a CIGAR string from a seed */
std::string CIGAR(const Seed &seed){
  // Alignments are always reported in forward direction wrt the reference sequence --> take care
  bool fw = (seed.start_pos >= 0); // Read was mapped to the forward sequence
  CigarVector cig;
  if (fw) {
    cig = seed.cigar_data;
  }
  else {
    // need to adjust the offsets --> find last offset != NO_MATCH
    auto rit = seed.cigar_data.rbegin();
    while((rit != seed.cigar_data.rend()) && (rit->offset == NO_MATCH)) {
      ++rit;
    }
    int loffset;
    if (rit != seed.cigar_data.rend())
      loffset = rit->offset;
    else
      loffset = 0;
    cig.reserve(seed.cigar_data.size());
    // reverse direction of 'cigar_data' and adjust offsets
    for(rit = seed.cigar_data.rbegin(); rit != seed.cigar_data.rend(); ++rit) {
      if (rit->offset != NO_MATCH)
	cig.emplace_back(rit->length,loffset - rit->offset);
      else
	cig.emplace_back(rit->length,NO_MATCH);
    }
  }
  
  // Now, calculate the CIGAR string
  std::string cigar;
  if (cig.size() > 0) {
    int previous_offset = 0; // the offset in the previous matching region
    for (uint32_t i = 0; i < cig.size(); i++) {
      // extend the cigar string accordingly
      if ( cig[i].offset != NO_MATCH ) {
	cigar.append(std::to_string(cig[i].length)+std::string("M"));
	previous_offset = cig[i].offset;
      }
      else {
	int offset_change = 0;
	int mm_len = cig[i].length -K+1;
	// correct the length of the mismatch region if there are insane positions
	if ( i+1 < cig.size() ) {
	  offset_change = cig[i+1].offset - previous_offset;
	  if (offset_change > mm_len) {
	    // offset change cannot be larger than the mismatch area
	    cig[i].length += (offset_change - mm_len);
	    cig[i+1].length -= (offset_change - mm_len);
	    mm_len = cig[i].length -K+1;
	  }
	}

	if ( offset_change == 0 ) {
	  // Sequence mismatch(es)
	  if ( i > 0 && i < cig.size()-1 ) {
	    cigar.append(std::to_string(cig[i].length -K+1)+std::string("X"));
	  }
	  else {
	    // report as soft-clipped at the beginning/end of the read
	    cigar.append(std::to_string(cig[i].length)+std::string("S"));
	  }
	}
	else {
	  // Number of insertions
	  int insertions = cig[i].length -K+1;
	  // Number of deletions
	  int deletions = - (offset_change - insertions);
	  // append to cigar string
	  if ( insertions > 0 )
	    cigar.append(std::to_string(insertions)+std::string("I"));
	  if ( deletions > 0 )
	    cigar.append(std::to_string(deletions)+std::string("D"));
	}
      }
    }
  }

  /*int count = K-1; // number of matches/mismatches in a row
    bool previous_match = true; // was the previous k-mer a match? -> recognize changes
    // change initial values if the read starts with a mismatch
    if(*(matches.begin()) == NO_MATCH) {
      count = K-1; 
      previous_match = false; 
    }
    int previous_offset = 0; // the offset in the previous matching region
    int pos = 0; // current position in the read
    for(auto it = matches.begin(); it != matches.end(); ++it){
      if ((*it) != NO_MATCH) {
	// The k-mer could be matched
	if (previous_match) {
	  // Check the sanity of this match, if found on backward seq
	  bool sane = true;
	  if (!fw){
	    sane = is_sane(it, matches);
	  }
	  if (sane) {
	    // Continue a match
	    count += 1;
	    previous_offset = (*it);
	  }
	  else {
	    // end the match here, treat this position as mismatch
	    //std::cout << "INSANE position!!!" << std::endl;
	    // start a new mismatch area here. Don't care about countinuing, this happens below
	    // Finish a match region
	    cigar.append(std::to_string(count)+std::string("M"));
	    // Start a new mismatch region
	    count = 1; 
	    previous_match = false;
	  }
	}
	// the previous k-mer was not matched
	else {
	  int offset_change = (*it) - previous_offset;
	  // check the sanity of the match
	  bool sane = ((count-K+1 - offset_change*(offset_change>0)) >= 0) && fw;
	  // add the criterium for backward matches here (as above)
	  sane = sane || (!fw && is_sane(it, matches));
	  if (sane) {
	    // Finish a mismatch region
	    if (offset_change == 0) {
	      if(count > K-1){
		// Sequence mismatch(es)
		cigar.append(std::to_string(count-K+1)+std::string("X"));
	      }
	    }
	    else {
	      // cure only on forward mapping reads
	      if (count > K-1){
		// Insertion compared to reference
		cigar.append(std::to_string(count-K+1)+std::string("I"));
	      }
	      if (count-offset_change > K-1){ 
		// In this case, a string of length count-offset_change was deleted from reference
		cigar.append(std::to_string(count-offset_change-K+1)+std::string("D"));
	      }
	    }
	    // Start a new match region. The first k-mer represents K matching positions
	    count = K;
	    previous_match = true;
	    previous_offset = (*it);
	  }
	  else {
	    // insane match --> treat as mismatch
	    count += 1;
	    previous_match = false;
	  }
	}
      }
      else {
	// The k-mer could NOT be matched
	if (!previous_match) {
	  // Continue a mismatch
	  count += 1;
	} else {
	  // Finish a match region
	  cigar.append(std::to_string(count)+std::string("M"));
	  // Start a new mismatch region
	  count = 1; 
	  previous_match = false;
	}
      }
      // increase the position counter
      pos++;
    }
    // Complete the cigar string
    if (previous_match) {
      // Fill cigar string with matches
      cigar.append(std::to_string(count)+std::string("M"));
    }
    else {
      // Fill cigar string with mismatches
      // TODO: probably report this as "soft-clipping"
      cigar.append(std::to_string(count)+std::string("X"));
    }
    }*/
  return cigar;
}
