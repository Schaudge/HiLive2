#include "bamout.h"

std::vector<uint32_t> generate_bam_cigar(Seed &seed) {
    // Check if read was mapped to the forward sequence
    bool fw = (seed.start_pos >= 0); 
    
    // create a cigar vector in forward direction
    CigarVector cig;
    if (fw || seed.cigar_data.size() <= 1) {
        cig = seed.cigar_data;
    }
    else {
        // get iterator to last element != NO_MATCH: either last or second last element
        auto rit = seed.cigar_data.rbegin();
        if (rit->offset == NO_MATCH) {
          ++rit;  
        }

        int loffset = rit->offset;

        cig.reserve(seed.cigar_data.size());
        // reverse direction of 'cigar_data' and adjust offsets
        for(rit = seed.cigar_data.rbegin(); rit != seed.cigar_data.rend(); ++rit) {
            if (rit->offset != NO_MATCH)
                cig.emplace_back(rit->length,loffset - rit->offset);
            else
                cig.emplace_back(rit->length,NO_MATCH);
        }
    }
    
    // iterate over forward cigar vector and convert into BAM CIGAR structure
    std::vector<uint32_t> bam_cigar;
    int last_offset = 0;
    for (auto it = cig.begin(); it != cig.end(); ++it) {
      if (it->offset == NO_MATCH) {
        if (it+1 != cig.end()) {
          // mismatch in the middle
          int offset_change = (it+1)->offset - last_offset;
          if (offset_change == 0) {
            // mismatch
            bam_cigar.push_back( bam_cigar_gen(it->length - K + 1, BAM_CDIFF) );
          }
          else if (offset_change > 0) {
            // insertion in read
            bam_cigar.push_back( bam_cigar_gen(offset_change, BAM_CINS) );
            if (it->length > offset_change + K - 1) {
              // insertion plus mismatch
              bam_cigar.push_back( bam_cigar_gen(offset_change + K - 1 - it->length, BAM_CDIFF) );
            }
          }
          else {
            // deletion in read
            bam_cigar.push_back( bam_cigar_gen(-offset_change, BAM_CDEL) );
            if (it->length > K - 1) {
              // deletion plus mismatch
              bam_cigar.push_back( bam_cigar_gen(K - 1 - it->length, BAM_CDIFF) );
            }
          }          
        }
        else {
          // mismatch at the end of the read
          bam_cigar.push_back( bam_cigar_gen(it->length, BAM_CDIFF) );
        }
      }
      else {
        bam_cigar.push_back( bam_cigar_gen(it->length, BAM_CMATCH) );
        last_offset = it->offset;
      }
    }
   
    return bam_cigar;
}


BAMOut::BAMOut(std::string filename, KixRun* index) {
    
    // create the header object using the data from the index
    header = bam_header_init();

    header->n_targets = int32_t(index->num_seq); // number of reference sequences

    header->l_text = 0; // length of the plain text in the header
    
    uint32_t tlen [header->n_targets];
    for (int32_t i = 0; i < header->n_targets; i++){
        tlen[i] = uint32_t(index->seq_lengths[i]);
    }
    header->target_len = tlen; // lengths of the reference sequences
    
    char** tnames = new char* [header->n_targets];
    for (int32_t i = 0; i < header->n_targets; i++){
        tnames[i] = new char [index->seq_names[i].size() + 1];
        memcpy(tnames[i], index->seq_names[i].c_str(), index->seq_names[i].size() + 1);
    }
    header->target_name = tnames; // names of the reference sequences
    
    // TODO: fill the other header fields with information
    header->text = NULL; // plain text
    
    header->sdict = NULL; // header dictionary
    
    
    // open the BAM file and write header
    bamfile = bam_open(filename.c_str(), "wb");
    
    if (bamfile == 0) {
        std::cerr << "Error: could not open BAM file: " << filename << std::endl;
        throw;
    }
    BAMname = filename;

    bam_header_write(bamfile, header);
    
    
    // start the monitor in separate thread
    written = 0;
    failed = 0;
    running = true;
    accept_alignments = true;
    monitorThread = std::thread(&BAMOut::monitor,this);
}



BAMOut::~BAMOut() {
    // make the write function to reject alignments
    accept_alignments = false;
    
    // wait while there is data in the queue
    while (true) {
        // get lock
        std::lock_guard<std::mutex> lk1(m1);
        std::lock_guard<std::mutex> lk2(m2);
        // check size of the queue
        if ((alnQ1.size() == 0) && (alnQ2.size() == 0)) {
            // exit the loop to continue
            break;
        }
        else {
            // wait a bit before checking again
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
    
    // stop the monitor
    running = false;
    
    std::cout << "Wrote " << written << " alignments to " << BAMname << ". " << failed << " failed." << std::endl;   
    
    // close the file
    bam_close(bamfile);

    // destroy BAM header
    bam_header_destroy(header);
    

}


uint32_t BAMOut::write(ReadAlignment* read, uint8_t lane_id, uint32_t tile_id, uint32_t read_id, std::vector<char>* bc) {
    // check if the module accepts accepts alignments
    if (!accept_alignments) {
        return 0;
    }
    
    char readname [255];
    uint32_t l_qname = sprintf(readname, "L%d_T%d_R%d", lane_id, tile_id, read_id) +1;
 
    uint32_t l_qseq = 0;
    std::vector<uint8_t> nuc_seq;
    std::vector<uint8_t> qual_seq;
    if ((bc != nullptr) && (bc->size() > 0)) {
        // now convert the basecalls into the BAM-desired format
        l_qseq = bc->size();
        
        nuc_seq.resize((l_qseq+1)/2); // 4-bit encoding for nucleotides ("nybble")
        qual_seq.resize(l_qseq);
        
        std::vector<char>::iterator b_it = bc->begin();
        std::vector<uint8_t>::iterator n_it = nuc_seq.begin();
        std::vector<uint8_t>::iterator q_it = qual_seq.begin();
        
        while( b_it != bc->end() ) { // iterate over raw basecalls vector
            *q_it = ((*b_it >> 2) & 63); // basecall quality: get bits 3-8
            if (*b_it == 0) {
                // write N in case of no BC
                *n_it = 15 << 4;
            }
            else {
                // convert illumina nucleotide encoding into BAM encoding
                *n_it = 1 << (*b_it & 3) << 4;  // we are lucky here that both illumina and samtools use the ordering ACGT
            }
            // increase only basecalls and quality iterator, nucleotides proceed at half speed
            ++b_it; ++q_it;
            if (b_it != bc->end()) {
                // same procedure again
                *q_it = ((*b_it >> 2) & 63);
                if (*b_it == 0) {
                    *n_it |= 15;
                }
                else {
                    *n_it |= (1 << (*b_it & 3));
                }
                ++b_it; ++n_it; ++q_it;
            }
        }
    }
    
    // iterate over all seed positions in this read and push to queue one by one
    for (uint32_t i = 0; i < read->seeds.size(); i++) {
        // create as pointer, will be deleted when written to file
        bam1_t* aln = new bam1_t();

        // set unused data to default values
        aln->core.bin = 0;  // bin calculated by bam_reg2bin()
        aln->core.mtid = -1;  // chromosome ID of next read in template, defined by bam_hdr_t
        aln->core.mpos = -1; // 0-based leftmost coordinate of next read in template
        aln->core.isize = 0;  // ? (not documented)
        
        // chromosome ID, defined by bam_hdr_t
        aln->core.tid = read->seeds[i]->gid;  
        
        // 0-based leftmost coordinate
        aln->core.pos = read->get_SAM_start_pos(i);
        
        // mapping quality
        aln->core.qual = read->get_SAM_quality(i);
        
        // length of the query name
        aln->core.l_qname = l_qname;
        
        // bitwise flag
        aln->core.flag = read->get_SAM_flags(i);  
        
        // length of the query sequence (read) 
        aln->core.l_qseq = l_qseq;
        
        // process CIGAR information
        Seed sd = *(read->seeds[i]);
        std::vector<uint32_t> cigar_data = generate_bam_cigar(sd);
        uint32_t cigar_length = cigar_data.size();
        
        // number of CIGAR operations
        aln->core.n_cigar = cigar_length;

        aln->l_data = l_qname + cigar_length*4 + (l_qseq+1)/2 + l_qseq;
        aln->m_data = aln->l_data;
        
        uint8_t * data = new uint8_t[aln->l_data];
        
        uint64_t pos = 0;
        memcpy(data+pos, readname, l_qname); pos += l_qname;
        memcpy(data+pos, cigar_data.data(), cigar_length*4); pos += cigar_length*4;
        memcpy(data+pos, nuc_seq.data(), (l_qseq+1)/2); pos += (l_qseq+1)/2;
        memcpy(data+pos, qual_seq.data(), l_qseq); pos += l_qseq;
        
        aln->data = data; // all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
        
        // send data to one of the queues
        while(true) {
            if ( m1.try_lock() ) {
                alnQ1.push(aln);
                m1.unlock();
                break;
            }
            if ( m2.try_lock() ) {
                alnQ2.push(aln);
                m2.unlock();
                break;
            }
        }
    }
    
    delete read;
    return 0;
}


void BAMOut::monitor() {
    
    while (running) {
        if (m1.try_lock()) {
            // write the whole queue content to file
            while (alnQ1.size() > 0) {
                bam1_t* aln = alnQ1.front();
                alnQ1.pop();
                if (bam_write1(bamfile, aln))
                    written++;
                else
                    failed++;
                delete[] aln->data;
                delete aln;                
            }
            m1.unlock();
        }
        if (m2.try_lock()) {
            // write the whole queue content to file
            while (alnQ2.size() > 0) {
                bam1_t* aln = alnQ2.front();
                alnQ2.pop();
                if (bam_write1(bamfile, aln))
                    written++;
                else
                    failed++;
                delete[] aln->data;
                delete aln;                
            }
            m2.unlock();
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    
}


