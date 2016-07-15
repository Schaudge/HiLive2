#ifndef BAMOUT_H
#define BAMOUT_H

#include "headers.h"
#include "definitions.h"
#include "alnread.h"

#include "sam.h"

//-------------------------------------------------------------------//
//------  The BAMOut class to write BAM files  ----------------------//
//-------------------------------------------------------------------//

class BAMOut {
    private:
        // use two internal queues to buffer incoming alignments
        std::queue<bam1_t*> alnQ1;
        std::queue<bam1_t*> alnQ2;
        
        // mutex to ensure that only one process can access the queue at once
        std::mutex m1;
        std::mutex m2;
        
        // pointer to the header data structure
        bam_hdr_t* header;
        
        // BAM file handle
        bamFile bamfile;
        
        // count the number of written and failed alignments
        uint64_t written;
        uint64_t failed;
        
        // store the filename for messages
        std::string BAMname;
        
        // switch to tell the write function when alignments are accepted for writing
        bool accept_alignments;
        
        // switch to tell the monitor to stop running
        bool running;
        
        // check the current status of the queue and write alignments if available
        std::thread monitorThread;
        void monitor ();

    public:
        // open new BAM file for writing. Create header from kmer index
        BAMOut (std::string filename, KixRun* index);
        
        // close BAM file and make sure that everything is written properly
        ~BAMOut ();
        
        // write an alignment. Take care of concurrent access. Returns bytes written
        uint32_t write (ReadAlignment* read, uint8_t lane_id, uint32_t tile_id, uint32_t read_id, std::vector<char>* bc);
        
}; // class BAMOut










#endif /* BAMOUT_H */