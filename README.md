HiLive - Live Mapping of Illumina reads
=======================================

Description
-----------

HiLive is a read mapping tool that maps Illumina HiSeq (or comparable) 
reads right in the moment when they are produced. This means, read mapping 
is finished as soon as the sequencer is finished.


Website
-------

The HiLive project website is https://gitlab.com/SimonHTausch/HiLive

There you can find the latest version of HiLive, source code, documentation,
and examples.


Installation
------------

If the binary distribution from the project website does not work for you,
you can still compile HiLive from source.

Make sure that the following dependences are installed:

 * cmake (>= 2.8)
 * boost (system, filesystem, program\_options)
 * zlib
 * lz4

---

You also need to download a specific seqan build. Navigate to a folder where
you want to store the seqan library and run:

    git clone https://github.com/seqan/seqan.git
    cd seqan && git checkout 9119aa6

---

Check out the HiLive source code from the project website and adjust the paths of the
seqan module in the file CMakeLists.txt in the hilive folder (line 31 and 32).
Then, compile HiLive with:

    cd [hilive-code]
    mkdir build && cd build
    cmake ..
    make

---

To compile HiLive with a different k-mer size than k=15 make the following
adjustment (here: k=10):

    cmake -DHiLive\_K 10 ..
    make


Usage
-----

HiLive has two components:

 * ``hilive-build``  builds the k-mer index of the reference genome
 * ``hilive``        the read mapper itself

---

#### Using hilive-build:

Building a k-mer index from FASTA file input.fa to output file input.fa.kix:

    hilive-build input.fa

Building an index from a large reference genome. Here is makes sense to use trimming,
i.e. removing k-mers from the index that occurr more than 1000 times (for example) in
the index. The index is written into the file trimmed.kix

    hilive-build -t 1000 -o trimmed.kix input.fa

---

#### Using hilive:

To map reads in a 100bp run using default settings:

    hilive /path/to/BaseCalls /path/to/index.kix 100 /path/to/outputFolder

Hint: to concatenate the resulting sam files of a lane to a single file using samtools you could use the following commands:

    for file in /path/to/outputFolder/L001/s_1_*.sam; do samtools view -bS $file > ${file%.sam}.bam; done
    samtools cat /path/to/outputFolder/L001/s_1_*.bam > /path/to/outputFolder/L001/allAlignments.bam
    samtools view -SH /path/to/outputFolder/L001/s_1_1101.sam > /path/to/outputFolder/L001/allAlignments.sam
    samtools view /path/to/outputFolder/L001/allAlignments.bam >> /path/to/outputFolder/L001/allAlignments.sam

---

#### Demultiplexing:

To map reads from multiplexed sequencing runs, you can provide HiLive with the barcode sequences from your Sample Sheet.
In default cases, barcode sequences are read after the (first) read, such that demultiplexing is carried out after the mapping is completed. 
If you use double indexing, please concatenate both indices in the correct order and provide them as one sequence. Please take care that the number of cycles is exactly the read length from your Sample Sheet plus that of your complete barcode sequence. All entered indices must be of the same length. To provide multiple indices, enter the -XXX argument for every barcode or barcode combination, e.g.:

	hilive /path/to/BaseCalls /path/to/index.kix 107 /path/to/outputFolder -XXX barcode1 -XXX barcode2 ...

Reads containing one of the given barcodes will be written to the resulting samfiles. The corresponding barcode sequence is stored in the BC-field of the samfile.

Hint: You can write the results of each sample into a separate files using the following command:

	for barcode in barcode1 barcode2 ... ; do grep "^@" allAlignments.sam > allAlignments_$barcode.sam && grep "BC:Z:$barcode" allAlignments.sam >> allAlignments_$barcode.sam; done

---

We recommend to adjust the numbers of threads used by HiLive with -n. If possible,
make use of all threads on the machine!

Please consult the project website for more details on the parameters!


License
-------

See the file LICENSE for licensing information.


Contact
-------

Please consult the HiLive project website for questions!

If this does not help, please feel free to consult:

 * Martin S. Lindner <martin (at) mail-lindner.de> (technical contact)
 * Bernhard Y. Renard <renardb (at) rki.de> (project head)

also see CONTRIBUTORS for a complete list of contributors and their contact information
