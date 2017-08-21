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

If you are using a Debian based system you can directly install the Debian 
package from the official [Debian repository](https://packages.debian.org/sid/hilive "HiLive Debian package")

If this does not work for you, you can still compile HiLive from source.

Make sure that the following dependencies are installed:

 * cmake (>= 2.8)
 * boost (system, filesystem, program\_options)
 * zlib
 * lz4

If using a local version of lz4 then adjust path in CMakeLists.txt line 32.

---

You also need to download download seqan.
Cloning the repository makes it possible to switch between different versions:

    git clone https://github.com/seqan/seqan.git

HiLive 1.0 was tested with SeqAn version 2.3.2.

---

Check out the HiLive source code from the project website and adjust the paths of the
seqan module in the file CMakeLists.txt in the hilive folder (line 41 and 42).
Then, compile HiLive with:

    cd [hilive-code]
    mkdir build && cd build
    cmake ..
    make

Usage
-----

HiLive has three components:

 * ``hilive-build``  builds the k-mer index of the reference genome
 * ``hilive``        the read mapper itself
 * ``hilive-out``    executable to produce output files

---

#### Using hilive-build:

Building a k-mer index from FASTA file input.fa to output file input.fa.kix with k-mer weight 15:

    hilive-build input.fa 15

Building an index from a large reference genome. Here is makes sense to use trimming,
i.e. removing k-mers from the index that occurr more than 1000 times (for example) in
the index. The index is written into the file trimmed.kix

    hilive-build -t 1000 -o trimmed.kix input.fa 15
    
For gapped k-mers, use the -p parameter to specify the gap positions.
For example, for a gap pattern of 1101110011, type:

    hilive-build -p 3 7 8 input.fa 7
    
With the current index structure, we strongly recommend to use a maximum k-mer weight of 15 because of huge disk space and memory requirements for large 
k-mers.

---

#### Using hilive:

To map reads in a 100bp run using default settings:

    hilive /path/to/BaseCalls /path/to/index.kix 100 /path/to/outputFolder

For an overview of additional parameters, type 
	
	hilive --help

To prevent errors during argument parsing we recommend to set optional parameters AFTER the positional options:

	hilive BC_DIR INDEX CYCLES OUTDIR [options]
	
However, if unexpected parsing errors occur, please try to specify all parameters with the "--"-syntax (e.g. --BCDIR /path/to/BaseCalls) instead of using positional arguments. This is also necessary when loading (some of the) positional arguments from a settings file instead of using the command line.

---

#### Using hilive-out:

To create a SAM or BAM alignment output from existing temporary files in HiLive, type:

	hilive-out --settings /path/to/temp/dir/hilive_settings.xml

This will output the alignment results of the last cycle based on the settings that were specified for the related HiLive run.
To produce output files for other cycles, e.g. 50, 70, 90, type:

	hilive-out --settings /path/to/temp/dir/hilive_settings.xml --output-cycles 50 70 90
	
Please note, that the temporary files for the respective cycles must be present in the temp folder.
This is only the case if the --keep-files parameter and/or the --output-cycles parameter for the respective cycles was activated in the HiLive run.


#### Demultiplexing:

To map reads from multiplexed sequencing runs, you can provide HiLive with the barcode sequences from your Sample Sheet.
In default cases, barcode sequences are read after the (first) read, such that demultiplexing is carried out after the mapping is completed. 
If you use double indexing, please concatenate both indices in the correct order and provide them as one sequence. Please take care that the number of cycles is exactly the read length from your Sample Sheet plus that of your complete barcode sequence. All entered indices must be of the same length. To provide multiple indices, enter the -XXX argument for every barcode or barcode combination, e.g.:

	hilive /path/to/BaseCalls /path/to/index.kix 107 /path/to/outputFolder -b barcode1 -b barcode2 ...

One output file will be produced for each barcode. To get alignments with undetermined barcodes, activate the --keep-all-barcodes parameter.
Dual barcodes must be delimited with "-" (e.g., -b ATCGTGAT-TAGTTAGC for a 2x8bp barcode).

---

Please consult the project website for more details on the parameters!


License
-------

See the file LICENSE for licensing information.


Contact
-------

Please consult the HiLive project website for questions!

If this does not help, please feel free to consult:

 * Technical support <hilive.team@gmail.com> (technical contact)
 * Bernhard Y. Renard <renardb (at) rki.de> (project head)

also see CONTRIBUTORS for a complete list of contributors and their contact information
