<!-- {#mainpage} -->

HiLive2 - Live Mapping of Illumina reads
========================================

- [Description](#description)
- [Website](#website)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)

Description
-----------

HiLive2 is a read mapping tool that maps Illumina HiSeq (or comparable) 
reads right in the moment when they are produced. This means, read mapping 
is finished as soon as the sequencer is finished. Alignment output can also 
be written for intermediate cycles providing insights into the data already 
during sequencing. HiLive2 supports paired-end reads and live demultiplexing, 
creating one output alignment file for each barcode and cycle.


Website
-------

The general HiLive project website is 
https://gitlab.com/rki_bioinformatics/HiLive.  
Here, you find information and links to HiLive and HiLive2.  

The direct link to HiLive2 is https://gitlab.com/LokaT/HiLive2.  
There you can find the latest version of HiLive2, source code, documentation and
example data.  

Please report bugs and problems by 
[opening an issue](https://gitlab.com/LokaT/HiLive2/issues) or mail us to 
hilive.team (at) gmail.com.


Installation
------------

HiLive2 can currently only be installed from source.  

The following dependencies are required:

 * cmake (>= 3.10)
 * boost (system, filesystem, program\_options, iostreams)
 * zlib
 * lz4
 * [SeqAn](http://packages.seqan.de/) (version 2.3.2)

---

Check out the HiLive2 source code from the project website and `make` the 
project.  
We recommend to have separated `source` and `build` directories:

```
mkdir HiLive2 && cd HiLive2
git clone https://gitlab.com/LokaT/HiLive2 source
mkdir build && cd build
cmake ../source -DSEQAN_PATH="PATH/TO/UNPACKED/SEQAN/"
make
```

If using local versions of boost or lz4, please specify the path with the `cmake` command using `-DLZ4_PATH="/PATH/TO/LZ4"` and `-DBOOST_PATH="/PATH/TO/BOOST"`, respectively.  


Usage
-----

HiLive2 has three components:

 * ``hilive-build``  builds the index of a reference genome
 * ``hilive``        the read mapper itself
 * ``hilive-out``    write output files of existing HiLive data

#### Using hilive-build:

Building a HiLive2 index from FASTA file input.fa:

```
hilive-build --input input.fa --out-prefix /path/to/index/directory/prefix
```

The command will produce several files with the declared prefix.  

#### Using hilive:

To map reads in a 100bp run (without barcodes) using default settings:

    hilive --bcl-dir /path/to/BaseCalls --index /path/to/prefix --reads 100R

For an overview of additional parameters, type 
	
	hilive --help

#### Using hilive-out:

hilive-out required a config file as input.  
This should be the config file that is automatically written in the temporary 
directory of the original HiLive2 run.  
To write output from existing temporary files in HiLive, type:

```
hilive-out --config /path/to/temp/hilive_config.ini
```

This command will report the alignment results of the last cycle based on the 
settings that were specified for the related HiLive2 run. To produce output 
files for other cycles, e.g. 50, 70, 90, type:

```
hilive-out --config /path/to/temp/hilive_config.ini --out-cycles 50,70,90
```

#### Tutorial
We provide more extensive explanations and examples in the [Tutorial](https://gitlab.com/lokat/hilive2/tutorial).


License
-------

See the file LICENSE for licensing information.


Contact
-------

Please consult the HiLive project website for questions!

If this does not help, please feel free to consult:
 * Technical support: hilive.team (at) gmail.com (technical contact)
 * Bernhard Y. Renard: renardb (at) rki.de (project head)

also see CONTRIBUTORS for a complete list of contributors and their contact 
information.
