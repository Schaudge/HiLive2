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
is finished as soon as the sequencer is finished.  
Alignment output can also be written for intermediate cycles providing 
insights into the data already during sequencing.  


Website
-------

The general HiLive project website is https://gitlab.com/rki_bioinformatics/HiLive  
Here, you find information and links to HiLive and HiLive2.  

The direct webpage of HiLive2 is https://gitlab.com/LokaT/HiLive2  
There you can find the latest version of HiLive2, source code, documentation and
example data.  

Please report bugs and problems as an issue at https://gitlab.com/LokaT/HiLive2


Installation
------------

HiLive2 can currently only be installed from source.  

The following dependencies are required:

 * cmake (>= 2.8)
 * boost (system, filesystem, program\_options)
 * zlib
 * lz4
 * [SeqAn](http://packages.seqan.de/) (version 2.3.2)

If using a local version of lz4 then adjust path in CMakeLists.txt line 32.  
Please also adjust the paths of the SeqAn library in CMakeLists.txt lines 41 and 42.  

---

Check out the HiLive source code from the project website and make the project.  
We recommend to have separated `source` and `build` directories:

```
mkdir HiLive2 && cd HiLive2
git clone https://gitlab.com/LokaT/HiLive2 source
mkdir build && cd build
cmake ../source
make
```

If using a local version of lz4 then adjust path in CMakeLists.txt line 32.
Please also adjust the paths of the SeqAn library in CMakeLists.txt lines 41 and 42.

Usage
-----

HiLive has three components:

 * ``hilive-build``  builds the index of a reference genome
 * ``hilive``        the read mapper itself
 * ``hilive-out``    write output files of existing HiLive data

#### Using hilive-build:

Building a HiLive index from FASTA file input.fa:

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
This should be the config file that is automatically written in the temporary directory of the original HiLive2 run.  
To write output from existing temporary files in HiLive, type:

	hilive-out --config /path/to/temp/hilive_config.ini

This command will report the alignment results of the last cycle based on the settings that were specified for the related HiLive run.
To produce output files for other cycles, e.g. 50, 70, 90, type:

	hilive-out --config /path/to/temp/dir/hilive_settings.xml --out-cycles 50 70 90
	
Please note, that the temporary files for the respective cycles must be present in the temp folder.
This is only the case if the --keep-files parameter and/or the --out-cycles parameter for the respective cycles was activated in the HiLive run.

#### Tutorial
We provide more extensive explanations and examples in the [Tutorial](https://gitlab.com/lokat/hilive2/tutorial).


License
-------

See the file LICENSE for licensing information.


Contact
-------

Please consult the HiLive project website for questions!

If this does not help, please feel free to consult:
 * Technical support <hilive.team (at) gmail.com> (technical contact)
 * Bernhard Y. Renard <renardb (at) rki.de> (project head)

also see CONTRIBUTORS for a complete list of contributors and their contact information
