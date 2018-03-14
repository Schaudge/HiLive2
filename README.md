HiLive2 - Live Mapping of Illumina reads
========================================


Content
-------

 - [Description](#description)
 - [Website](#website)
 - [Installation](#installation)
 - [Usage](#usage)
  - [hilive-build](#using-hilive-build)
  - [hilive](#using-hilive)
  - [hilive-out](#using-hilive-out)
  - [Demultiplexing](#demultiplexing)
  - [Use a config file](#config-files)
  - [List of parameters](#list-of-parameters)
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

 * cmake (>= 2.8)
 * boost (system, filesystem, program\_options)
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
cmake ../source
make
```

If using a local version of lz4 you will need to adjust path in CMakeLists.txt 
(line 32).  
Please also adjust the paths to the SeqAn library in CMakeLists.txt (lines 41 
and 42).


Usage
-----

HiLive2 has three components:

 * ``hilive-build``  builds the index of a reference genome
 * ``hilive``        the read mapper itself
 * ``hilive-out``    write output files of existing HiLive data

---

#### Using hilive-build:

Building a HiLive2 index from FASTA file input.fa:

```
hilive-build --input input.fa --out-prefix /path/to/index/directory/prefix
```

The command will produce several files with the declared prefix.  

---

#### Using hilive:

To map reads in a 100bp run (without barcodes) using default settings:

    hilive --bcl-dir /path/to/BaseCalls --index /path/to/prefix --reads 100R

For an overview of additional parameters, type 
	
	hilive --help

---

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

Please note, that the temporary files for the respective cycles must be present 
in the temp folder. This is only the case if the --keep-files parameter and/or 
the --out-cycles parameter for the respective cycles was activated for the 
original HiLive2 run.

---

#### Demultiplexing:

To map reads from multiplexed sequencing runs, you can provide HiLive2 with the 
barcode sequences from your Sample Sheet.  
In default cases, barcode sequences are read after the (first) read, such that 
demultiplexing is carried out after the mapping is completed.  
If you use double indexing, please concatenate both indices in the correct 
order and provide them as one sequence. Separate different barcodes by "," and 
dual barcodes by "-":

```
hilive --bcl-dir /path/to/BaseCalls --index /path/to/prefix --reads 101R,8B,8B,101R --barcodes ACGTACGT-TGCATGCA,GCTAGCTA-TGTGCACA
```

---

#### Config files:

HiLive2 comes with a `--config` option to load parameters from a configuration 
file:

```
hilive --config config.ini
```

All options can be defined in a config file instead of using command line 
parameters. You can also combine parameters from a config file and the command 
line:

```
hilive --config config.ini --bcl-dir ./BaseCalls
```

Thus, it is possible to use configuration files as a template and specify the 
options of a single run with command line parameters. With this concept in 
mind, please note that command line input has priority over options obtained 
from the configuration file.  
  
The configuration file must be given in ".ini" file format.  
In this format, each parameter can only be specified once. Thus, please 
separate multitoken arguments with the "," character. This is the same syntax 
for multitokens as required if you specify arguments via the command line.  
  
##### Example: 

config_hg19_miseq.ini:
```
index=/home/user/index/hilive2/hg19
reads=101R,8B,8B,101R
out-format=BAM
num-threads=38
lanes=1
max-tile=2119
align-mode=fast
keep-all-files=true
```

Run HiLive2:
```
hilive --config config_hg19_miseq.ini --bcl-dir /path/to/BaseCalls 
```
Please note that, in this example, `--bcl-dir` must be specified when running 
HiLive2 since it is required but not set in the configuration file.  

---

#### List of parameters

| *REQUIRED OPTIONS* | |
| ------------------ | --- |
| **-b [--bcl-dir] ARG** | Illumina's BaseCalls directory which contains the sequence information of the reads. |
| **-i [--index] ARG** | Path to the HiLive index. |
| **-r [--reads] ARG[,ARG,..]** | Length and types of the read segments. |

---

| *GENERAL OPTIONS* | |
| ------------------ | --- |
| **-h [--help]** | Print the help message and exit. |
| **-l [--license]** | Print the license information and exit. |
| **-c [--config] ARG** | Path to a configuration file in ".ini" format. |
| **--runinfo ARG** | Path to Illuminas `runInfo.xml` file. |
| **--continue ARG** | Continue an interrupted HiLive run from a specified cycle. |
  
---
  
| *SEQUENCING OPTIONS* | |
| ------------------ | --- |
| **-b [--bcl-dir] ARG** | Illumina's BaseCalls directory which contains the sequence information of the reads. [REQUIRED]|
| **-l [--lanes] ARG[,ARG,..]** | Specify the lanes to be considered for read alignment. [Default: 1-8] |
| **-t [--tiles] ARG[,ARG,..]** | Specify the tiles to be considered for read alignment. [Default:[1-2][1-3][01-16]; 96 tiles] |
| **-T [--max-tile] ARG** | Specify the highest tile number. [Default: 2316] |
| **-r [--reads] ARG[,ARG,..]** | Length and types of the read segments. [REQUIRED]|
| **-B [--barcodes] ARG[,ARG,..]** | Barcode(s) of the sample(s) to be considered for read alignment. |
| **--run-id ARG** | ID of the sequencing run.  |
| **--flowcell-id ARG** | ID of the flowcell. |
| **--instrument-id ARG** | ID of the sequencing machine. |

---

| *REPORT OPTIONS* | |
| ------------------ | --- |
| **-o [--out-dir] ARG** | Path to the directory that is used for the output files. [Default: ./out] |
| **-f [--out-format] ARG** | Format of the output file(s). [Default: BAM] |
| **-O [--out-cycles] ARG[,ARG,..]** |  Cycles for that alignment output is written. |
| **-M [--out-mode] ARG** | The output mode. [Default: ANYBEST] |
| **--report-unmapped** |  Activate reporting unmapped reads. [Default: false] |
| **--extended-cigar** | Activate extended CIGAR format for the alignment output files. [Default: false] |
| **--force-resort** | Always sort temporary alignment files before writing output. [Default: false] |
| **--max-softclip-ratio ARG** | Maximal relative length of the front softclip during output. [Default: 0.2] |

---

| *ALIGNMENT OPTIONS* | |
| ------------------ | --- |
| **-i [--index] ARG** | Path to the HiLive index. [REQUIRED] |
| **-m [--align-mode] ARG** |  Alignment mode to balance speed and accuracy. [Default: balanced] |
| **--anchor-length ARG** |  Length of the alignment anchor (or initial seed) [Default: set by the selected alignment mode] |
| **--error-interval ARG** |  The interval to tolerate more errors during alignment. [Default: 'anchor-length'/2] |
| **--seeding-interval ARG** |  The interval to create new seeds. [Default: 'anchor-length'/2] |
| **--barcode-errors ARG[,ARG,..]** |  The number of errors that are tolerated for the barcode segments. [Default: 1] |
| **--align-undetermined-barcodes** |  Align reads with undetermined barcodes. [Default: false] |
| **--min-basecall-quality ARG** | Minimum basecall quality for a nucleotide to be considered as a match. [Default: 1] |
| **--keep-invalid-sequences** | Keep sequences of invalid reads. [Default: false] |

---

| *SCORING OPTIONS* | |
| ------------------ | --- |
| **-s [--min-as] ARG** | Minimum alignment score. [Default: Depends on the alignment mode] |
| **--match-score** | Score for a match. [Default: 0] |
| **--mismatch-penalty ARG** | Penalty for a mismatch. [Default: 6] |
| **--insertion-opening-penalty ARG** | Penalty for insertion opening. [Default: 5] |
| **--insertion-extension-penalty ARG** | Penalty for insertion extension. [Default: 3] |
| **--deletion-opening-penalty ARG** | Penalty for deletion opening. [Default: 5] |
| **--deletion-extension-penalty ARG** | Penalty for deletion extension. [Default: 3] |
| **--max-gap-length ARG** | Maximal permitted consecutive gap length. [Default: 3] |
| **--softclip-opening-penalty ARG** | Penalty for softclip opening during output. [Default: 'mismatch-penalty]' |
| **--softclip-extension-penalty ARG** |  Penalty for softclip extension (only relevant during output). [Default: 'mismatch-penalty'/'anchor-length'] |

---

| *TECHNICAL OPTIONS* | |
| ------------------ | --- |
| **--temp-dir ARG** | Temporary directory to store the alignment files and hilive_config.ini. [Default: ./temp] |
| **-k [--keep-files] ARG[,ARG,..]** | Keep intermediate alignment files for these cycles. [Default: Keep files of output cycles] |
| **-K [--keep-all-files]** |  Keep all intermediate alignment files. [Default: false] |
| **--block-size ARG** | Block size for the alignment input/output stream in Bytes. [Default: 64M] |
| **--compression ARG** |  Compression of temporary alignment files. [Default: LZ4] |
| **-n [--num-threads] ARG** |  Number of threads to spawn (including output threads). [Default: 1] |
| **-N [--num-out-threads] ARG** | Maximum number of threads to use for output. [Default: 'num-threads'/2] |


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
