HiLive2 - Live Mapping of Illumina reads (Tutorial)
===================================================

Content
-------

 - [Introduction](#introduction)
  - [Installation](#installation)
  - [Example data](#example-data)
 - [Get started](#get-started)
  - [Index building](#index-building)
  - [Run HiLive2](#run-hilive2)
 - [HiLive2 options](#hilive2-options)
 - [Using a configuration file](#using-a-configuration-file)
 - [Using the hilive-out executable](#using-the-hilive-out-executable)
 - [Additional remarks](#additional-remarks)


Introduction
------------

HiLive2 is a real-time read mapping software. Input data are directly obtained from an Illumina sequencing machine. Thus, your machine must have access to the output directory of the sequencer.  
We provide several scripts for creating test data for HiLive2, e.g. `fastq2bcl.py` to convert Illumina reads in FASTQ format to the Illumina BaseCalls format. Please find these scripts on the [HiLive2 repository](https://gitlab.com/lokat/hilive2/scripts).  

### Installation
Please make sure you installed HiLive2 as described in the README.md of the program.  
To test HiLive2, please try to open the help of the three executables:

```
/path/to/hilive -h
/path/to/hilive-out -h
/path/to/hilive-build -h
```

### Example data
In the provided sample data, you find the genome of Human immunodeficiency virus 1 (HIV1) in FASTA format. The genome was obtained from NCBI under accession number NC_001802.1 (https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1).  
From this genome, we simulated 100 paired-end reads of length 2x100bp with the mason read simulator (Holtgrewe, 2010). We manually added barcode sequences of length 8bp between both reads and converted them to Illumina BaseCalls format using our `fastq2bcl.py` script.

Get started
-----------

At the beginning of the tutorial, please change your working directory to the tutorial data:

```
cd /path/to/HiLive2/tutorial
ls
```

The `ls` command should list four files/directories:

```
BaseCalls
config.ini
hiv1.fa
README.md
```

### Index building
To use HiLive2, you first need to build an index of your reference genome(s) of interest. In our example, this is the reference genome of HIV1.  
Use the `hilive-build` executable to build your index. Please note, that the output directory must exist:

```
mkdir index
hilive-build -i hiv1.fa -o index/hiv1
```

This should only take a few seconds and the command line output should look like this:

```
________________________________________________________________________________

HiLive Index Builder v2.0 - Build Index for Realtime Alignment of Illumina Reads
________________________________________________________________________________
                                                                                
Create FM-index from file hiv1.fa ...
Index was successfully built to file(s) index/hiv1                              
```

### Run HiLive2
After [buiding your index](#index-building), you can run HiLive2 with the prepared input data. The data are in Illumina BaseCall format and placed in the `BaseCalls` directory. The structure of the data in our example is the same as for real data. Depending on the actual sequencing machine (HiSeq, MiSeq, etc.), the data is split into different lanes and tiles. In our example, there is only one lane (L001) and one tile (1101).  
To run HiLive2 with default parameters, type:

```
hilive -b BaseCalls -i index/hiv1 -r 100R --lanes 1 --tiles 1101 
```

The command line output should look similar to this:

```
__________________________________________________

HiLive v2.0 - Realtime Alignment of Illumina Reads
__________________________________________________

Running HiLive with       1 thread(s).
BaseCalls directory:      BaseCalls
Temporary directory:      ./temp
BAM output directory:     ./out
Lanes:                    1
K-mer index:              index/hiv1
Read lengths:             100R
Min. alignment score:     -24
Mapping mode:             ANYBEST
Anchor length:            11

Loading Index ...
Waiting for the first cycle to finish...
Initializing Alignment files...
First cycle complete. Starting alignment.
Creating 1 threads.
Task [Lane 1 Tile 1101 Cycle 1.1]: Found 0 seeds.
Task [Lane 1 Tile 1101 Cycle 1.2]: Found 0 seeds.
Task [Lane 1 Tile 1101 Cycle 1.3]: Found 0 seeds.

[...]

Task [Lane 1 Tile 1101 Cycle 1.100]: Found 103 seeds.
Finished all alignments.
Waiting for output tasks...
Finished output tasks.
All threads joined.
Total run time: 31 s
``` 

HiLive2 created two new directories: `out` and `temp`. In the `temp` directory, you find temporary alignment files of the last cycle (100). Additionally, there is a file `hilive_config.ini` which stores the parameter settings used by HiLive2 for the given run. The `out` directory contains the actual output files in BAM format. There is one cycle for each cycle and barcode. For the example command with default parameters, this is only one single file `hilive_out_cycle100_undetermined.bam`.
