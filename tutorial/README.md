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
  - [Important general options](#important-general-options)
  - [Paired-end sequencing](#paired-end-sequencing) 
  - [Live demultiplexing](#live-demultiplexing)
  - [Output cycles](#output-cycles)
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
> /path/to/hilive -h
> /path/to/hilive-out -h
> /path/to/hilive-build -h
```

### Example data
In the provided sample data, you find the genome of Human immunodeficiency virus 1 (HIV1) in FASTA format. The genome was obtained from NCBI under accession number NC_001802.1 (https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1).  
From this genome, we simulated 100 paired-end reads of length 2x100bp with the mason read simulator (Holtgrewe, 2010). We manually added barcode sequences of length 8bp between both reads and converted them to Illumina BaseCalls format using our `fastq2bcl.py` script.

Get started
-----------

At the beginning of the tutorial, please change your working directory to the tutorial data:

```
> cd /path/to/HiLive2/tutorial
> ls
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
> mkdir index
> hilive-build -i hiv1.fa -o index/hiv1
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
> hilive -b BaseCalls -i index/hiv1 -r 100R --lanes 1 --tiles 1101 
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


HiLive2 options
---------------

### Important general options:

This sections lists the options of HiLive2 that we expect to be most important when using HiLive2.  
This short introductions does not contain information about options concerning [paired-end sequencing](#paired-end-sequencing), [live demultiplexing](#live-demultiplexing) and [output cycles](#output-cycles) as they are described in separated section to provide more details.  

##### Alignment mode
HiLive2 comes with three different alignment modes: `FAST`, `BALANCED` or `ACCURATE`.  
By default, HiLive2 is started in `BALANCED` mode. Use the `--align-mode` option to make HiLive2 run faster (`FAST`) or more sensitive (`ACCURATE`):

```
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R --align-mode FAST
```

The selected alignment mode can affect several algorithmic parameters, as the anchor length or seeding interval.

##### Report mode
HiLive2 comes with several report modes:

 * `ANYBEST` reports one best alignment  
 * `ALLBEST` reports all best alignments
 * `BESTN#` reports the best up to # alignments (also suboptimal)
 * `ALL` reports all alignments found by the algorithm
 * `UNIQUE` reports only unique alignments (only aligning to one position)

Use the `--out-mode` option to select the respective report mode:

```
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R --out-mode BESTN5
```

##### Lanes and tiles
Different Illumina sequencing machines (and run modes/protocols) have a different number of lanes and tiles.  
These can be set with the `--lanes` and `--tiles` option as a comma-separated list:

```
# This command will not work with the example data since lane 2 and tile 1102 are missing.
> hilive -b BaseCalls -i index/hiv1 --lanes 1,2 --tiles 1101,1102 -r 100R
```

Since the number of tiles can be very large for real sequencing applications (e.g., 96 tiles for Illumina HiSeq in High Output mode), HiLive2 also comes with the `--max-tile` option.
By typing the maximum tile number that occurs for your sequencing run, all other numbers will be computed automatically according to Illumina nomenclature.  
The following command will include the tiles 1101-1116, 1201-1216, 2101-2116 and 2201-2216:

```
# This command will not work with the example data since several tiles are missing.
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --max-tile 2216 -r 100R
```

Alternatively, the number of lanes and tiles can be determined automatically from Illumina's `runInfo.xml` file by using the `--runinfo` option of HiLive2 (not shown in this tutorial).

##### Output and temp directories
HiLive2 writes output and temp files in two different directories.  
By default, this is `./out` and `./temp` for the output and temporary files, respectively. The directories are created in your current working location.  
Use the `-o [--out-dir]` and `--temp-dir` parameters to change the output and temporary directory, respectively. Not existing directories are automatically created.  
  
NOTE: Existing files from previous runs will be overwritten. We strongly recommend to use separated output and temporary directories for each sequencing run!

##### Output format
You may change the output format from `BAM` (default) to `SAM` by using `--out-format SAM` or `-f SAM`.

##### Multithreading
Please specify the maximal number of threads with the `-n [--num-threads]` option.  
The default number of threads is 1. However, we recommend to use one thread for each analyzed tile.


### Paired-end sequencing

To align both reads rather than only the first 100bp, the segment order must be specified with the `-r [--reads]` option.  
In Illumina sequencing, barcodes of paired-end sequencing mostly occur in the middle of both reads. For the given data, we simulated 2x100bp reads and 2x4bp barcodes in between them.
This corresponds to the segment structure 100R-4B-4B-100R. In HiLive2, specify the segments order as a comma-separated list:

```
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R,4B,4B,100R 
```

Alignments for both reads are written to the same BAM file as defined in the SAM/BAM specification.  
However, HiLive2 only reports independent alignments for both reads rather than proper pairs.


### Live demultiplexing

With live demultiplexing, it is possible to distinguish different samples of the same sequencing run.  
HiLive2 supports arbitrary combinations of the segment structure which also allows to use dual barcodes as shown in this example (2x4bp).  
Before you continue, please remove all previously created data:

```
> rm -r ./out
> rm -r ./temp
```

Our example data set contains reads with two slightly different Barcodes:
```
ACAG-TCGA
|| |-| ||
ACGG-TGGA
```

With default parameters, two errors per barcode fragment are tolerated. Thus, all reads fulfill the criteria and will be aligned and reported.  
The desired barcode sequence(s) are declared with the `-B [--barcodes]` option.  
The number of permitted barcodes per fragment can be set with `--barcode-errors`.  
While both parameters are comma-separated lists in general, the tolerated number of errors can also be set with a single number to set this value for all fragments. 
Dual barcodes are separated with a "-" character. 
To also align and report undetermined barcodes, set the `--align-undetermined-barcodes` option.  
There are 75 reads with the first and 25 reads with the second barcode. 
For each read, there are two reported alignments (one per segment).  
The following examples show how the number of reads in the respective output files behaves for the given data: 

```
# Default barcode options
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R,4B,4B,100R -B ACAG-TCGA --align-undetermined-barcodes
[...]
> samtools view out/hilive_out_cycle208_ACAG-TCGA.bam | wc -l
200
> samtools view out/hilive_out_cycle208_undetermined.bam | wc -l
0

# Strict barcode options (not tolerating errors)
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R,4B,4B,100R -B ACAG-TCGA --barcode-errors 0 --align-undetermined-barcodes
[...]
> samtools view out/hilive_out_cycle208_ACAG-TCGA.bam | wc -l
150
> samtools view out/hilive_out_cycle208_undetermined.bam | wc -l
50

# Tolerate 0 errors in the first barcode segment and 1 error in the second, don't align reads with undetermined barcode
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R,4B,4B,100R -B ACAG-TCGA --barcode-errors 0,1
[...]
> samtools view out/hilive_out_cycle208_ACAG-TCGA.bam | wc -l
150
> samtools view out/hilive_out_cycle208_undetermined.bam | wc -l
[E::hts_open_format] fail to open file 'out/hilive_out_cycle208_undetermined.bam'
samtools view: failed to open "out/hilive_out_cycle208_undetermined.bam" for reading: No such file or directory
0

# Specify both barcodes, don't align reads with undetermined barcode. 
# In this example, all reads match the first barcode with the default number of tolerated error, thus all reads are in the output file of ACAG-TCGA
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R,4B,4B,100R -B ACAG-TCGA,ACGG-TGGA
[...]
> samtools view out/hilive_out_cycle208_ACAG-TCGA.bam | wc -l
200
> samtools view out/hilive_out_cycle208_ACGG-TGGA.bam | wc -l
0

# If the tolerated number of error is adapted, the barcodes are assigned as expected:
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R,4B,4B,100R -B ACAG-TCGA,ACGG-TGGA --barcode-errors 0
[...]
> samtools view out/hilive_out_cycle208_ACAG-TCGA.bam | wc -l
150
> samtools view out/hilive_out_cycle208_ACGG-TGGA.bam | wc -l
50

```


### Output cycles
When performing real-time read mapping, it is often valuable to get intermediate mapping results when the sequencing machine is still running rather than only producing results after the last sequencing cycle.  
Therefore, it is possible to write SAM/BAM output files for intermediate sequencing cycles using the `-O [--out-cycles]` option.  
In the following example, output files are written in the cycles 50, 100, 158 and 208:

```
> hilive -b BaseCalls -i index/hiv1 --lanes 1 --tiles 1101 -r 100R,4B,4B,100R -B ACAG-TCGA -O 50,100,158,208
[...]
> ls -lah ./out
hilive_out_cycle50_ACAG-TCGA.bam
hilive_out_cycle100_ACAG-TCGA.bam
hilive_out_cycle158_ACAG-TCGA.bam
hilive_out_cycle208_ACAG-TCGA.bam
```

