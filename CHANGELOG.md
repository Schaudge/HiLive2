HiLive - Live Mapping of Illumina reads
=======================================

Changelog
-----------

v1.1
-----

General functionality:
 * Continue a hilive run in a specified cycle if the temporary file is available (--continue ARG)
 * Keep temporary alignment files of specified cycles (--keep-files ARG)
 * New option to keep all temporary alignment files (--keep-all-files)

Output:
 * Multithreaded SAM/BAM output
 * Different alignment modes are supported for intermediate sequencing cycles
 * New output mode UNIQUE to report only unique alignments
 * Output mode can be changed when using hilive-out executable without losing information
 * Option to specify a minimum alignment score ratio for output (--min-as-ratio).

v1.0
-----

General functionality:
 * Support of gapped k-mers
 * Support of arbitrary sequence structures (e.g. for paired-end sequencing)
 * Support of dual barcodes
 * Extensive algorithmic optimizations

Usability:
 * K-mers are defined for each index individually instead of defining at compile time
 * K-mer definition is loaded from the index when running HiLive
 * Read sequencing information from Illuminas RunInfo.xml
 * Load HiLive settings from a settings file
  
Index:
 * Slightly changed index. Index files of versions <1.0 are no longer supported.
 
 Output:
 * Completely revised output structure (one file per barcode; no separation in lanes or tiles)
 * Support BAM output
 * Real-time SAM/BAM output for intermediate sequencing cycles
 * Optional: Extended CIGAR format
 * Executable "hilive-out" to create SAM/BAM output from existing temporary alignment files

v0.3
-----

 * Live demultiplexing
 * Code optimization

v0.2
-----

 * Code optimization
 
v0.1
-----
 
 * Initial version  