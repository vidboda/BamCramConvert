# BamCramConvert
Script to convert bam2cram and reverse using find, samtools (requires UNIX find and samtools...), and optionally a slighlty modified [forked](https://github.com/beboche/bam2cram-check) version of [bam2cram-check](https://github.com/wtsi-hgi/bam2cram-check). BamCramConvert can also apply [Crumble](https://github.com/jkbonfield/crumble) compression to your files!

Usage : bash bcc.sh

	Mandatory arguments :
  
		* -d|--directory	<path to search dir>: root dir for find command    
		* -s|--size		<File size to search (man find -size)>: ex: +200000000k or +20G will search for files greater than 20Go; see man find -size argument    
		* -mt|--modif-time	<File last modif to search (man find -mtime): ex: +180 will search for files older than 6 months; see man find -mtime argument    
		* -f|--file-type	<bam|cram>: file type to find and convert from (bam will search for bam files and convert to cram)
    
	Optional arguments :
  
		* -rm|--remove		:removes original file and index (in case of full conversion success implies bam2cram-check) - default: false
		* -st|--samtools	<path to samtools> - default: try to locate in PATH
		* -fa|--ref-fasta	<path to ref genome .fa>: path to a fasta file reference genome (the directory containing the fasta file must also contain samtools index) - default: /usr/local/share/refData/genome/hg19/hg19.fa
		* -c|--check		: uses bam2cram-check (slightly modified) to check the conversion - implicitely included with -rm - if fails and -rm: rm canceled) - requires python >3.5 and samtools > 1.3
		* -p|--python3		<path to python3> - used in combination with -c or -rm: needed to run submodule bam2cram-check - default: /usr/bin/python3 - python version must be > 3.5
		* -uc|--use-crumble     : uses crumble to compress the converted BAM/CRAM file - Note: a file that already contains "crumble" in its name will not be converted again
		* -cp|--crumble-path    <path to crumble> - used in combination with -uc: needed to run crumble - default: try to locate in PATH
		* -v | --verbosity 	<integer> : decrease or increase verbosity level (ERROR : 1 | WARNING : 2 | INFO : 3 | COMMAND [default] : 4 | DEBUG : 5)
    
	General arguments :
  
		* -sl|--slurm   : when running in SLURM environnment, generates srun commands - default: false
		* -th|--threads : number of threads to be used for samtools -@ option (0 => 1 total thread, 1 => 2 total threads...)
		* -h		: show this help message and exit
		* -t		: test mode (dont execute command just print them) - default: false

Installation:

see INSTALL.md

Use cases:

* convert BAM files in directory dir (recursively) into CRAM, BAM size min 100Mo and older than 7 days, do not remove original BAM files, Reference genome is hg38 (must match with BAM files to be converted, non matching files based on chr1 length will be ignored - supports hg19 and hg38 only):

```bash
bash bcc.sh -d path/to/dir/ -mt +7 -s +100M -f bam -fa /path/to/hg38.fa
```

* the same in SLURM environnment using 4 threads for samtools commands:

```bash
bash bcc.sh -d path/to/dir/ -mt +7 -s +100M -f bam -fa /path/to/hg38.fa -sl -th 3
```

* the same in dry run mode and providing optional samtools path (otherwise is searched in PATH) - it is highly recommanded to use the dry run mode before launching a command, just to be sure of what will be done:

```bash
bash bcc.sh -d path/to/dir/ -mt +7 -s +100M -f bam -fa /path/to/hg38.fa -sl -th 3 -st /special/place/samtools -t
```

* convert CRAM files in directory dir (recursively) into BAM, CRAM size min 5Go, from today, remove original CRAM files (implies bam2cram-check):

```bash
bash bcc.sh -d patho/to/dir/ -mt -1 -s +5G -f cram -fa /path/to/hg38.fa -rm
```

* the same as above without removing the original files and now we explicitely apply bam2cram-check to the new files (and provide optional python3 (>3.5) path):

```bash
bash bcc.sh -d patho/to/dir/ -mt -1 -s +5G -f cram -fa /path/to/hg38.fa -c -p /usr/bin/python3
```

* convert BAM files in directory dir (recursively) into CRAM, BAM size min 100Mo and older than 7 days, apply crumble to the new CRAM and provide optional crumble path (otherwise is searched in PATH), remove original file:

```bash
bash bcc.sh -d path/to/dir/ -mt +7 -s +100M -f bam -fa /path/to/hg38.fa -uc -cp /special/place/crumble -rm
```




