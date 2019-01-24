# BamCramConvert
Script to convert bam2cram and reverse using find, samtools (requires UNIX find and samtools...), and optionally a slighlty modified [forked](https://github.com/beboche/bam2cram-check) version of [bam2cram-check](https://github.com/wtsi-hgi/bam2cram-check).

Usage : bash bcc.sh

	Mandatory arguments :
  
		* -d|--directory	<path to search dir>: root dir for find command    
		* -s|--size		<File size to search (man find -size)>: ex: +200000000k will search for files greater than 20Go; see man find -size argument    
		* -mt|--modif-time	<File last modif to search (man find -mtime): ex: +180 will serach for files older than 6 months; see man find -mtime argument    
		* -f|--file-type	<bam|cram>: file type to find and convert from (bam will search for bam files and convert to cram)
    
	Optional arguments :
  
		* -rm|--remove		:removes original file and index (in case of full conversion success) - default: false
		* -st|--samtools	<path to samtools> - default: /usr/local/bin/samtools
		* -fa|--ref-fastq	<path to ref genome .fa>: path to a fasta file reference genome (the directory containing the fasta file must also contain samtools index) - default:/usr/local/share/refData/genome/hg19/hg19.fa
		* -c|--check		: uses bam2cram-check (slightly modified) to check the conversion - implicitely included with -rm - if fails and -rm: rm canceled) - requires python >3.5 and samtools > 1.3
		* -p|--python3		<path to python3> - used in combination with -c or -rm: needed to run submodule bam2cram-check - default: /usr/bin/python3
    
	General arguments :
  
		* -sl|--slurm   : when running in SLURM environnment, generates srun commands - default: false
		* -th|--threads : number of threads to be used for samtools -@ option (0 => 1 total thread, 1 => 2 total threads...)
		* -h		: show this help message and exit
		* -t		: test mode (dont execute command just print them) - default: false
