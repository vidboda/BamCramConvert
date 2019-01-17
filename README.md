# BamCramConvert
script to convert bam2cram and reverse using find and samtools (requires UNIX find and samtools...)

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
    
	General arguments :
  
		* -sl|--slurm    : when running in SLURM environnment, generates srun commands - default: false
		* -h		: show this help message and exit
		* -t		: test mode (dont execute command just print them) - default: false
