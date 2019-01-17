#!/bin/bash

###########################################################################
#########							###########
######### bam2cram.sh						###########
######### @uthor : D Baux	david.baux<at>inserm.fr		###########
######### Date : 17/01/2019					###########
#########							###########
###########################################################################

###########################################################################
###########
########### 	Script to automate BAM/CRAM conversions
###########
###########################################################################

VERSION=1.0
USAGE="
Program: bam2cram
Version: ${VERSION}
Contact: Baux David <david.baux@inserm.fr>
Usage: bash bam2cram.sh -d <path to search dir> -s <File size to search (man find -size)> -mt <File last modif to search (man find -mtime)> -f <bam|cram> -st <path to samtools> -fa <path to ref genome .fa>
	"

if [ $# -eq 0 ]; then
	echo "${USAGE}"
	echo "Error Message : Arguments not provided"
	echo ""
	exit 1
fi

usage ()
{
	echo 'This script convert bam files to cram  and reverse.';
	echo 'Usage : bash bam2cram.sh';
	echo '	Mandatory arguments :';
	echo '		* -d|--directory	<path to search dir>: root dir for find command';
	echo '		* -s|--size		<File size to search (man find -size)>: ex: +200000000k will search for files greater than 20Go; see man find -size argument';
	echo '		* -mt|--modif-time	<File last modif to search (man find -mtime): ex: +180 will serach for files older than 6 months; see man find -mtime argument';
	echo '		* -f|--file-type	<bam|cram>: file type to find and convert from (bam will search for bam files and convert to cram';
	echo '	Optional arguments :';
	echo '		* -rm|--remove		:rm converted file';
	echo '		* -st|--samtools	<path to samtools>: default: /usr/local/bin/samtools';
	echo '		* -fa|--ref-fastq	<path to ref genome .fa>: path to a fasta fils ref genome (the directory containing the fasta file must also contain samtools index)';
	echo '	General arguments :';
	echo '		* -sl|--slurm    : when running in SLURM environnment, generates srun commands - cannot be used in combination with -c 0';
	echo '		* -h		: show this help message and exit';
	echo '		* -t		: test mode (dont execute command just print them)';
	echo '';
	exit
}

# --Option processing
SLURM_MODE=false
TEST_MODE=false
RM=false
SAMTOOLS=/usr/local/bin/samtools
REF_FASTA=/usr/local/share/refData/genome/hg19/hg19.fa

# --Parse command line
while [ "$1" != "" ];do
	case $1 in
		-d | --directory )	shift
			DIR=$1
			;;
		-s | --size )	shift
			SIZE=$1
			;;
		-mt | --modif-time )	shift
			MTIME=$1
			;;
		-f | --file-type )	shift
			FILE_TYPE=$1
			;;
		-st | --samtools )	shift
			SAMTOOLS=$1
			;;
		-fa | --ref-fasta )	shift
			REF_FASTA=$1
			;;
		-rm | --remove ) RM=true
			;;
		-sl | --slurm )	SLURM_MODE=true
			;;
		-t | --test )	TEST_MODE=true
			;;
		-h | --help )	usage
			exit
			;;
		* )	usage
			exit 1
	esac
	shift
done

##### Control options
if [[ -z "${SIZE}" || -z "${MTIME}" || -z "${FILE_TYPE}" || -z "${DIR}" ]];then
	echo "size: ${SIZE} - mtime: ${MTIME} - file type: ${FILE_TYPE} - dir: ${DIR}"; 
	usage
	exit 1
fi
if [[ ! -d "${DIR}" ]];then
	echo "${DIR} is not a directory";
	usage
	exit 1
fi
if [ "${FILE_TYPE}" != "bam" -a "${FILE_TYPE}" != "cram" ];then
	echo "bad file type ${FILE_TYPE}, should be bam or cram";
	usage
	exit 1
fi
if [[ ! -f "${SAMTOOLS}" || ! -f "${REF_FASTA}" ]];then
	echo "samtools path ${SAMTOOLS} or ref fasta path ${REF_FASTA} is false";
	usage
	exit 1
fi

SLURM=''
if [ "${SLURM_MODE}" == true ];then
	SLURM="srun -N1 -c1 "
fi
TEST_PREFIX=''
TEST_SUFFIX=''
if [ "${TEST_MODE}" == true ];then
	TEST_PREFIX='echo "test cmd: '
	TEST_SUFFIX='"'
fi



convert () {
	FILE_LIST=$(find "${DIR}" -xdev -name "*.${FILE_TYPE}" -mtime "${MTIME}" -type f -size "${SIZE}" -exec ls "{}" \;)
	for FILE in ${FILE_LIST}; do
		OUT=$(echo "${FILE}" | sed "s/${FILE_TYPE}/${CONVERT_TYPE}/g")
		echo "Launching: ${SLURM} ${SAMTOOLS} view -T ${REF_FASTA} ${CONVERT_OPT} -o ${OUT} ${FILE}";
		${TEST_PREFIX} ${SLURM} ${SAMTOOLS} view -T ${REF_FASTA} ${CONVERT_OPT} -o ${OUT} ${FILE} ${TEST_SUFFIX}
		 if [ $? -eq 0 ];then 
			 echo "Launching: ${SLURM} ${SAMTOOLS} index ${OUT} ${OUT}${CONVERT_SUFFIX_INDEX}"
			 ${TEST_PREFIX} ${SLURM} ${SAMTOOLS} index ${OUT} ${OUT}${CONVERT_SUFFIX_INDEX} ${TEST_SUFFIX}
			 if [ $? -eq 0 -a "${RM}" == true ];then
				 RM_FILE=${FILE%.${FILE_TYPE}}
				 echo "Deleting file:  rm ${RM_FILE}${FILE_SMALL_SUFFIX}*" 
				 ${TEST_PREFIX} ${SLURM} rm ${RM_FILE}${FILE_SMALL_SUFFIX}* ${TEST_SUFFIX}
			 fi
		 fi
	 done
}

if [ "${FILE_TYPE}" == 'bam' ];then
	CONVERT_TYPE='cram'
	CONVERT_OPT='-C'
	CONVERT_SUFFIX_INDEX='.crai'
	FILE_SMALL_SUFFIX='.ba'
else
	CONVERT_TYPE='bam'
	CONVERT_OPT='-b'
	CONVERT_SUFFIX_INDEX='bai'
	FILE_SMALL_SUFFIX='cra'
fi

convert

exit

#BAM_LIST=$(find NextSeq/ -xdev -name "*.bam" -mtime +180  -type f  -size +20000000k -exec ls "{}" \;)
#for i in ${BAM_LIST}; do
#	CRAM=$(echo $i | sed "s/bam/cram/g")
#	echo "Launching: /usr/local/bin/samtools view -T /usr/local/share/refData/genome/hg19/hg19.fa -C -o ${CRAM} ${i}"
#	srun -N1 -c1 /usr/local/bin/samtools view -T /usr/local/share/refData/genome/hg19/hg19.fa -C -o ${CRAM} ${i}
#	echo "Launching: /usr/local/bin/samtools index ${CRAM} ${CRAM}.crai"
#	srun -N1 -c1 /usr/local/bin/samtools index ${CRAM} ${CRAM}.crai
#	if [ $? -eq 0 ];then
#		FILENAME=${i%.bam}
#		rm ${FILENAME}.ba*
#	fi	
#done

#BAM_LIST=$(find NextSeq/ -xdev -name "*.bam" -mtime +180  -type f  -size +20000000k -exec ls "{}" \;)
#for i in ${BAM_LIST}; do
#	rm ${i}
#	rm ${i}.bai
#done

