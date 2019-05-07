#!/bin/bash

###########################################################################
#########							###########
######### bamCramConverter: bcc						###########
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
Program: bamCramConverter: bcc
Version: ${VERSION}
Contact: Baux David <david.baux@inserm.fr>
Usage: bash bcc.sh -d <path to search dir> -s <File size to search (man find -size)> -mt <File last modif to search (man find -mtime)> -f <bam|cram> -st <path to samtools> -fa <path to ref genome .fa>
	"

if [ $# -eq 0 ]; then
	echo "${USAGE}"
	echo "Error Message : Arguments not provided"
	echo ""
	exit 1
fi

usage ()
{
	echo 'This script convert bam files to cram and reverse.'
	echo 'Usage : bash bcc.sh'
	echo '	Mandatory arguments :'
	echo '		* -d|--directory	<path to search dir>: root dir for find command'
	echo '		* -s|--size		<File size to search (man find -size)>: ex: +200000000k will search for files greater than 20Go; see man find -size argument'
	echo '		* -mt|--modif-time	<File last modif to search (man find -mtime): ex: +180 will search for files older than 6 months; see man find -mtime argument'
	echo '		* -f|--file-type	<bam|cram>: file type to find and convert from (bam will search for bam files and convert to cram'
	echo '	Optional arguments :'
	echo '		* -rm|--remove		: removes original file and index (in case of full conversion successi implies bam2cram-check) - default: false'
	echo '		* -st|--samtools	<path to samtools> - default: try to locate in PATH'
	echo '		* -fa|--ref-fasta	<path to ref genome .fa>: path to a fasta file reference genome (the directory containing the fasta file must also contain samtools index) - default: /usr/local/share/refData/genome/hg19/hg19.fa'
	echo '		* -c|--check		: uses bam2cram-check (slightly modified) to check the conversion - implicitely included with -rm - if fails and -rm: rm canceled) - requires python >3.5 and samtools > 1.3'
	echo '		* -p|--python3		<path to python3> - used in combination with -c or -rm: needed to run submodule bam2cram-check - default: /usr/bin/python3 - python version must be > 3.5'
	echo '		* -uc|--use-crumble	: uses crumble to compress the converted BAM/CRAM file - Note: a file that already contains "crumble" in its name will not be converted again'
	echo '		* -cp|--crumble-path	<path to crumble> - used in combination with -uc: needed to run crumble - default: try to locate in PATH'
	echo "        	* -v | --verbosity 	<integer> : decrease or increase verbosity level (ERROR : 1 | WARNING : 2 | INFO : 3 | COMMAND [default] : 4 | DEBUG : 5)"
	echo '	General arguments :'
	echo '		* -sl|--slurm   	: when running in SLURM environnment, generates srun commands - default: false'
	echo '		* -th|--threads 	: number of threads to be used for samtools -@ option (0 => 1 total thread, 1 => 2 total threads...)'
	echo '		* -h			: show this help message and exit'
	echo '		* -t			: test mode (dont execute command just print them) - default: false'
	echo ''
	exit
}

# --Option processing
SLURM_MODE=false
THREADS=0
TEST_MODE=false
RM=false
CHECK=false
SAMTOOLS=$(which samtools)
REF_FASTA='/usr/local/share/refData/genome/hg19/hg19.fa'
GENOME_VERSION='hg19'
CHR1_LEN=249250621
PYTHON3='/usr/bin/python3'
CRUMBLE=$(which crumble)
CRUMBLE_MOTIF='crumble'
#CRUMBLE=/usr/local/bin/crumble
# -- Script log 

VERBOSITY=4


# -- Log variables 

ERROR=1
WARNING=2
INFO=3
CMD=4
DEBUG=5

# -- Log functions 

error() { log ${ERROR} "[error] : $1" ; }
warning() { log ${WARNING} "[warning] : $1" ;}
cmd() { log ${CMD} "[command] :\n$1" ; }
info() { log ${INFO} "[info] : $1" ; }
debug() { log ${DEBUG} "[debug] : $1" ; }

# -- Print log 

echoerr() { echo -e "$@" 1>&2 ; }

log() {

	if [ ${VERBOSITY} -ge $1 ]
		then
			echoerr "[`date +'%Y-%m-%d %H:%M:%S'`] - $2"
	fi
}





# -- Parse command line
while [ "$1" != "" ];do
	case $1 in
		-c | --check)	CHECK=true			
			;;
		-d | --directory)	shift
			DIR=$1
			;;
		-s | --size)	shift
			SIZE=$1
			;;
		-mt | --modif-time)	shift
			MTIME=$1
			;;
		-f | --file-type)	shift
			FILE_TYPE=$1
			;;
		-st | --samtools)	shift
			SAMTOOLS=$1
			;;
		-fa | --ref-fasta)	shift
			REF_FASTA=$1
			;;
		-p | --python3)	shift
			PYTHON3=$1
			;;
		-uc | --use-crumble)	USE_CRUMBLE=true
			;;
		-cp | --crumble-path)	shift
			CRUMBLE=$1
			;;
		-rm | --remove) RM=true
			;;
		-sl | --slurm)	SLURM_MODE=true
			;;
		-th | --threads)	shift
			THREADS=$1
			;;
		-t | --test)	TEST_MODE=true
			;;
		-v | --verbosity) shift 
			# Check if verbosity level argument is an integer before assignment 
			if ! [[ "$1" =~ ^[0-9]+$ ]]
			then 
				error "\"$1\" must be an integer !"
				echo " "
				help 
			else 
				VERBOSITY=$1
				((VERBOSITYCOUNTER++))
			fi 
			;;
		-h | --help)	usage
			exit
			;;
		* )	usage
			exit 1
	esac
	shift
done


info "BamCramConverter version ${VERSION}:"

##### Control options
if [[ -z "${SIZE}" || -z "${MTIME}" || -z "${FILE_TYPE}" || -z "${DIR}" ]];then
	#echo "ERROR - [`date +'%Y-%m-%d %H:%M:%S'`] - size: ${SIZE} - mtime: ${MTIME} - file type: ${FILE_TYPE} - dir: ${DIR}"
	error " size: ${SIZE} - mtime: ${MTIME} - file type: ${FILE_TYPE} - dir: ${DIR}"
	usage
	exit 1
fi
if [[ ! -d "${DIR}" ]];then
	error " ${DIR} is not a directory"
	usage
	exit 1
fi
if [ "${FILE_TYPE}" != "bam" -a "${FILE_TYPE}" != "cram" ];then
	error " bad file type ${FILE_TYPE}, should be bam or cram"
	usage
	exit 1
fi
if [[ ! -x "${SAMTOOLS}" || ! -r "${REF_FASTA}" ]];then
	error " samtools path ${SAMTOOLS} or ref fasta path ${REF_FASTA} is false"
	usage
	exit 1
fi
if [[ ! "${THREADS}" =~ ^[0-9]+$ ]];then
	error " -th option must be an integer"
	usage
	exit 1
fi
if [[ ! -x "${PYTHON3}" ]] && [[ "$RM{}" == 'true' || "${CHECK}" == 'true' ]];then
	error " Could not find Python 3 executable ${PYTHON3}"
	usage
	exit 1
fi
if [[ ! -x "${CRUMBLE}" && "${USE_CRUMBLE}" == 'true' ]];then
	error " Could not find Crumble executable ${CRUMBLE}"
	usage
	exit 1
fi


##### Modif variables
### RM modifies CHECK
if [[ "${RM}" == true ]];then
	CHECK=true
fi
### samtools threads and slurm (samtools threads are  0-based
SAMTOOLS_MULTI=''
if [ "${THREADS}" -gt 0 ];then
	SAMTOOLS_MULTI="-@ ${THREADS}"
fi
SLURM_THREADS=$((THREADS+1))
SLURM=''
SLURM_MULTI=''
if [ "${SLURM_MODE}" == true ];then
	SLURM="srun -N1 -c1 "
	SLURM_MULTI="${SLURM}"
	if [ "${THREADS}" -gt 0 ];then
		SLURM_MULTI="srun -N1 -c${SLURM_THREADS} "
	fi
fi


##### Debug
debug "DIR:${DIR}"
debug "FILE TYPE TO BE CONVERTED:${FILE_TYPE}"
debug "MIN FILE SIZE:${SIZE}"
debug "LAST MODIFICATION TIME:${MTIME} DAYS"
debug "SAMTOOLS:${SAMTOOLS}"
debug "THREADS:${THREADS}"
debug "REMOVE ORIGINAL FILE:${RM} - if true, automatically set next variable to true"
debug "CHECK CONVERSION:${CHECK}"
debug "PYTHON 3:${PYTHON3}"
debug "USE CRUMBLE:${USE_CRUMBLE}"
debug "CRUMBLE:${CRUMBLE}"
debug "FASTA REF:${REF_FASTA}"
debug "USE SLURM:${SLURM_MODE}"
debug "TEST MODE:${TEST_MODE}"


check () {
	OUT_DIR=$(dirname "${1}")
	if [ ${1##*\.} == 'bam' ];then
		BAM="${1}"
		CRAM="${2}"
	else
		CRAM="${1}"
		BAM="${2}"
	fi
	cmd "\"${PYTHON3}\" bam2cram-check/main.py -b \"${BAM}\" -c \"${CRAM}\" -r \"${3}\" -s -e \"${OUT_DIR}/bam2cramcheck_error.txt\" --log \"${OUT_DIR}/bam2cramcheck_log.txt\""
	if [ "${TEST_MODE}" == false ];then
		"${PYTHON3}" bam2cram-check/main.py -b "${BAM}" -c "${CRAM}" -r "${3}" -s -e "${OUT_DIR}/bam2cramcheck_error.txt" --log "${OUT_DIR}/bam2cramcheck_log.txt"
		if [ $? -eq 0 ];then
			if grep -q 'There were no errors and no differences between the stats for the 2 files' "${OUT_DIR}/bam2cramcheck_log.txt";then
				echo "bam2cram-check successfull"
			else
				echo "bam2cram-check failed: check ${OUT_DIR}/bam2cramcheck_error.txt"
			fi
		else
			echo "bam2cram-check failed: check ${OUT_DIR}/bam2cramcheck_error.txt"
		fi
	else
		echo "bam2cram-check successfull"
	fi
}

convert () {
	FILE_LIST=$(find "${DIR}" -xdev -name "*.${FILE_TYPE}" -mtime "${MTIME}" -type f -size "${SIZE}" -exec ls "{}" \;)
	for FILE in ${FILE_LIST}; do
		#checking genome version - as genome version is not always explicitely cited in BAM headers (see LRM outputs), we wall check crh1 length
		#hg19 chr1 : 249250621
		#hg38 chr1 : 248956422
		#then bcc must be aware of the REF_FASTA version : will look for 19 or 38 in the name
		#must be updated with new assemblies
		TEST_REF=$(echo "${REF_FASTA}" | grep hg38)
		debug "BEFORE ASSIGNMENT REF:${REF_FASTA}-TESTREF:${TEST_REF}-GENOME:${GENOME_VERSION}-CHR1:${CHR1_LEN}"
		if [ "${TEST_REF}" != '' ];then
			GENOME_VERSION=hg38
			CHR1_LEN=248956422
		fi
		TEST_REF=$("${SAMTOOLS}" view -H "${FILE}" | grep "LN:${CHR1_LEN}")
		debug "AFTER ASSIGNMENT REF:${REF_FASTA}-TESTREF:${TEST_REF}-GENOME:${GENOME_VERSION}-CHR1:${CHR1_LEN}"
		if [ "${TEST_REF}" == '' ];then
			warning " ${FILE} genome version does not match ${GENOME_VERSION} based on chr1 length" # -${TEST_REF}-${GENOME_VERSION}-${CHR1_LEN}"
			warning " ${FILE} won't be converted"
			continue
		fi
		info "Starting conversion of \"${FILE}\" into ${CONVERT_TYPE}"
		BASE_DIR=$(dirname "$FILE")
		FILE_NAME=$(basename "$FILE")
		OUT_FILE=$(echo "${FILE_NAME}" | sed "s/\.${FILE_TYPE}/\.${CONVERT_TYPE}/")
		OUT="${BASE_DIR}/${OUT_FILE}"
		PER_FILE_USE_CRUMBLE="${USE_CRUMBLE}"
		#look for 'crumble' motif in file name not to recrumble a file
		if [[ "${USE_CRUMBLE}" == 'true' ]] && [[ "${FILE_NAME}" =~ ${CRUMBLE_MOTIF} ]];then
			PER_FILE_USE_CRUMBLE=false
			warning "Crumble is activated but file name ${FILE_NAME} already contains 'crumble' motif therefore crumble will be deactivated for this file"
		fi
		if [ -f "${OUT}" ]; then
			#file already exists
			warning "${OUT} already exists - conversion aborted"
			continue
		fi
		info "Launching samtools view"
		cmd "${SLURM_MULTI} \"${SAMTOOLS}\" view -T \"${REF_FASTA}\" ${CONVERT_OPT} ${SAMTOOLS_MULTI} -o \"${OUT}\" \"${FILE}\""
		if [ "${TEST_MODE}" == false ];then
			${SLURM_MULTI} "${SAMTOOLS}" view -T "${REF_FASTA}" ${CONVERT_OPT} ${SAMTOOLS_MULTI} -o "${OUT}" "${FILE}"
		fi
		if [ $? -eq 0 ];then
			info "Conversion successfull - Launching samtools index"
			cmd "${SLURM} \"${SAMTOOLS}\" index \"${OUT}\" \"${OUT}${CONVERT_SUFFIX_INDEX}\""
			if [ "${TEST_MODE}" == false ];then
				${SLURM} "${SAMTOOLS}" index "${OUT}" "${OUT}${CONVERT_SUFFIX_INDEX}"
			fi
			if [ $? -eq 0 ];then
				if [ "${CHECK}" == true ];then
					info "Indexing sucessfull - Checking files"
					CHECK_VALUE=$(check "${FILE}" "${OUT}" "${REF_FASTA}")
					info "${CHECK_VALUE}"
				fi
				if [ "${RM}" == true ];then
					#info "Indexing sucessfull - Checking files"
					#CHECK_VALUE=$(check "${FILE}" "${OUT}" "${REF_FASTA}")
					if [[ "${CHECK_VALUE}" =~ "bam2cram-check successfull" ]];then
						#info "${CHECK_VALUE}"
						RM_FILE=${FILE%.${FILE_TYPE}}
						info "Removing \"${RM_FILE}${FILE_SMALL_SUFFIX}\"*" 
						if [ "${TEST_MODE}" == false ];then
							${SLURM} rm "${RM_FILE}${FILE_SMALL_SUFFIX}"*
						fi
					else
						#echo "${CHECK_VALUE}"
						warning "bam2cram-check failed: rm canceled"
					fi
				elif [ "${CHECK}" != true ];then
					info "Indexing sucessfull"
					#no check so we pass CHECK_VALUE to success
					CHECK_VALUE='bam2cram-check successfull'
				fi
				#elif [ "${CHECK}" == true ];then
				#	info "Indexing sucessfull - checking files:"
				#	CHECK_VALUE=$(check "${FILE}" "${OUT}" "${REF_FASTA}")
				#	info "${CHECK_VALUE}"
				#crumble
				if [[ "${PER_FILE_USE_CRUMBLE}" == true && "${CHECK_VALUE}" =~ "bam2cram-check successfull" ]];then
					CRUMBLE_OUT=$(echo "${OUT}" | sed "s/\.${CONVERT_TYPE}$/\.${CRUMBLE_MOTIF}\.${CONVERT_TYPE}/")
					info "Launching crumble"
					cmd "${SLURM_MULTI} \"${CRUMBLE}\" -O \"${CONVERT_TYPE},nthreads=${SLURM_THREADS}\" \"${OUT}\" \"${CRUMBLE_OUT}\""
					if [ "${TEST_MODE}" == false ];then
						#${TEST_PREFIX} ${SLURM_MULTI} ${CRUMBLE} -O ${CONVERT_TYPE},nthreads=${SLURM_THREADS} ${OUT} ${CRUMBLE_OUT} ${TEST_SUFFIX}
						${SLURM_MULTI} "${CRUMBLE}" -O "${CONVERT_TYPE},nthreads=${SLURM_THREADS}" "${OUT}" "${CRUMBLE_OUT}"
					fi
					if [ $? -eq 0 ];then
						info "Crumble compression successfull - Launching samtools index"
						cmd "${SLURM} \"${SAMTOOLS}\" index \"${CRUMBLE_OUT}\" \"${CRUMBLE_OUT}${CONVERT_SUFFIX_INDEX}\""
						if [ "${TEST_MODE}" == false ];then
							#${TEST_PREFIX} ${SLURM} ${SAMTOOLS} index ${CRUMBLE_OUT} ${CRUMBLE_OUT}${CONVERT_SUFFIX_INDEX} ${TEST_SUFFIX}
							${SLURM} "${SAMTOOLS}" index "${CRUMBLE_OUT}" "${CRUMBLE_OUT}${CONVERT_SUFFIX_INDEX}"
						fi
						if [ $? -eq 0 ];then
							info "Indexing sucessfull - removing \"${OUT}\""
							if [ "${TEST_MODE}" == false ];then
								#${TEST_PREFIX} ${SLURM} rm ${OUT} ${OUT}${CONVERT_SUFFIX_INDEX} ${TEST_SUFFIX}
								${SLURM} rm "${OUT}" "${OUT}${CONVERT_SUFFIX_INDEX}"
							fi
						else
							warning "Error while indexing ${CRUMBLE_OUT} - ${OUT} won't be erased"
						fi
					else
						warning "Error while compressing ${OUT} - ${OUT} won't be erased"
					fi
				fi	
			else
				warning "Error while indexing ${OUT}"
			fi
		else
			warning "Error while converting ${FILE} into ${CONVERT_TYPE}"
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
	CONVERT_SUFFIX_INDEX='.bai'
	FILE_SMALL_SUFFIX='.cra'
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

