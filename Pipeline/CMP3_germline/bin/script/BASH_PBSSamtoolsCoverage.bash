#! /bin/bash

input="$1"
output_folder="$2"

########################################################
##### PARAMETERS and Function
#######################################################
CONTAINER=/hpcshare/genomics/flanduzzi/SamBcfBed_Tools.sif
BIND_INFOLDER=/home/mydata
BIND_OUTFOLDER=/home/myresult
WORK=/hpcshare/genomics/ASL_NEU/work

function extractname() {
    fullname=$1
    filename=$(basename -- "$fullname")
    name="${filename%%.*}"
    echo $name
}

########################################################
##### INPUT - Check
#######################################################
if [ ! -f ${input} ] ; then
    echo "ERROR: ${input} file not found."
    exit 1
else
        echo "INPUT BAM: $input"
fi

if [ ! -d ${output_folder} ] ; then
    echo "ERROR: ${output_folder} directory not found."
    exit 1
else
        echo "OUTPUT FOLDER: ${output_folder}"
fi

### Extract FOLDER and NAME from PATH
input_bam=$( basename -- ${input} )
input_folder=$( dirname ${input} )

########################################################
##### OUTPUT
########################################################
name=$( extractname $input_bam )
output=coverage_${name}.txt
echo "OUTPUT TXT: ${output_folder}/${output}"

cat <<EOS | qsub - 
##!/bin/bash
#PBS -N samtools_coverage
#PBS -l select=1:ncpus=1:mem=4GB
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -q cpunodes

#command line della pipeline Somatic
echo "COMMAND: singularity exec --bind ${input_folder}:${BIND_INFOLDER} --bind ${output_folder}:${BIND_OUTFOLDER} --bind ${WORK} ${CONTAINER} samtools coverage ${BIND_INFOLDER}/${input_bam} -o ${BIND_OUTFOLDER}/${output}"
StartingTime="\$( date )"

echo "Starting EXECUTION: \$( date )"
time singularity exec --bind ${input_folder}:${BIND_INFOLDER} --bind ${output_folder}:${BIND_OUTFOLDER} ${CONTAINER} samtools coverage ${BIND_INFOLDER}/${input_bam} --output ${BIND_OUTFOLDER}/${output}
#--region chrM:1-16569 
echo "Finish EXECUTION:\$( date )"
EOS
