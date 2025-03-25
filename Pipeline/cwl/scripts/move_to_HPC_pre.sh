
#!/bin/bash



#samples="test1,test2,test3"
samples="&SAMPLE_LIST"
# Split the string into an array using ',' as the delimiter
IFS=',' read -ra sample_array <<< "$samples"



if [ &THREE_ACTORS = true ]; then
    ssh &HPC_USERNAME@&HPC_HOSTNAME "if [ ! -d &HPC_FASTQS_FOLDER ]; then mkdir -p &HPC_FASTQS_FOLDER ; fi "
    if [ &SSH_OC_HPC = true ]; then
        # If there is an SSH key and SSH_OC_HPC is true
        for test in "${sample_array[@]}"; do
            ssh &OC_USERNAME@&OC_HOSTNAME "scp &OC_FASTQ_DIR/*${test}_R1_001.fastq.gz &OC_FASTQ_DIR/*${test}_R2_001.fastq.gz &HPC_USERNAME@&HPC_HOSTNAME:&HPC_FASTQS_FOLDER"
        done
    else
        for test in "${sample_array[@]}"; do
            scp &OC_USERNAME@&OC_HOSTNAME:&OC_FASTQ_DIR/*${test}_R1_001.fastq.gz &PIPELINE_PATH
            scp &OC_USERNAME@&OC_HOSTNAME:&OC_FASTQ_DIR/*${test}_R2_001.fastq.gz &PIPELINE_PATH
        done
        sleep 2
        for test in "${sample_array[@]}"; do
            scp &PIPELINE_PATH/*${test}_R1_001.fastq.gz &HPC_USERNAME@&HPC_HOSTNAME:&HPC_FASTQS_FOLDER 
            scp &PIPELINE_PATH/*${test}_R2_001.fastq.gz &HPC_USERNAME@&HPC_HOSTNAME:&HPC_FASTQS_FOLDER 
        done
        sleep 2
        for test in "${sample_array[@]}"; do
            rm &PIPELINE_PATH/*${test}_R1_001.fastq.gz 
            rm &PIPELINE_PATH/*${test}_R2_001.fastq.gz
        done
    fi
else
    #No operation to be performed
    :
fi


