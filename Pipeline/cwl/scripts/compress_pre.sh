
#!/bin/bash

#DONE

# Define the comma-separated string of samples
#samples="test1,test2,test3"

samples="&SAMPLE_LIST"
# Split the string into an array using ',' as the delimiter
IFS=',' read -ra sample_array <<< "$samples"

if [ &MOVE_FASTQS = true ]; then
    if [ &THREE_ACTORS = true ]; then
        # Loop over the samples
        for test in "${sample_array[@]}"; do
            # Construct the file names and execute the command
            singularity exec &OC_WORK_DIR/spring/container_spring.sif spring -c -i &OC_FASTQ_DIR/*${test}_R1_001.fastq.gz &OC_FASTQ_DIR/*${test}_R2_001.fastq.gz -o &OC_WORK_DIR/${test}_spring -g
            sleep 2
        done

        rm &OC_WORK_DIR/streamflow_CMP3.*
        rm -r &OC_WORK_DIR/streamflow_res
        rm &OC_WORK_DIR/*_spring
    else
        # Loop over the samples
        for test in "${sample_array[@]}"; do
            # Construct the file names and execute the command
            singularity exec &HPC_WORK_DIRECTORY/VariantsStream/Pipeline/CMP3_germline/bin/def/spring/container_spring.sif spring -c -i &HPC_FASTQS_FOLDER/*${test}_R1_001.fastq.gz &HPC_FASTQS_FOLDER/*${test}_R2_001.fastq.gz -o &HPC_FASTQS_FOLDER/${test}_spring -g
            sleep 2
        done
    fi
fi
