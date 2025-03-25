
#!/bin/bash


# Define the comma-separated string of samples
#samples="test1,test2,test3"
samples="&SAMPLE_LIST"
# Split the string into an array using ',' as the delimiter
IFS=',' read -ra sample_array <<< "$samples"


cd &PIPELINE_PATH/Pipeline


if [ &MOVE_FASTQS = true ]; then
    if [ ! -d &PIPELINE_PATH/Pipeline/compressed_&SAMPLE_LIST ]; then mkdir compressed_&SAMPLE_LIST; fi
    cd compressed_&SAMPLE_LIST
    if [ &THREE_ACTORS = true ]; then
        for test in "${sample_array[@]}"; do
            scp -r &OC_USERNAME@&OC_HOSTNAME:&OC_WORK_DIR/${test}_spring .
            sleep 2
            ssh &HPC_USERNAME@&HPC_HOSTNAME "rm &OC_WORK_DIR/${test}_spring"
            sleep 5
        done
    else
        for test in "${sample_array[@]}"; do
            scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_FASTQS_FOLDER/${test}_spring .
            sleep 2
            ssh &HPC_USERNAME@&HPC_HOSTNAME "rm &HPC_FASTQS_FOLDER/${test}_spring"
            sleep 5
        done
    fi
else
    # No operations are performed when MOVE_FASTQ is not true
    :
fi






