# #!/bin/bash

 rm -r &HPC_OUTPUT_DIR/streamflow_res
 rm &HPC_OUTPUT_DIR/streamflow_CMP3.* 
 rm &HPC_OUTPUT_DIR/nf_del.*

 if [ &THREE_ACTORS = true ]; then
    rm -r &HPC_FASTQS_FOLDER/*
 fi
