#!bin/bash

cd &PIPELINE_PATH
cd ..
#cd streamflow_vc/streamflow_short/short/bin/def
scp -r VariantsStream &HPC_USERNAME@&HPC_HOSTNAME:&HPC_WORK_DIRECTORY

#cd parabricks_4.0
#scp clara-parabricks4_zless.sif srawash@fe02.franklin.iit.local:/home/srawash/streamflow_vc/streamflow_short/short/bin/def/parabricks_4.0

if [ &THREE_ACTORS = true ]; then
ssh &OC_USERNAME@&OC_HOSTNAME 'if [ ! -d &OC_WORK_DIR/spring ]; then mkdir &OC_WORK_DIR/spring; fi'
scp &PATH_BIN_LOCAL/bin/def/spring/container_spring.sif &OC_USERNAME@&OC_HOSTNAME:&OC_WORK_DIR/spring
ssh &HPC_USERNAME@&HPC_HOSTNAME 'if [ ! -d &HPC_FASTQS_FOLDER ]; then mkdir &HPC_FASTQS_FOLDER; fi'
fi

