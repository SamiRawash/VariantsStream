cd &PIPELINE_PATH/Pipeline
if [ ! -d &PIPELINE_PATH/Pipeline/results_&SAMPLE_LIST ]; then mkdir results_&SAMPLE_LIST; fi
cd results_&SAMPLE_LIST
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/00_Coverage .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/01_SnpEff .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/02_cosmic_annotation .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/04_filter_impact .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/05_filter_pass .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/06_cancer_gene_census .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/07_InterVar .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/08_Tables_Germline .
scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/09_IncidentalFindings_Tables .

samples=&SAMPLE_LIST
# Split the string into an array using ',' as the delimiter
IFS=',' read -ra sample_array <<< "$samples"

#if [ &MOVE_FASTQS = true ]; then
for test in "${sample_array[@]}"; do
    scp -r &HPC_USERNAME@&HPC_HOSTNAME:&HPC_OUTPUT_DIR/streamflow_res/${test} .
done
#fi

