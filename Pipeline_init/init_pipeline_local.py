import os
import sys


# Get the directory of the current script
script_directory = os.path.dirname(os.path.realpath(__file__))

# Get the parent directory
pipeline_path = os.path.dirname(script_directory)

# Executes the config file
config_file = pipeline_path + "/Pipeline_init/config_user.txt"  
with open(config_file, "r") as file:
    code = file.read()
exec(code)

# Define parabricks version and, consequently, its path
if parabricks_version == "4.0":
    parabricks_path = '/parabricks_4.0/clara-parabricks_4.sif'
    parabricks_path_def = '/parabricks_4.0/clara-parabricks4_zless.def'
elif parabricks_version == "4.1":
    parabricks_path = '/parabricks_4.1/clara-parabricks_4.sif'
    parabricks_path_def = '/parabricks_4.1/clara-parabricks4_zless.def'
else:
    print('The Parabricks version is not available, please choose between "4.0" and "4.1" (which requires gpu a100).')
    sys.exit()

# Define the streamflow script to run, depending on whether the config is 2 or 3 actors
if three_actors == "true":
    command = "streamflow run ./Pipeline/streamflow.yml"
elif three_actors == "false":
    command = "streamflow run ./Pipeline/streamflow_2act.yml"
else:
    print(f'The paramenter "three_actors" in the config file can only take value "true" or "false", you typed {three_actors}.')


path_bin_local = pipeline_path + '/Pipeline/CMP3_germline'
path_bin_HPC = hpc_work_dir + '/VariantsStream/Pipeline/CMP3_germline'

input_file_path = pipeline_path + '/Pipeline/streamflow_pre.yml'
output_file_path = pipeline_path + '/Pipeline/streamflow.yml'

input_paths = ['/Pipeline/streamflow_pre.yml',
               '/Pipeline/streamflow_2act_pre.yml',
               '/Pipeline/environments/pbs_delete_pre.pbs',
               '/Pipeline/environments/pbs_nextflow_pre.pbs',
               '/Pipeline/cwl/config_cwl_pre.yml',
               '/Pipeline/cwl/config_cwl_2act_pre.yml',
               '/Pipeline/cwl/scripts/delete_pre.sh','/Pipeline/cwl/scripts/copy_results_back_pre.sh',
               '/Pipeline/CMP3_germline/bin/hpc_set_germ_pre.config',
               '/Pipeline/CMP3_germline/bin/container_germ_pre.config',
               '/Pipeline/CMP3_germline/nextflow_pre.config',
               '/Pipeline/CMP3_germline/Makefile_pre',
               '/Pipeline/CMP3_germline/bin/def/MELT/Melt_container_pre.def',
               '/Pipeline/CMP3_germline/bin/def/Filtering/TSV_FilterAnnotation_MantaMelt_AnnotSV3.2.3_SVAfotate_pre.def',
               '/Pipeline/CMP3_germline/bin/def/sv_cons_general_utilities/sv_cons_general_utilities_pre.def',
               '/Pipeline/CMP3_germline/bin/def/VCF_AnnotationFiltering/VCF_FilterAnnotation_v9_pre.def',
               '/Pipeline_init/streamflow_pre.yml',
               '/Pipeline_init/environments/indexing_Cosmic_pre.pbs',
               '/Pipeline_init/environments/download_databases_pre.pbs',
               '/Pipeline_init/environments/add_permissions_pre.pbs',
               '/Pipeline_init/cwl/config_cwl_pre.yml',
               '/Pipeline_init/cwl/script/indexing_Cosmic_sh_pre.sh',
               '/Pipeline_init/cwl/script/copy_to_HPC_sh_pre.sh',
               '/Pipeline_init/cwl/script/add_permissions_sh_pre.sh',
               '/Pipeline_init/cwl/script/download_db_pre.pbs',
               '/Pipeline/cwl/scripts/compress_pre.sh',
               '/Pipeline/environments/pbs_compress_pre.pbs',
               '/Pipeline/cwl/scripts/move_to_HPC_pre.sh',
               '/Pipeline/cwl/scripts/move_compressed_to_hosp_pre.sh']

output_paths = ['/Pipeline/streamflow.yml',
                '/Pipeline/streamflow_2act.yml',
                '/Pipeline/environments/pbs_delete.pbs',
               '/Pipeline/environments/pbs_nextflow_nosample.pbs',
                '/Pipeline/cwl/config_cwl_no_sample.yml',
                '/Pipeline/cwl/config_cwl_2act_no_sample.yml',
               '/Pipeline/cwl/scripts/delete.sh','/Pipeline/cwl/scripts/copy_results_back_nosample.sh',
               '/Pipeline/CMP3_germline/bin/hpc_set_germ.config',
               '/Pipeline/CMP3_germline/bin/container_germ.config',
               '/Pipeline/CMP3_germline/nextflow.config',
               '/Pipeline/CMP3_germline/Makefile',
               '/Pipeline/CMP3_germline/bin/def/MELT/Melt_container.def',
               '/Pipeline/CMP3_germline/bin/def/Filtering/TSV_FilterAnnotation_MantaMelt_AnnotSV3.2.3_SVAfotate.def',
               '/Pipeline/CMP3_germline/bin/def/sv_cons_general_utilities/sv_cons_general_utilities.def',
               '/Pipeline/CMP3_germline/bin/def/VCF_AnnotationFiltering/VCF_FilterAnnotation_v9.def',
               '/Pipeline_init/streamflow.yml',
               '/Pipeline_init/environments/indexing_Cosmic.pbs',
               '/Pipeline_init/environments/download_databases.pbs',
               '/Pipeline_init/environments/add_permissions.pbs',
               '/Pipeline_init/cwl/config_cwl.yml',
               '/Pipeline_init/cwl/script/indexing_Cosmic_sh.sh',
               '/Pipeline_init/cwl/script/copy_to_HPC_sh.sh',
               '/Pipeline_init/cwl/script/add_permissions_sh.sh',
               '/Pipeline_init/cwl/script/download_db.pbs',
               '/Pipeline/cwl/scripts/compress_nosample.sh',
               '/Pipeline/environments/pbs_compress.pbs',
               '/Pipeline/cwl/scripts/move_to_HPC_nosample.sh',
               '/Pipeline/cwl/scripts/move_compressed_to_hosp_nosample.sh']




def substitution(content):
    content = content.replace('&PIPELINE_PATH', pipeline_path)
    content = content.replace('&HPC_USERNAME', hpc_username)
    content = content.replace('&HPC_HOSTNAME', hpc_hostname)
    content = content.replace('&HPC_OUTPUT_DIR', hpc_output_dir)
    content = content.replace('&PROJECT_NAME', project_name)
    content = content.replace('&EMAIL', email)
    content = content.replace('&CONDA_MODULE', conda_module)
    content = content.replace('&SINGULARITY_MODULE', singularity_module)
    content = content.replace('&HPC_QUEUE_PARABRICKS', hpc_queue_parabricks)
    content = content.replace('&PARABRICKS_PATH_DEF', parabricks_path_def)
    content = content.replace('&PARABRICKS_PATH', parabricks_path)
    content = content.replace('&HPC_QUEUE', HPC_queue)
    content = content.replace('&HPC_WORK_DIRECTORY', hpc_work_dir)
    content = content.replace('&HPC_FASTQS_FOLDER', hpc_fastqs_folder)
    content = content.replace('&PATH_BIN_LOCAL', path_bin_local)
    content = content.replace('&PATH_BIN_HPC', path_bin_HPC)
    content = content.replace('&OC_USERNAME', oc_username)
    content = content.replace('&OC_HOSTNAME', oc_hostname)
    content = content.replace('&OC_WORK_DIR', oc_work_dir)
    content = content.replace('&OC_FASTQ_DIR', oc_fastq_dir)
    content = content.replace('&THREE_ACTORS', three_actors)
    content = content.replace('&MOVE_FASTQS', move_fastqs)
    content = content.replace('&SSH_OC_HPC', ssh_oc_hpc)
    content = content.replace('&COMMAND', command)

    return content


for i in range(len(input_paths)):
    input_file_path =  pipeline_path + input_paths[i]
    output_file_path =  pipeline_path + output_paths[i]
    with open(input_file_path, 'r') as input_file:
        content = input_file.read()

    content = substitution(content)

    # Write the modified content back to the file
    with open(output_file_path, 'w') as output_file:
        output_file.write(content)

    print(f"Substitution complete. Updated file written to {output_file_path}")




