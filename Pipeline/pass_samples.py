import os
import argparse

def substitution(content, sample_list):
    content = content.replace('&SAMPLE_LIST', sample_list)
    return content

def main(args):
    # Get the directory of the current script
    script_directory = os.path.dirname(os.path.realpath(__file__))

    # Get the parent directory
    pipeline_path = os.path.dirname(script_directory)

    input_file_path = [pipeline_path + '/Pipeline/cwl/config_cwl_no_sample.yml',
                       pipeline_path + '/Pipeline/cwl/config_cwl_2act_no_sample.yml',
                       pipeline_path + '/Pipeline/cwl/scripts/compress_nosample.sh',
                       pipeline_path + '/Pipeline/cwl/scripts/copy_results_back_nosample.sh',
                       pipeline_path + '/Pipeline/environments/pbs_nextflow_nosample.pbs',
                       pipeline_path + '/Pipeline/cwl/scripts/move_to_HPC_nosample.sh',
                       pipeline_path + '/Pipeline/cwl/scripts/move_compressed_to_hosp_nosample.sh']
    output_file_path = [pipeline_path + '/Pipeline/cwl/config_cwl.yml',
                        pipeline_path + '/Pipeline/cwl/config_cwl_2act.yml',
                        pipeline_path + '/Pipeline/cwl/scripts/compress.sh',
                        pipeline_path + '/Pipeline/cwl/scripts/copy_results_back.sh',
                        pipeline_path + '/Pipeline/environments/pbs_nextflow.pbs',
                       pipeline_path + '/Pipeline/cwl/scripts/move_to_HPC.sh',
                       pipeline_path + '/Pipeline/cwl/scripts/move_compressed_to_hosp.sh']
    
    for i in range(len(input_file_path)):

        with open(input_file_path[i], 'r') as input_file:
            content = input_file.read()

        content = substitution(content, args.sample_list)

        # Write the modified content back to the file
        with open(output_file_path[i], 'w') as output_file:
            output_file.write(content)

        print(f"Substitution complete. Updated file written to {output_file_path[i]}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Substitute sample_list in a YAML file.")
    parser.add_argument("sample_list", type=str, help="Comma-separated list of sample IDs")

    args = parser.parse_args()
    main(args)

