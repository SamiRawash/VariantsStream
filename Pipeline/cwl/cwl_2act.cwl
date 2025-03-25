cwlVersion: v1.2
class: Workflow

inputs:
  nf_script: File
  nf_config: File
  command_C: string
  command_run: string
  del_file: File
  command1: string
  command2: string
  command3: string
  command4: string
  command5: string
  command6: string
  copy_results_back: File

outputs:
  end:
    type: string
    outputSource: delete/end  
  
  



steps:
#-------------------------------------#


# This step runs the variant calling pipeline with Nextflow

  pbs_submission:
    run:
      class: CommandLineTool
      baseCommand: nextflow
      inputs:
        command_C:
          type: string
          inputBinding:
            position: 1
        nf_config:
          type: File
          inputBinding:
            position: 2
        command_run:
          type: string
          inputBinding:
            position: 3
        nf_script:
          type: File
          inputBinding:
            position: 4
        command1:
          type: string
          inputBinding:
            position: 5
        command2:
          type: string
          inputBinding:
            position: 6
        command3:
          type: string
          inputBinding:
            position: 7
        command4:
          type: string
          inputBinding:
            position: 8
        command5:
          type: string
          inputBinding:
            position: 9
        command6:
          type: string
          inputBinding:
            position: 10
      outputs:
        hold_pbs:
          type: string

    in:
      command_C: command_C
      nf_config: nf_config
      nf_script: nf_script
      command1: command1
      command2: command2
      command3: command3
      command4: command4
      command5: command5
      command6: command6
      command_run: command_run

    out: [hold_pbs]

#-----------------------------------#=


# This is a waiting step, to undeploy the first environment and returning the output files
# to the client before cleaning the cluster

  sleep:
    run:
      class: CommandLineTool
      baseCommand: [bash]
      inputs:
        hold_pbs:
          type: string
        copy_results_back:
          type: File
          inputBinding:
            position: 1
      outputs:
        hold_sleep:
          type: string

    in:
      hold_pbs: pbs_submission/hold_pbs
      copy_results_back: copy_results_back
    out: [hold_sleep]

#-----------------------------------#=


# Deletes the files from the HPC

  delete:
    run:
      class: CommandLineTool
      baseCommand: [bash]
      inputs:
        hold_sleep:
          type: string
        del_file:
          type: File
          inputBinding:
            position: 1
      outputs:
        end:
          type: string

    in:
      hold_sleep: sleep/hold_sleep
      del_file: del_file
    out: [end]
