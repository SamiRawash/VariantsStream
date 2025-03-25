cwlVersion: v1.2
class: Workflow

inputs:
  copy_to_HPC_sh: File
  command_f: string
  makefile_singularity: File
  download_db: File
  command_make_container: string
  indexing_Cosmic_sh: File
  add_permissions_sh : File


outputs:
  add_permissions_done:
    type: string
    outputSource: add_permissions/add_permissions_done

steps:
#-------------------------------------#

#This step builds the singularities
  build_singularity:
    run:
      class: CommandLineTool
      baseCommand: make
      inputs:
        command_f:
          type: string
          inputBinding:
            position: 1
        makefile_singularity:
          type: File
          inputBinding:
            position: 2
        command_make_container:
          type: string
          inputBinding:
            position: 3
      outputs:
        build_singularity_done:
          type: string


    in:
      command_f: command_f
      makefile_singularity: makefile_singularity
      command_make_container: command_make_container
    out: [build_singularity_done]



#-----------------------------------#


# This is a waiting step, before copying project to hpc

  sleep_before_copying:
    run:
      class: CommandLineTool
      baseCommand: [sleep, '1']
      inputs:
        build_singularity_done:
          type: string
      outputs:
        sleep_before_copying_done:
          type: string
#      stdout: hold_HPC2

    in:
      build_singularity_done: build_singularity/build_singularity_done
    out: [sleep_before_copying_done]



#-----------------------------------#

# This step builds the sif files and sends them to the HPC

  copy_to_HPC:
    run:
      class: CommandLineTool
      baseCommand: bash
      inputs:
        copy_to_HPC_sh:
          type: File
          inputBinding:
            position: 1
        sleep_before_copying_done:
          type: string
      outputs:
       copy_to_HPC_done:
         type: string

    in:
      copy_to_HPC_sh: copy_to_HPC_sh
      sleep_before_copying_done: sleep_before_copying/sleep_before_copying_done
    out: [copy_to_HPC_done]

#-----------------------------------#=

# This is a waiting step, after copying project HPC

  sleep_after_copying:
    run:
      class: CommandLineTool
      baseCommand: [sleep, '1']
      inputs:
        copy_to_HPC_done:
          type: string
      outputs:
        sleep_after_copying_done:
          type: string
#      stdout: hold_HPC2

    in:
      copy_to_HPC_done: copy_to_HPC/copy_to_HPC_done
    out: [sleep_after_copying_done]



#-----------------------------------#

#This step downloads the databases
  download_databases:
    run:
      class: CommandLineTool
      baseCommand: qsub
      inputs:
        download_db:
          type: File
          inputBinding:
            position: 1
        sleep_after_copying_done:
          type: string
      outputs:
        download_databases_done:
          type: string


    in:
      download_db: download_db
      sleep_after_copying_done: sleep_after_copying/sleep_after_copying_done
    out: [download_databases_done]


#-----------------------------------#

#This step indexes the Cosmic DB
  indexing_Cosmic:
    run:
      class: CommandLineTool
      baseCommand: bash
      inputs:
        indexing_Cosmic_sh:
          type: File
          inputBinding:
            position: 1
        download_databases_done:
          type: string
      outputs:
        indexing_Cosmic_done:
          type: string


    in:
      indexing_Cosmic_sh: indexing_Cosmic_sh
      download_databases_done: download_databases/download_databases_done
    out: [indexing_Cosmic_done]

#-----------------------------------#

#This step adds permissions in get_pass_variants
  add_permissions:
    run:
      class: CommandLineTool
      baseCommand: bash
      inputs:
        add_permissions_sh:
          type: File
          inputBinding:
            position: 1
        indexing_Cosmic_done:
          type: string
      outputs:
        add_permissions_done:
          type: string


    in:
      add_permissions_sh: add_permissions_sh
      indexing_Cosmic_done: indexing_Cosmic/indexing_Cosmic_done
    out: [add_permissions_done]















# cwlVersion: v1.2
# class: Workflow

# inputs:
#   #singularities_sh: File
#   #HPC_download_py: File
#   #HPC_download_py2: File
#   copy_to_HPC_sh: File
#   command_f: string
#   makefile_singularity: File
#   command_make_container: string
#   command_make_database: string
#   indexing_Cosmic_sh: File
#   add_permissions_sh : File

# #outputs: [download_databases_done, indexing_Cosmic_done, add_permissions_done]

# outputs:
#   indexing_Cosmic_done:
#     type: string
#     outputSource: indexing_Cosmic/indexing_Cosmic_done
#   add_permissions_done:
#     type: string
#     outputSource: add_permissions/add_permissions_done

# steps:
# #-------------------------------------#

# #This step builds the singularities
#   build_singularity:
#     run:
#       class: CommandLineTool
#       baseCommand: make
#       inputs:
#         command_f:
#           type: string
#           inputBinding:
#             position: 1
#         makefile_singularity:
#           type: File
#           inputBinding:
#             position: 2
#         command_make_container:
#           type: string
#           inputBinding:
#             position: 3
#       outputs:
#         build_singularity_done:
#           type: string


#     in:
#       command_f: command_f
#       makefile_singularity: makefile_singularity
#       command_make_container: command_make_container
#     out: [build_singularity_done]



# #-----------------------------------#


# # This is a waiting step, before copying project to hpc

#   sleep_before_copying:
#     run:
#       class: CommandLineTool
#       baseCommand: [sleep, '1']
#       inputs:
#         build_singularity_done:
#           type: string
#       outputs:
#         sleep_before_copying_done:
#           type: string
# #      stdout: hold_HPC2

#     in:
#       build_singularity_done: build_singularity/build_singularity_done
#     out: [sleep_before_copying_done]



# #-----------------------------------#

# # This step builds the sif files and sends them to the HPC

#   copy_to_HPC:
#     run:
#       class: CommandLineTool
#       baseCommand: bash
#       inputs:
#         copy_to_HPC_sh:
#           type: File
#           inputBinding:
#             position: 1
#         sleep_before_copying_done:
#           type: string
#       outputs:
#        copy_to_HPC_done:
#          type: string

#     in:
#       copy_to_HPC_sh: copy_to_HPC_sh
#       sleep_before_copying_done: sleep_before_copying/sleep_before_copying_done
#     out: [copy_to_HPC_done]

# #-----------------------------------#=

# # This is a waiting step, after copying project HPC

#   sleep_after_copying:
#     run:
#       class: CommandLineTool
#       baseCommand: [sleep, '1']
#       inputs:
#         copy_to_HPC_done:
#           type: string
#       outputs:
#         sleep_after_copying_done:
#           type: string
# #      stdout: hold_HPC2

#     in:
#       copy_to_HPC_done: copy_to_HPC/copy_to_HPC_done
#     out: [sleep_after_copying_done]



# #-----------------------------------#

# #This step downloads the databases
#   download_databases:
#     run:
#       class: CommandLineTool
#       baseCommand: make
#       inputs:
#         command_f:
#           type: string
#           inputBinding:
#             position: 1
#         makefile_singularity:
#           type: File
#           inputBinding:
#             position: 2
#         command_make_database:
#           type: string
#           inputBinding:
#             position: 3
#         sleep_after_copying_done:
#           type: string
#       outputs:
#         download_databases_done:
#           type: string


#     in:
#       command_f: command_f
#       makefile_singularity: makefile_singularity
#       command_make_database: command_make_database
#       sleep_after_copying_done: sleep_after_copying/sleep_after_copying_done
#     out: [download_databases_done]


# #-----------------------------------#

# #This step indexes the Cosmic DB
#   indexing_Cosmic:
#     run:
#       class: CommandLineTool
#       baseCommand: bash
#       inputs:
#         indexing_Cosmic_sh:
#           type: File
#           inputBinding:
#             position: 1
#         download_databases_done:
#           type: string
#       outputs:
#         indexing_Cosmic_done:
#           type: string


#     in:
#       indexing_Cosmic_sh: indexing_Cosmic_sh
#       download_databases_done: download_databases/download_databases_done
#     out: [indexing_Cosmic_done]

# #-----------------------------------#

# #This step adds permissions in get_pass_variants
#   add_permissions:
#     run:
#       class: CommandLineTool
#       baseCommand: bash
#       inputs:
#         add_permissions_sh:
#           type: File
#           inputBinding:
#             position: 1
#         sleep_after_copying_done:
#           type: string
#       outputs:
#         add_permissions_done:
#           type: string


#     in:
#       add_permissions_sh: add_permissions_sh
#       sleep_after_copying_done: sleep_after_copying/sleep_after_copying_done
#     out: [add_permissions_done]









