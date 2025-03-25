#!bin/bash

singularity exec --bind &PATH_BIN_HPC/bin/resource/cosmic:/home/dbs \
&PATH_BIN_HPC/bin/def/sv_cons_general_utilities/sv_cons_general_utilities.sif \
tabix -f &PATH_BIN_HPC/bin/resource/cosmic/CosmicCodingMuts.vcf.gz