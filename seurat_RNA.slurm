#!/bin/bash

##General options
#SBATCH --account=nn4605k
#SBATCH --job-name="Seurat_RNA"
#SBATCH --partition=bigmem
#SBATCH --time=0-10:00:00


##CPU
#SBATCH --mem-per-cpu=100G
#SBATCH --ntasks=1

set -e
set -u


##Load modules
module load R-bundle-CRAN/2023.12-foss-2023a 

#Run the script
cd /cluster/projects/nn4605k/peter/from_junbai/to_peter_data/scripts
Rscript Seurat_RNA.R

#Unload the modules
module unload R-bundle-CRAN/2023.12-foss-2023a

#End of the script

exit 0


