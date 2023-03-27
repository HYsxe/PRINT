#!/bin/bash
#SBATCH -n 16                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p bigmem   # Partition to submit to
#SBATCH --mem=200G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o run_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e run_%j.err  # File to which STDERR will be written, %j inserts jobid

# Arguments
genome="hg38"
cores=16

cp ./rawData/T133.hepG2.ATAC.hg38.fragments.tsv all.frags.tsv
gzip all.frags.tsv

# Load required modules
module load R/4.1.0-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.1.0:$R_LIBS_USER
module load macs2/2.1.1.20160309-fasrc02

# Peak Calling
if [ ! -f peaks.bed ]; then

  macs2 callpeak --nomodel -t ./rawData/T133.hepG2.ATAC.hg38.fragments.tsv \
        --outdir ./peakCalling -n HepG2 -f BEDPE \
        --nolambda --keep-dup all --call-summits
  
  # Filter and resize peaks
  singularity exec --cleanenv --env R_LIBS_USER=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/apps/R_4.1.0 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_13.sif Rscript ../../code/filterPeaks.R $genome
  
fi

# Remove fragments from low quality cells
if [ ! -f all.frags.filt.gz ]; then
  echo "Filtering out low quality cells"
  singularity exec --cleanenv --env R_LIBS_USER=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/apps/R_4.1.0 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_13.sif Rscript ../../code/filterFrags.R
fi
