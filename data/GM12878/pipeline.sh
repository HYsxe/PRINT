#!/bin/bash
#SBATCH -n 16                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test   # Partition to submit to
#SBATCH --mem=80G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o run_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e run_%j.err  # File to which STDERR will be written, %j inserts jobid

# Specify reference genome
genome="hg19"

# Get input bam file
bam=GM.trial52.atac.hg19.rmdup.cutoff.bam

# Arguments
fname=$bam
is_bap=false
cores=16

# Load required modules
module load R/4.1.0-fasrc01
module load macs2/2.1.1.20160309-fasrc02
module load bedtools2/2.26.0-fasrc01
module load intel/2017c impi/2017.4.239 SAMtools/1.9

# Peak Calling
if [ ! -f peaks.bed ]; then

  # Call peaks using MACS2
  mkdir peakCalling
  experiment=`echo $fname | sed 's/.bam//g'`
  macs2 callpeak \
    -t $bam \
    -f BAMPE \
    -n $experiment \
    --outdir peakCalling \
    --keep-dup all \
    --nolambda --nomodel \
    --call-summits 
  
  # Filter and resize peaks
  singularity exec --cleanenv --env R_LIBS_USER=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/apps/R_4.1.0 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_13.sif Rscript ../../code/filterPeaks.R $genome
fi

# If input bam is generated using BAP, the read names will not contain cell barcodes
# Therefore, we rename the reads by adding the cell barcodes to read names
if [ "$is_bap" = true  ]; then

  echo "Input is BAP data"
  index_bam=`echo $fname | sed 's/.bam/.bam.bai/g'`
  renamed_bam=`echo $fname | sed 's/.bam/.renamed.bam/g'`
  
  if [ ! -f $index_bam ]; then
    echo "Indexing the input bam file"
    samtools index $bam -@16
  fi

  if [ ! -f $renamed_bam ]; then
    echo "Renaming reads in the input bam file"
    python3 ../../code/renameRead.py $bam DB $renamed_bam
  fi
  
  fname=$renamed_bam
  
fi

# Convert bam file to fragments files
sorted_bam=`echo $fname | sed 's/.bam/.name_st.bam/g'`
frags_name=`echo $fname | sed 's/.bam/.frags.gz/g'`
echo $frags_name

if [ ! -f $sorted_bam ]; then
  echo "Name sorting bam file"
  samtools sort -n -@ $cores -o $sorted_bam $fname
fi

if [ ! -f all.frags.tsv.gz ]; then
  echo "Converting bam file to fragments file"
  bedtools bamtobed -i $sorted_bam -bedpe | \
  sed 's/_/\t/g' | \
  awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6-5,$8}else if($10=="-"){print $1,$2-5,$6+4,$8}}' |\
  sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n -k3,3n | \
  uniq -c | \
  awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | \
  gzip > $frags_name

  # Merge fragments files
  zcat *frags.gz | gzip > all.frags.tsv.gz
fi

if [ ! -f all.frags.filt.tsv.gz ]; then
  echo "Filtering out low quality cells"
  singularity exec --cleanenv --env R_LIBS_USER=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/apps/R_4.1.0 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_13.sif Rscript ../../code/filterFrags.R
fi

if [ $genome = "hg19" ]; then
  echo "Lifting over from hg19 to hg38"
  singularity exec --cleanenv --env R_LIBS_USER=/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/apps/R_4.1.0 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_13.sif Rscript ../../code/liftOverFrags.R
fi
