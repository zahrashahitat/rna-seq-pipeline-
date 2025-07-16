#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --output=logs/STARalign_%A_%a.log
#SBATCH --error=logs/STARalign_%A_%a.err
#SBATCH --array=1-63
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# Load Conda and activate your STAR environment
source /users/k21017883/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

# Set paths
GENOME_DIR="/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/genome/STAR"
OUTPUT_DIR="/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/aligned"
FILE_DIR="/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/trimmed"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Read sample name from files.txt
cd $FILE_DIR
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" files.txt)

# Run STAR alignment
STAR \
  --runThreadN $SLURM_CPUS_PER_TASK \
  --genomeDir $GENOME_DIR \
  --readFilesCommand zcat \
  --readFilesIn ${file}_R1_001_val_1.fq.gz ${file}_R2_001_val_2.fq.gz \
  --outFileNamePrefix $OUTPUT_DIR/${file}_ \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts

