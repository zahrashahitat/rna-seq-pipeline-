#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=logs/star_index.out
#SBATCH --error=logs/star_index.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --mail-user=zahra.shahitat@kcl.ac.uk
#SBATCH --mail-type=END,FAIL

# 1. Load Conda and activate your RNA-seq environment
source /users/k21017883/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

# 2. Set input and output directories
GENOME_DIR="/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/genome"
STAR_DIR="$GENOME_DIR/STAR"

# 3. Create output directory if it doesn't exist
mkdir -p $STAR_DIR

# 4. Run STAR genome indexing
STAR --runThreadN $SLURM_CPUS_PER_TASK \
     --runMode genomeGenerate \
     --genomeDir $STAR_DIR \
     --genomeFastaFiles $GENOME_DIR/genome.fa \
     --sjdbGTFfile $GENOME_DIR/annotation.gtf \
     --sjdbOverhang 99

