#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --output=logs/salmon_%A_%a.out
#SBATCH --error=logs/salmon_%A_%a.err
#SBATCH --array=0-62
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=zahra.shahitat@kcl.ac.uk
#SBATCH --mail-type=END,FAIL

#Activate conda environment 
source /users/k21017883/miniconda3/etc/profile.d/conda.sh
conda activate salmon_env

#Define directories
IN_DIR=/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/trimmed
OUT_DIR=/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/salmon_quant
INDEX=/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/salmon_index

mkdir -p $OUT_DIR
 

# Get list of paired-end FASTQ files
R1_FILES=($IN_DIR/*_R1*.fastq.gz)
R2_FILES=($IN_DIR/*_R2*.fastq.gz)

# Select file for this array job
R1_FILE=${R1_FILES[$SLURM_ARRAY_TASK_ID]}
R2_FILE=${R2_FILES[$SLURM_ARRAY_TASK_ID]}


# Get all R1 and R2 files into Bash arrays
R1_FILES=(/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/trimmed/*_R1_001_val_1.fq.gz)
R2_FILES=(/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/trimmed/*_R2_001_val_2.fq.gz)

# Use SLURM array index to pick sample
R1_FILE=${R1_FILES[$SLURM_ARRAY_TASK_ID]}
R2_FILE=${R2_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract sample name from filename (adjust to match your pattern)
SAMPLE_NAME=$(basename "$R1_FILE" | cut -d "_" -f1,2)

# Run Salmon (only if output doesn't already exist)
if [ ! -d "$OUT_DIR/$SAMPLE_NAME" ]; then
  salmon quant \
    -i $INDEX \
    -l A \
    -1 "$R1_FILE" \
    -2 "$R2_FILE" \
    --validateMappings \
    -p 4 \
    -o "$OUT_DIR/$SAMPLE_NAME"
else
  echo "Skipping $SAMPLE_NAME – already exists."
fi

