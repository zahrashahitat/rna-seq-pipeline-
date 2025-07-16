#!/bin/bash
#SBATCH --job-name=trim_rna
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --array=1-63
#SBATCH --output=logs/trim_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zahra.shahitat@kcl.ac.uk

# Activate conda environment
source /users/k21017883/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

# Add TrimGalore to PATH (manually installed version)
export PATH=/users/k21017883/TrimGalore-0.6.10:$PATH

# Load environment variables
source /scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/paths.sh

# Move to FASTQ directory
cd "$FASTQ_DIR" || { echo "FASTQ directory $FASTQ_DIR not found!"; exit 1; }

# Get sample name
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../files.txt)

if [[ -z "$sample" ]]; then
  echo "Sample name is empty for task ID $SLURM_ARRAY_TASK_ID!"
  exit 1
fi

# Define input files
r1_file="${sample}_R1_001.fastq.gz"
r2_file="${sample}_R2_001.fastq.gz"

# Check if input files exist
if [[ ! -f "$r1_file" ]]; then
  echo "Error: $r1_file not found!"
  exit 1
fi

if [[ ! -f "$r2_file" ]]; then
  echo "Error: $r2_file not found!"
  exit 1
fi

# Create output directory
mkdir -p ../trimmed

# Run Trim Galore
echo "Running trim_galore for sample $sample..."
trim_galore --paired "$r1_file" "$r2_file" -o ../trimmed

# Check success
if [[ $? -ne 0 ]]; then
  echo "Trim_galore failed for sample $sample!"
  exit 1
else
  echo "Trimming completed successfully for sample $sample."
fi

x
