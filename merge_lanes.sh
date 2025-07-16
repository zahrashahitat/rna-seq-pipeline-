#!/bin/bash
#SBATCH --job-name=merge_lanes
#SBATCH --output=logs/merge_%A_%a.out
#SBATCH --error=logs/merge_%A_%a.err
#SBATCH --array=1-16
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --mail-user=zahra.shahitat@kcl.ac.uk
#SBATCH --mail-type=END,FAIL

# Activate conda
source /users/k21017883/miniconda3/etc/profile.d/conda.sh
conda activate samtools_fix

# Set paths
BAM_DIR="/scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/aligned"
cd $BAM_DIR

# Get sample
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/prj/bcn_pd_pesticides/SHSY5Y/RNA/files.txt)
echo "Merging lanes for: $sample"

# Sort both BAMs before merging
samtools sort -@ 4 -o ${sample}_L001_sorted.bam ${sample}_L001_Aligned.toTranscriptome.out.bam
samtools sort -@ 4 -o ${sample}_L002_sorted.bam ${sample}_L002_Aligned.toTranscriptome.out.bam

# Merge the sorted BAMs
samtools merge -@ 4 ${sample}_merged_Aligned.toTranscriptome.out.bam \
    ${sample}_L001_sorted.bam ${sample}_L002_sorted.bam

# Clean up
rm ${sample}_L001_sorted.bam ${sample}_L002_sorted.bam





