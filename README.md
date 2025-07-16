RNA-seq Analysis Workflow
This repository offers SLURM batch scripts and documentation for executing a comprehensive RNA-seq analysis pipeline. It includes steps for data preprocessing, genome alignment, transcript-level quantification, and output refinement for downstream analysis tasks.

Workflow Summary:
The pipeline comprises these core stages:

1. Quality Control and Trimming
Script: 1_trim_rna.slurm
Tool: Trim Galore

- Conducts paired-end read trimming

- Applies a Phred quality score threshold of 30

- Outputs saved in trimmed/ directory

- Uses SLURM array jobs to handle multiple samples efficiently

2. STAR Genome Index Generation
Script: 3_STAR_genome_index.sh
Tool: STAR

- Constructs a genome index using GRCh38 and GENCODE annotations

- Set with --sjdbOverhang 99 for 100bp sequencing reads

- Output path: STAR/

3. Read Alignment with STAR
- Script: 3_STAR_alignment.sh
- Tool: STAR

- Aligns quality-trimmed reads against the reference genome

- Produces both genome-aligned and transcriptome BAM files

- Includes --quantMode TranscriptomeSAM GeneCounts for quantification

4. Merging  BAM Files
Script: merge_lanes.sh
Tool: Samtools

- Integrates sorted BAM files from different flowcell lanes

- Optimises for consistent read alignment

- Removes temporary files to conserve disk space

5. Transcript Quantification with Salmon
Script: run_salmon.sh
Tool: Salmon

- Uses quasi-mapping against GRCh38 transcriptome index

- Applies --validateMappings and corrects for GC bias

- Results stored in salmon_quant/ folder

Required Inputs:
- files.txt: Sample identifiers for SLURM array execution

- Reference genome (.fa) and annotation file (.gtf)

- Raw FASTQ data placed in raw_fastq/

Output Structure:
- trimmed/: Contains cleaned and trimmed reads

- STAR/: Location of STAR genome index files

- aligned/: Output BAMs from STAR alignment

- salmon_quant/: Salmon quantification results

Tool References
 The following tools used in this pipeline:

Trim Galore version 0.6.10


STAR version 2.7.10.b

Salmon version 1.10.3 

Samtools version 1.9


This pipeline was developed as part of an MSc project in Genomic Medicine at King's College London by Zahra Shahitat.
