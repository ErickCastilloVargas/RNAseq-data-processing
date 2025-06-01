# RNA-seq Processing Pipeline

> **Project status**: Work in Progress 🚧
>
> This version is not yet complete. Some features are still pending implementation and it may contain bugs.

This repository provides a Nextflow-based workflow for processing RNA-seq data on HPC clusters (Slurm), from raw FASTQ files through to gene-level count matrices. The pipeline includes comprehensive quality control (QC) at each step.

### Features

1. **Initial QC (Raw Data)**
   - Inspect raw FASTQ files for sequencing quality.
   - Determine if adapter trimming is necessary.
   - Tools: FastQC, MultiQC.

2. **Data Processing Workflow**
   1. **Adapter Trimming (optional)**  
      - Tool: Trimmomatic  
      - Run FastQC and MultiQC again on trimmed reads to verify improvement.
   2. **Sequence Alignment**  
      - Tool: STAR  
      - Align reads to the reference genome.
   3. **BAM Sorting**  
      - Tool: SAMtools  
      - Sort aligned BAM files by genomic coordinates.
   4. **Alignment QC**  
      - Tool: RNASeQC  
      - Evaluate alignment quality (e.g., mapping rates, coverage).
   5. **Transcript Quantification**  
      - Tool: RSEM  
      - Generate gene-level abundance estimates.

## Prerequisites

- **Nextflow** (≥20.10)  
- **Java** (≥1.8)  
- **Slurm** on your HPC cluster  
- Reference genome FASTA and GTF (collapsed GTF if using RSEM)  
- Adapter sequence file (if trimming is required) 

## Usage

1. **Organize Your FASTQ Files**  
   Place all FASTQ files for your samples into a single directory (hereafter, `sampledir/`).

2. **Run Initial QC Workflow**  
   Inspect raw reads and determine whether adapter trimming is needed:
   ```bash
   nextflow run nf/workflows/raw_data_QC.nf \
     --sampleDir sampledir \
     --outDir <OUTDIR>
   ```
   * `--sampleDir`: Path to the directory containing FASTQ files
   * `--outDir`: Output directory for QC reports
   * `--adapters`: FASTA file of adapter sequences

3. Run Data Processing Workflow
  After evaluating QC reports, run the full processing pipeline (include `--adapters` if trimming is required):
  ```bash
  nextflow run nf/workflows/data_processing.nf \
    --genome        <PATH_TO_GENOME_FASTA> \
    --gtf           <PATH_TO_GTF> \
    --gtfCollapsed  <PATH_TO_COLLAPSED_GTF> \
    --sampleDir     sampledir \
    --outDir        <OUTDIR> \
    --adapters      <PATH_TO_ADAPTER_FILE>   # (optional)
  ```
  * `--genome`: Reference genome FASTA
  * `--gtf`: Annotation GTF file
  * `--gtfCollapsed`: Collapsed GTF (required by RSEM)
  * `--sampleDir`: Directory with raw or trimmed FASTQ files
  * `--outDir`: Base output directory for all pipeline results
  * `--adapters`: (Optional) FASTA of adapter sequences for trimming

## Output Structure

All results are organized under `<OUTDIR>/`. Each subfolder contains the relevant files:

* **`trimmomatic_adapter_trimming/`**  
  Trimmed FASTQ files (if trimming was enabled).

* **`fastQC_reports_post_trimming/`**  
  FastQC reports for trimmed FASTQ files.

* **`multiQC_reports_post_trimming/`**  
  Aggregated MultiQC report summarizing all post-trimming QC.

* **`indexes/`**  
  STAR and RSEM indices generated from the reference genome and GTF.

* **`STAR_alignment/`**  
  Aligned BAM files (genome and transcriptome alignments) for each sample.

* **`SAMtools_sort_by_coordinates/`**  
  Coordinate-sorted BAM files for each sample.

* **`QC_post_alignment/`**  
  RNASeQC reports evaluating the quality of each aligned BAM.

* **`RSEM_transcript_quantification/`**  
  Gene-level quantification files (e.g., `.gene.results`).

* **`logs/`**  
  Log files and reports from every pipeline step.

## Notes & Tips

* **Resource Allocation**:
  If needed, adjust CPU and memory parameters in the Nextflow configuration (`nextflow.config`) file based on your HPC cluster specifications.

* **Resume Capability**:
  If the pipeline stops unexpectedly, rerun the same Nextflow command to resume from the last successful step using the `-resume` flag in the command.