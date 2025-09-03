# RNA-seq Data Processing Pipeline

This repository provides a Nextflow-based workflow for processing RNA-seq pair-end data on HPC clusters (Slurm), from raw FASTQ files to gene-level count matrices. The pipeline includes comprehensive quality control (QC) before and after the alignments.

### Features

1. **Initial QC (Raw Data)**
   - Inspect raw FASTQ files for sequencing quality.
   - Determine if adapter trimming is necessary.
   - Tools: FastQC, MultiQC.

2. **Data Processing Workflow**
  1. **Infer strandness**
      - Tool: Salmon
  1. **Adapter Trimming (optional)**  
      - Tool: Trimmomatic 
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

- **Nextflow** (≥20.10).
- **Java** (≥1.8).  
- **Slurm** on your HPC cluster.  
- Reference data: genome FASTA, annotation GTF and transcriptome FASTA files.
- Adapter sequence file (if trimming is required).

> You can download the reference data from GENCODE. It is important to ensure that the reference genome, GTF annotation, and transcriptome files all come from the same genome version.

## Usage

1. **Organize Your FASTQ Files**  
  Place all FASTQ files for your samples into a single directory (hereafter, `sampledir/`).

> ⚠️ **Important**: Make sure that the files follow the naming convention SRR12345678_1.fastq and SRR12345678_2.fastq. Files may also be provided in compressed .gz format.

2. **Run Initial QC Workflow**  
  Inspect raw reads and determine whether adapter trimming is needed:
  ```bash
  nextflow run nf/workflows/raw_data_QC.nf \
     --sampleDir sampledir \
     --outDir <OUTDIR> \
     -with-report
  ```
  * `--sampleDir`: Path to the directory containing the samples FASTQ files
  * `--outDir`: Output directory for QC reports
  * `--adapters`: FASTA file of adapter sequences
  * `-with-report`: Indicate Nextflow to create an HTML execution report. It is a single document which includes many usefulmetrics about the workflow execution. The report is organised in the three main sections: Summary, Resources and Tasks. 

3. Run Data Processing Workflow
  After evaluating QC reports, run the full processing pipeline (include `--adapters` if trimming is required):
  ```bash
  nextflow run nf/workflows/data_processing.nf \
    --genome <PATH_TO_GENOME_FASTA> \
    --gtf <PATH_TO_GTF> \
    --transcriptome  <PATH_TO_TRANSCRIPTOME> \
    --sampleDir sampledir \
    --outDir <OUTDIR> \
    --adapters <PATH_TO_ADAPTER_FILE> \
    -with-report
  ```
  * `--genome`: Reference genome FASTA file
  * `--gtf`: Annotation GTF file
  * `--transcriptome`: Transcriptome FASTA file
  * `--sampleDir`: Directory with raw or trimmed FASTQ files
  * `--outDir`: Base output directory for all workflow results
  * `--adapters`: (Optional) FASTA of adapter sequences for trimming

> If you want to keep the alignment files, you can use the `--getAlignments` option, which will save them in CRAM format in the `outDIR/Alignment_files` folder, created automatically.

> If the workflow fails, after fixing the error you can resume it from the last successful step by adding the `-resume` flag to the nextflow command above.

## Output Structure

All results are organized under `<OUTDIR>/`. The directory structure is as follows:

* **`Alignment_files/`** (optional): Contains the alignment files for each sample, stored in CRAM format. These files can be used for downstream analyses or re-processed if needed.

* **`Gene_level_quant/`**: Contains the gene-level quantification tables (e.g., raw counts, normalized counts) generated after alignment and read summarization.

* **`Reports/`**: Includes the QC reports generated with MultiQC. This folder stores both raw data QC (before trimming) and processed data QC (after trimming). The latter includes not only the FastQC metrics of the trimmed data, but also the strandness inference and the RNA-SeQC metrics from the alignment.

## Notes & Tips

* **Resource Allocation**: If necessary, adjust the CPU and memory parameters in the `nextflow.config` file according to your HPC cluster specifications and analysis requirements.

* When launching the Nextflow commands, make sure you are in the pipeline directory.

* After a successful run, if you are sure everything is correct, you can use the `nextflow clean -f` command to erase the intermediate files and free up disk space.
