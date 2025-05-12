# RNAseq-processing workflow

## 1. Project Overview

This repository contains a comprehensive RNA sequencing (RNA-seq) data processing workflow implemented in Nextflow. The workflow provides an end-to-end solution for analyzing RNA-seq data, from raw sequencing reads to meaningful biological insights.

**Purpose**

The RNAseq-processing workflow standardizes and automates the complex process of RNA-seq analysis, ensuring reproducibility and scalability while following best practices in the field. It handles the computationally intensive tasks of read alignment, quality control, quantification, and differential expression analysis in a streamlined manner.

**Key Features** 

* Modular design: Each analysis step is implemented as a separate process that can be executed independently or as part of the complete workflow
* Extensive quality control: Incorporates FastQC and MultiQC for raw data assessment and RNA-SeQC for post-alignment validation
* Flexible read processing: Optional adapter trimming using Trimmomatic before alignment
* High-performance alignment: Utilizes STAR for efficient and accurate read mapping to reference genomes
* Robust quantification: Implements RSEM for transcript abundance estimation
* Statistical analysis: Integrates with DESeq2 for differential expression analysis
* HPC compatibility: Optimized for execution on high-performance computing clusters using SLURM scheduler
 
This workflow is designed for researchers and bioinformaticians working with RNA-seq data who need a reliable, efficient, and standardized approach to process and analyze their sequencing results, regardless of sample size or computational environment.
