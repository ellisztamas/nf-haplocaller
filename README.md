# Nextflow pipeline for calling SNPs based on GATK best practices

## Table of contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Running the pipeline](#running-the-pipeline)
4. [Dependencies](#dependencies)
5. [Credits](#credits)
6. [Citation](#citation)

## Overview

Here is a high-level summary of the pipeline:

* Convert BAM to fastq.
* Fastqc on the input files
* Trim Galore! on the input files to trim reads and repeat quality control on trimmed reads.
* Align reads to a reference genome with BWA-MEM
* Sort, index and return statistics with samtools
* Remove duplicate reads with picard
* qualimap on deduplicated reads
* Not run: GATK to realign the reads at the positions where there are indels (this was deprecated in GATK 4)
* Base recalibration with GATK tools BaseRecalibrator, ApplyBQSR and AnalyzeCovariates
* SNP calls with gatk HaplotypeCaller.
* Multiqc to summarise the various QC checks.

See `main.nf` for full details.

## Installation

This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub rbpisupati/nf-haplocaller
The pipeline runs with singularity container (based on environment.yml file included with package).

```bash
git clone https://github.com/Gregor-Mendel-Institute/nf-haplocaller.git
```

## Running the pipeline

### Basic usage

Assuming you have:

* cloned the repo into directory `library`
* a set of BAM files to process in directory `data`
* a FASTA file gibing a reference genome to align to in directory `data`
* a valid, active installation of NextFlow

then a minimal command to run the pipeline is:

```bash
nextflow run library/nf-haplocaller/main.nf \
    --reads "data/*bam" \
    --fasta "data/TAIR10_wholeGenome.fasta" \
    --outdir output_folder \
    -profile cbe
```

This will take a long time, so it is recommended to run this in a detatchable window, such as tmux.

### Options

Here is a full list of arguments and options

#### Input

* `--reads`: Path to input files. This will usually include a wildcard to include
all files matching a pattern, and be enclosed in double quotes ("").
* `--fasta`: Optional path to a reference fasta file to align reads to.
* `--file_ext`: File type of the input files. Options are "bam", 'fastq' and 'aligned_bam'.
* `--singleEnd`: Flag for whether data are single- or paired end. Defaults to false.
* `-profile`: Give a nextflow profile to allow the pipeline to talk to the job scheduling system on your machine. Valid inputs are:
    * `mendel` for PBS systems
    * `cbe` for SLURM systems
    * `singularity`
    * `local` to run on a local machine

#### Processing

* `--saveTrimmed`: If true, keep trimmed data. Defaults to false.
* `--notrim`: If true, skip trimming reads. Defaults to false
<!-- params.illumina = true I couldn't tell what this does - TJE-->
* `--clip_r1` Integer number of bases to trim from the 5\` end of read 1.
* `--clip_r2` Integer number of bases to trim from the 5\` end of read 2.
* `--three_prime_clip_r1` Integer number of bases to trim from the 3\` end of read 1.
* `--three_prime_clip_r2` Integer number of bases to trim from the 3\` end of read 2.

#### Output

* `--project`: Project name
* `--outdir`: Path to directory for the results.
* `--cohort`: Optional. Specify a group of samples to lump together into a single output file.
<!-- params.`run_name` = false Not sure what this is - TJE-->
<!-- params.clusterOptions = false Not sure what this is: TJE -->
* `--email`: Optional email address to contact when the pipeline finishes.
* `--plaintext_email` = If true, send the notification email in plain text.
* `-w`: Path to working directory. Defaults to the current working directory. Note that `w` is preceded by only one hyphen.

## Credits

- Rahul Pisupati (rahul.pisupati[at]gmi.oeaw.ac.at)
- Fernando Rabanal (fernando.rabanal@tuebingen.mpg.de)

## Citation

Please cite the paper below if you use this pipeline.

> Pisupati, R. *et al.*. Verification of *Arabidopsis* stock collections using SNPmatch, a tool for genotyping high-plexed samples.  *Nature Scientific Data*  **4**, 170184 (2017).
[doi:10.1038/sdata.2017.184](https://www.nature.com/articles/sdata2017184)
