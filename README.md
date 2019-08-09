# Nextflow pipeline for GATK best practices, SNP calling

## Installation

This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub rbpisupati/nf-haplocaller
The pipeline runs with singularity container (based on environment.yml file included with package).

```bash
git clone https://github.com/Gregor-Mendel-Institute/nf-haplocaller.git
```

## Running the pipeline

```bash
nextflow run nf-haplocaller/main.nf --reads "*bam" --fasta "TAIR10_wholeGenome.fasta" --outdir output_folder -profile mendel/cbe
```

Add `--notrim false` to include trimming of the reads. Also you can change the trimming parameters. You can to provide a VCF file under `--known-sites` parameter to perform BQSR step for all the samples.

Also you can run the pipeline in either PBS (`mendel`) or slurm (`cbe`) batch environment. You can also run the pipeline on your local system providing `-profile local`.

## Credits

- Rahul Pisupati (rahul.pisupati[at]gmi.oeaw.ac.at)
- Fernando Rabanal (fernando.rabanal@tuebingen.mpg.de)

## Citation
Cite the paper below if you use this pipeline.
Pisupati, R. *et al.*. Verification of *Arabidopsis* stock collections using SNPmatch, a tool for genotyping high-plexed samples.  *Nature Scientific Data*  **4**, 170184 (2017).
[doi:10.1038/sdata.2017.184](https://www.nature.com/articles/sdata2017184)
