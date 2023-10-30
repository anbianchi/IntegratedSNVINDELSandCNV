#!/bin/bash
#$ -S /bin/bash
#$ -pe mpi 50
#$ -cwd
#$ -o /NFSHOME/abianchi/log/std_$JOB_ID.out
#$ -e /NFSHOME/abianchi/log/err_$JOB_ID.out
#$ -l h=compute-0-6.local

source ~/.bashrc
conda activate approach

#snakemake --use-conda --cores 80
snakemake --use-conda --cores 80 --rerun-incomplete
#snakemake --use-conda --cores 40 --delete-temp-output --printshellcmds --summary
#snakemake --use-conda --cores 40 --dry-run --debug-dag
#snakemake --use-conda --cores 80 --rerun-incomplete --latency-wait 15
#snakemake --use-conda --cores 80 --rerun-incomplete --unlock


