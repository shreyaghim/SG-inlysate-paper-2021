#!/bin/bash

# This script submits a job to the HPC job scheduler to build HISAT2 index.
#
# projdir is the project directory
# fafile is the genome fasta file

projdir=/Shared/pezzulolab/inlysate
fafile=HSapiens_GRCh38.p13_RefSeq.fa
ht2index=HSapiens_GRCh38.p13_RefSeq_hisat2ind

qsub -q CCOM -N build_index -wd $projdir/raw/ref -e $projdir/raw/ref -o $projdir/raw/ref -pe smp 16 $projdir/bin/build_hisat2_index.job $fafile $ht2index

