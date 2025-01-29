#!/bin/bash

# This script submits a job to the HPC job scheduler to build HISAT2 index.
#
# projdir is the project directory
# outputdir is the output directory
# ht2index is the HISAT2 index base

projdir=/Shared/pezzulolab/inlysate
outputdir=20200217_HISAT2_run
ht2index=/Shared/pezzulolab/inlysate/raw/ref/HSapiens_GRCh38.p13_RefSeq_hisat2ind

if [ ! -d $projdir/exp/$outputdir ]; then
  mkdir $projdir/exp/$outputdir
fi

sample_dir=(`ls -d $projdir/raw/*/`)
l_sample_dir=${#sample_dir[*]}

for i in $(seq 0 $((l_sample_dir-1)))
  do 
    fastq_TF=$((`ls ${sample_dir[i]} | grep 'fastq' | wc -l`))
    if [ ! $fastq_TF = 0 ]; then
      sample_name_i=(`echo "${sample_dir[i]}" | awk '{split($1,temp,"raw/"); print temp[2]}' | awk '{split($1,temp2,"/"); print temp2[1]}'`)
      if [ ! -d $projdir/exp/$outputdir/HISAT2_$sample_name_i ]; then
        qsub -q CCOM -N HS2_$sample_name_i -wd $projdir/exp/$outputdir -e $projdir/exp/$outputdir -o $projdir/exp/$outputdir -pe smp 16 $projdir/bin/run_hisat2.job $projdir $ht2index $sample_name_i
      fi
    fi
  done


