#!/bin/bash

# This script submits all jobs to the HPC job scheduler.
#
# projdir is the project directory
# inputdir is the input directory 
# outputdir is the output directory 

projdir=/Shared/pezzulolab/inlysate
inputdir=$projdir/exp/20200217_HISAT2_run
#outputdir=$projdir/exp/20200219_featurecounts_run
outputdir=$projdir/exp/20210411_featurecounts_run
usrname=(`whoami`)

if [ ! -d $outputdir ]; then
  mkdir $outputdir
fi

sampledir=(`ls -d $inputdir/*/`)
l_sample_dir=${#sampledir[*]}

for i in $(seq 0 $((l_sample_dir-1)))
  do 
    fastq_TF=$((`ls ${sampledir[i]} | grep '.sam' | wc -l`))
    if [ ! $fastq_TF = 0 ]; then
      sample_name_i=(`ls ${sampledir[i]} | grep ".sam" | awk '{split($1,temp,".sam"); print temp[1]}'`)
      if [ ! -d $outputdir/featurecounts_$sample_name_i ]; then
        qsub -q CCOM,UI -N fc_job_$sample_name_i -wd /localscratch/Users/$usrname -e /localscratch/Users/$usrname -o /localscratch/Users/$usrname -pe smp 16 $projdir/bin/run_featurecounts.job $projdir $inputdir $outputdir $sample_name_i 
      fi
    fi
  done

