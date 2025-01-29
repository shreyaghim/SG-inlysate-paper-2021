
#!/bin/bash

# This script will run the FastQC software on all of the raw sequencing data sets 
# (in gz compressed fastq format)

projectdir=/Shared/pezzulolab/inlysate
outputdir=$projectdir/exp/20200214_FastQC_results
filesuffix=L00
filenames=(`ls $projectdir/raw/*/* | grep "$filesuffix"`)


l_filenames=${#filenames[*]}

if [ ! -d $outputdir ]; then
  mkdir $outputdir
fi

for i in $(seq 0 $((l_filenames-1)))
  do
    qsub -q all.q -N fqc_job_$i -wd $outputdir -e $outputdir -o $outputdir -pe smp 4 $projectdir/bin/run_fastqc.job $outputdir ${filenames[i]}
  done

