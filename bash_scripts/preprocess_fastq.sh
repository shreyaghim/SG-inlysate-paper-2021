
#!/bin/bash

# This script will submit the raw sequencing files to be pre-processed

projectdir=/Shared/pezzulolab/inlysate
outputdir=$projectdir/raw
filesuffix=R1
infilenames=(`ls $projectdir/raw/*/* | grep "$filesuffix" | grep -v "trim" | awk '{gsub(".fastq.gz","");print $0}'`)

l_infilenames=${#infilenames[*]}

for i in $(seq 0 $((l_infilenames-1)))
  do
    qsub -q all.q -pe smp 4 -e $outputdir -o $outputdir -wd $outputdir $projectdir/bin/preprocess_fastq.job $projectdir ${infilenames[i]}
  done

