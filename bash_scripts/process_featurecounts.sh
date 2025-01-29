
#!/bin/bash

module load R

# This script will process the results for runs of featureCounts

projectdir=/Shared/pezzulolab/inlysate
outputdir=$projectdir/exp/20210411_featurecounts_run
fnames=(`ls $outputdir/*/*summary`)
samp_names=(`ls $outputdir/ | awk '{split($1,b,"featurecounts_"); print b[2]}'`)

echo "Sample" > $outputdir/var_names
awk 'NR > 1 && NR < 14{print $1}' ${fnames[0]} >> $outputdir/var_names
cat $outputdir/var_names > $outputdir/all_run_info.txt

n=${#fnames[*]}
for i in $(seq 0 $((n-1)));
  do
    echo ${samp_names[i]} > $outputdir/temp1;
    awk 'NR > 1 && NR < 14{print $2}' ${fnames[i]} >> $outputdir/temp1;
    cat $outputdir/all_run_info.txt > $outputdir/temp2;
    paste $outputdir/temp2 $outputdir/temp1 > $outputdir/all_run_info.txt;
  done

rm $outputdir/temp1 $outputdir/temp2 $outputdir/var_names 

Rscript --no-save --no-restore $projectdir/bin/process_featurecounts.R > $outputdir/temp.Rout $outputdir

rm $outputdir/temp.Rout


