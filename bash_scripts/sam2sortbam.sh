
#!/bin/bash

# This script submits one job for every HISAT2 sam file
# to the cluster to convert sam to bam and sort.

projdir=/Shared/pezzulolab/inlysate
HS2dir=/Shared/pezzulolab/inlysate/exp/20200217_HISAT2_run

cd $HS2dir

filedir=(`ls -d */ | awk '{split($1,b,"/"); print b[1]}'`)
l_filedir=${#filedir[*]}

for i in $(seq 0 $((l_filedir-1)))
  do
    if [ ! -f $HS2dir/${filedir[i]}/${filedir[i]}"_sorted.bam" ]; then
      qsub -N sam2sortbam_${filedir[i]} -q all.q -pe smp 4 -wd $HS2dir/${filedir[i]} -e $HS2dir/${filedir[i]} -o $HS2dir/${filedir[i]} $projdir/bin/sam2sortbam.job $HS2dir ${filedir[i]}
    fi
  done

