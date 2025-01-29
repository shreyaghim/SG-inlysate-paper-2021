
#!/bin/bash

# This script submits a job that merges bam files and indexes the merged file

projdir=/Shared/pezzulolab/inlysate
samdir=/Shared/pezzulolab/inlysate/exp/20200217_HISAT2_run
outdir=/Shared/pezzulolab/inlysate/exp/20200217_HISAT2_run
qsub -q CCOM -N mrg_bam -pe smp 56 -wd $outdir -e $outdir -o $outdir $projdir/bin/merge_bam.job $samdir $outdir

