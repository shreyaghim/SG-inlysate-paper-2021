
#!/bin/bash

fastqc_dir=/Shared/pezzulolab/inlysate/exp/20200214_FastQC_results
outdir=/Shared/pezzulolab/inlysate/exp/20200217_MultiQC_results
multiqc -o $outdir $fastqc_dir
