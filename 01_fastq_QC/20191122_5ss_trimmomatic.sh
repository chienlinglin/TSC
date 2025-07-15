#!usr/bin/sh
projectPath=/home2/cllin/run3_run4/
readPath=${projectPath}5ss_raw/
outputPath=${projectPath}5ss_qc/
output_unPath=${outputPath}unpair/
qcoutputPath=${projectPath}5ss_qc/
recordPath=${projectPath}5ss_qc/summary/
listRaw=(${readPath}*R1_001.fastq.gz)

mkdir ${output_unPath} ${recordPath}


for R1 in ${listRaw[@]}
do 
  R2=${readPath}$(basename ${R1%%R1_001.fastq.gz})R2_001.fastq.gz
  trimmomatic \
           PE \
          -threads 16 \
          -summary ${recordPath}$(basename ${R1%%.fastq.gz}).summary \
          $R1 \
          $R2 \
          ${outputPath}$(basename ${R1%%R1_001.fastq.gz})R1.QC.fastq.gz \
          ${output_unPath}$(basename ${R1%%R1_001.fastq.gz})R1.unpair.QC.fastq.gz \
          ${outputPath}$(basename ${R2%%R2_001.fastq.gz})R2.QC.fastq.gz \
          ${output_unPath}$(basename ${R2%%R2_001.fastq.gz})R2.unpair.QC.fastq.gz \
          ILLUMINACLIP:/home/cllin/splice_seq/NextEra.fa:0:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 HEADCROP:3 MINLEN:50
done

fastqc -t 16 ${outputPath}*.fastq.gz -o ${qcoutputPath}
