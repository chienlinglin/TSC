#!/usr/bin/sh
projectPath=/home/N420/lab_data/cllin/run3_run4/
readPath=${projectPath}/5ss_qc/
indexPath=${projectPath}/
outputPath=${projectPath}/ang_chu_test/5ss/output/
unAlignPath=${projectPath}/ang_chu_test/5ss/unAlign/
recordPath=${projectPath}/ang_chu_test/5ss/record/
list=(${readPath}*R1.QC.fastq.gz)

mkdir ${projectPath}/ang_chu_test/5ss ${outputPath} ${unAlignPath} ${recordPath}

for U1 in ${list[@]}
do 
  U2=${readPath}$(basename ${U1%%R1.QC.fastq.gz})R2.QC.fastq.gz
  hisat2 \
      -p 30 \
      -x ${indexPath}/reordered_5ssindex/reorder5_index \
      -1 $U1 \
      -2 $U2 \
      --mp 10,4 \
      --rdg 6,4 \
      --rfg 6,4 \
      -k 1 \
      --non-deterministic \
      --pen-noncansplice 100 \
      -S ${outputPath}$(basename ${U1%%_L001_R1.QC.fastq.gz}).sam \
      --un-conc ${unAlignPath}$(basename ${U1%%_L001_R1.QC.fastq.gz}).fastq \
      >& ${recordPath}$(basename ${U1%%_L001_R1.QC.fastq.gz}).record
done


