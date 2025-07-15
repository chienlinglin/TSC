#!/usr/bin/sh
projectPath=/home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/
bamPath=${projectPath}output/
countPath=${projectPath}junction/
juncPath=${bamPath}spliced/
unalignPath=${projectPath}unAlign/

mkdir ${juncPath} ${countPath}
 

#sam to bam(q60 ,sort)
unbamList=(${bamPath}*.sam)
for unbam in ${unbamList[@]}
do
samtools view -q 60 -b $unbam | samtools sort > ${bamPath}$(basename ${unbam%%.sam})_q60_sorted.bam
samtools index ${bamPath}$(basename ${unbam%%.sam})_q60_sorted.bam
echo "$(basename $unbam) DONE"
done

#combine 2 batches of reads
#samtools merge /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-1_t_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-1_S4_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-1_S12_q60_sorted.bam
#echo "5-1 merge DONE"
#samtools merge /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-2_t_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-2_S5_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-2_S13_q60_sorted.bam
#echo "5-2 merge DONE"
#samtools merge /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-3_t_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-3_S6_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/output/5-3_S14_q60_sorted.bam
#echo "5-3 merge DONE"

#keep the pair end reads if aligning to same reference
bamList=(${bamPath}*_q60_sorted.bam)
for bam in ${bamList[@]}
do
samtools view -h $bam | awk '$7 ~ /=/ || $1 ~ /^@/' | samtools view -bS > ${bamPath}$(basename ${bam%%.bam})_pairend.bam
samtools index ${bamPath}$(basename ${bam%%.bam})_pairend.bam
echo "$(basename $bam) pairend filter DONE"
done

#filter mt info by bedtools intersect
sortbamList=(${bamPath}*pairend.bam)
for sobam in ${sortbamList[@]}
do
samtools index $sobam
bedtools intersect -abam $sobam -b /home/N420/lab_data/cllin/run3_run4/hi_5ss_mt_info.bed > ${bamPath}$(basename ${sobam%%.bam})_mt_info.bam
bedtools intersect -v -abam $sobam -b /home/N420/lab_data/cllin/run3_run4/hi_5ss_mt_info.bed > ${bamPath}$(basename ${sobam%%.bam})_no_mt_info.bam
echo "filter mt info $(basename $sobam) DONE"
done

#with or no mt_info bam file separate spliced, unspliced
bamList=(${bamPath}*_mt_info.bam)
for bam in ${bamList[@]} 
do
samtools view -h $bam | awk '$6 !~ /N/ || $1 ~ /^@/' | samtools view -bS > ${juncPath}$(basename ${bam%%.bam})_N.bam
samtools index ${juncPath}$(basename ${bam%%.bam})_N.bam
done

for bam in ${bamList[@]} 
do
samtools view -h $bam | awk '$6 ~ /N/ || $1 ~ /^@/' | samtools view -bS > ${juncPath}$(basename ${bam%%.bam})_S.bam
samtools index ${juncPath}$(basename ${bam%%.bam})_S.bam
echo "$(basename $bam) DONE"
done

#Index
tbamList=(${bamPath}*.bam)
for tbam in ${bamList[@]}
do
  samtools \
     index \
      $tbam
  echo "$(basename $tbam) DONE"
done

#Count total by bedtool
bedtools multicov -q 60 -bams ${bamPath}*S*_q60_sorted.bam -bed /home/N420/lab_data/cllin/run3_run4/5ss-1_new_lib.bed > ${countPath}2022_5ss_q60_total.count
echo "6 file bedtools multicov DONE"
bedtools multicov -q 60 -bams ${bamPath}*pairend.bam -bed /home/N420/lab_data/cllin/run3_run4/5ss-1_new_lib.bed > ${countPath}2022_5ss_q60_pairend.count
echo "merged 3 file total bedtools multicov DONE"
bedtools multicov -q 60 -bams ${bamPath}*_mt_info.bam -bed /home/N420/lab_data/cllin/run3_run4/5ss-1_new_lib.bed > ${countPath}2022_5ss_q60_pairend_mt_info.count
echo "mt info bedtools multicov DONE"

#Count each sample's spliced and unspliced by bedtools, separated by awk
list=(${juncPath}*_N.bam)
for U1 in ${list[@]}
do 
U2=${juncPath}$(basename ${U1%%_N.bam})_S.bam
bedtools multicov -bams $U1 $U2 -bed /home/N420/lab_data/cllin/run3_run4/5ss-1_new_lib.bed > ${countPath}$(basename ${U1%%_N.bam})_2022_5ss_spliced_unspliced.count
  echo "$(basename $U1) DONE"
done

