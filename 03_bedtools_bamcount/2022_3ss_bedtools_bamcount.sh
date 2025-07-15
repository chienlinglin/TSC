#!/usr/bin/sh
projectPath=/home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss
unalignPath=${projectPath}/unAlign/
bamPath=${projectPath}/output/
#add new dir
juncPath=${bamPath}/spliced/
countPath=${projectPath}/junction/

mkdir ${juncPath} ${countPath}


#sam to bam(q60 ,sort)
unbamList=(${bamPath}*.sam)
for unbam in ${unbamList[@]}
do
samtools view -q 60 -b $unbam | samtools sort > ${bamPath}$(basename ${unbam%%.sam})_q60_sorted.bam
samtools index ${bamPath}$(basename ${unbam%%.sam})_q60_sorted.bam
echo "$(basename $unbam) DONE"
done
##########################################################################################
#combine 2 batches of reads
#samtools merge /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-1_t_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-1_S1_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-1_S9_q60_sorted.bam
#echo "3-1 merge DONE"
#samtools merge /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-2_t_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-2_S2_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-2_S10_q60_sorted.bam
#echo "3-2 merge DONE"
#samtools merge /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-3_t_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-3_S3_q60_sorted.bam /home/N420/lab_data/cllin/run3_run4/ang_chu_test/3ss/output/3-3_S11_q60_sorted.bam
#echo "3-3 merge DONE"
##########################################################################################
#keep the pair end reads if aligning to same reference
bamList=(${bamPath}*_q60_sorted.bam)
for bam in ${bamList[@]}
do
samtools view -h $bam | awk '$7 ~ /=/ || $1 ~ /^@/' | samtools view -bS > ${bamPath}$(basename ${bam%%.bam})_pairend.bam
samtools index ${bamPath}$(basename ${bam%%.bam})_pairend.bam
echo "$(basename $bam) pairend filter DONE"
done
##########################################################################################
#filter mt info by bedtools intersect
sortbamList=(${bamPath}*pairend.bam)
for sobam in ${sortbamList[@]}
do
samtools index $sobam
bedtools intersect -abam $sobam -b /home/N420/lab_data/cllin/run3_run4/hi_3ss_mt_info.bed > ${bamPath}$(basename ${sobam%%.bam})_mt_info.bam
bedtools intersect -v -abam $sobam -b /home/N420/lab_data/cllin/run3_run4/hi_3ss_mt_info.bed > ${bamPath}$(basename ${sobam%%.bam})_no_mt_info.bam
echo "filter mt info $(basename $sobam) DONE"
done
##########################################################################################
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
##########################################################################################
#Index
tbamList=(${bamPath}.bam)
for tbam in ${bamList[@]}
do
  samtools \
     index \
      $tbam
  echo "$(basename $tbam) DONE"
done
##########################################################################################
#Count total by bedtool
bedtools multicov -q 60 -bams ${bamPath}*S*_q60_sorted.bam -bed /home/N420/lab_data/cllin/run3_run4/3ss_new_lib.bed > ${countPath}2022_3ss_q60_total.count
echo "6 file bedtools multicov DONE"
bedtools multicov -q 60 -bams ${bamPath}*pairend.bam -bed /home/N420/lab_data/cllin/run3_run4/3ss_new_lib.bed > ${countPath}2022_3ss_q60_pairend.count
echo "merged 3 file total bedtools multicov DONE"
bedtools multicov -q 60 -bams ${bamPath}*_mt_info.bam -bed /home/N420/lab_data/cllin/run3_run4/3ss_new_lib.bed > ${countPath}2022_3ss_q60_mt_info.count
echo "mt info bedtools multicov DONE"

#Count each sample's spliced and unspliced by bedtools, separated by awk
list=(${juncPath}*_N.bam)
for U1 in ${list[@]}
do 
U2=${juncPath}$(basename ${U1%%_N.bam})_S.bam
bedtools multicov -bams $U1 $U2 -bed /home/N420/lab_data/cllin/run3_run4/3ss_new_lib.bed > ${countPath}$(basename ${U1%%_N.bam})_2022_3ss_spliced_unspliced.count
  echo "$(basename $U1) DONE"
done


