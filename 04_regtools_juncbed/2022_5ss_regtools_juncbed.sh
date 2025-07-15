#!/usr/bin/sh
projectPath=/home/N420/lab_data/cllin/run3_run4/ang_chu_test/5ss/
bamPath=${projectPath}output/
countPath=${projectPath}junction/
juncPath=${bamPath}spliced/
unalignPath=${projectPath}unAlign/

#junction extract by regtools
bamSList=(${juncPath}*_S.bam)

for bamS in ${bamSList[@]}
do
regtools junctions extract -a 0 -m 20 -s 0 $bamS -o ${countPath}$(basename ${bamS%%.bam}).juncbed
  echo "$(basename ${bamS%%.bam}.juncbed) DONE"
done

#annotate junction by regtools
juncbedList=(${countPath}*_S.juncbed)

for juncbed in ${juncbedList[@]}
do
regtools junctions annotate ${juncbed} /home/N420/lab_data/cllin/run3_run4/reordered_5ss-1_form_mod.fasta /home/N420/lab_data/cllin/run3_run4/5ss_exon.gtf -o ${countPath}$(basename ${juncbed})_exon.annotation
echo "$(basename $juncbed) exon DONE"
done

for juncbed in ${juncbedList[@]}
do
regtools junctions annotate ${juncbed} /home/N420/lab_data/cllin/run3_run4/reordered_5ss-1_form_mod.fasta /home/N420/lab_data/cllin/run3_run4/5ss_intron.gtf -o ${countPath}$(basename ${juncbed})_intron.annotation
echo "$(basename $juncbed) intron DONE"
done

#unalign sequence
unalignList=(${unalignPath}*.fastq)
for unalign in ${unalignList[@]}
do
fastx_collapser -i $unalign | fasta_formatter -t -o ${unalignPath}$(basename ${unalign%%.fastq}).unalign.txt
echo "$(basename $unalign) DONE"
done
