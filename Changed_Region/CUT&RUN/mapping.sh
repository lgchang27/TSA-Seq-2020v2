#!/bin/bash

dir="trim_bam"
meta="metadata.tsv"

index="bowtie2_index"
threads=20

cut -f 2,3,4,17,18 $meta |awk '{
	if($5=="1"){
		print $0
	}
}' > sample.txt

mkdir -p bam
rm -rf bam/*

input="sample.txt"

for i in `sed 's/\t/_/g' $input`
do
	echo $i
	j=${i%_*}
	pair2=${j##*_}
	j=${j%_*}
	pair1=${j##*_}
	j=${j%_*}
	Experiment_Accession=${j##*_}
	j=${j%_*}
	Experiment_Set_Accession=${j}

	echo $Experiment_Set_Accession $Experiment_Accession $pair1 $pair2 

	mkdir -p bam/${pair1}
	
	bowtie2 -p 20 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x $index -1 ${dir}/${pair1}.paired.fastq.gz -2 ${dir}/${pair2}.paired.fastq.gz > bam/${pair1}/${pair1}.bam
	samtools sort bam/${pair1}/${pair1}.bam -o bam/${pair1}/${pair1}.sort.bam
	samtools rmdup bam/${pair1}/${pair1}.sort.bam bam/${pair1}/${pair1}.dedup.bam

done

#echo "Done" > ${input}.log




