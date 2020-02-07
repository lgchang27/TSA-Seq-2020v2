#!/bin/bash

dir="/YOUR_DIR_fastq"
input="INPUT_FILE"
trim_dir="/YOUR_DIR/trimmomatic/Trimmomatic-0.39"
adapter="/YOUR_DIR/Trimmomatic-0.39/adapters"

mkdir -p trim_bam
rm -rf trim_bam/*

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

	java -jar $trim_dir/trimmomatic-0.39.jar PE -threads 1 ${dir}/${pair1}.fastq.gz ${dir}/${pair2}.fastq.gz trim_bam/${pair1}.paired.fastq.gz trim_bam/${pair1}.unpaired.fastq.gz trim_bam/${pair2}.paired.fastq.gz trim_bam/${pair2}.unpaired.fastq.gz ILLUMINACLIP:$adapter/TruSeq3-PE-2.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

done


