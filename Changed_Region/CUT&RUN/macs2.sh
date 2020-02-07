#!/bin/bash

data_list="data.list.txt"
igg_list="IgG.list.txt"

dir="bam_dir"
hg38_size="hg38.chrom.sizes"


for i in `cut -f 1 $data_list`
do
	cell=`grep $i $data_list|cut -f 2`
	target=`grep $i $data_list|cut -f 3`
	input=`awk -v cell=$cell '{if($2==cell){print $1}}' $igg_list`
	echo $i $cell $target $input

	mkdir -p $cell"__"$target

	macs2 callpeak -t ${dir}/${i}/${i}.bam -c ${dir}/${input}/${input}.bam -B --nomodel --extsize 100 --SPMR -g hs --outdir $cell"__"$target -n $cell"__"$target
	macs2 bdgcmp -t $cell"__"$target/$cell"__"${target}_treat_pileup.bdg -c $cell"__"$target/$cell"__"${target}_control_lambda.bdg -o $cell"__"$target/$cell"__"${target}_FE.bdg -m FE
	macs2 bdgcmp -t $cell"__"$target/$cell"__"${target}_treat_pileup.bdg -c $cell"__"$target/$cell"__"${target}_control_lambda.bdg -o $cell"__"$target/$cell"__"${target}_logLR.bdg -m logLR -p 0.00001

	sortBed -i $cell"__"$target/$cell"__"${target}_FE.bdg > $cell"__"$target/$cell"__"${target}_FE.sort.bdg
	bedGraphToBigWig $cell"__"$target/$cell"__"${target}_FE.sort.bdg $hg38_size $cell"__"$target/$cell"__"${target}_FE.sort.bw

	bedToBigBed -as=bigNarrowPeak.as -type=bed6+4 $cell"__"$target/$cell"__"${target}_peaks.narrowPeak $hg38_size $cell"__"$target/$cell"__"${target}_peaks.narrowPeak.bb


done




