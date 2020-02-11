# TSA-Seq analysis

Python 2.7

Required packages: bx-python, numpy, scipy, progressbar, pytabix, pyfasta

Tested in Linux

## TSA-Seq normalization

TSA-Seq normalization pipeline and softwares are at https://github.com/zocean/Norma. This pipeline will take raw sequencing fastq files (pulldown and input), align reads to genome, remove PCR duplicates, generate TSA-Seq enrichment scores (20kb bin) as .wig and .bw files (e.g. TSA-Seq_20kb.wig, TSA-Seq_20kb.bw)

Figures 1E (top), 2A (top), Supplementary Figure 3C were generated from the output .bw files from this set of codes.

## TSA-Seq data smoothing

This code is used to smooth 20-kb binned TSA-Seq enrichment scores. (folder: TSA-Seq-2020/Smooth)

## Distance and residuals

This set of codes is used for distance prediction from smoothed TSA-Seq enrichment scores and for calculation of distance residuals between different TSA-Seq conditions. (folder: TSA-Seq-2020/Distance)

## SPAD

This set of codes is used to call Speckle Associated Domains (SPADs), compare SPADs in different cell lines, and correlate SPADs with gene expression. (folder: TSA-Seq-2020/SPAD)

## Changed regions

This set of codes is used for cell type pair-wise comparison to identify changed regions in SON TSA-Seq mapping and to correlate these regions with gene expression. (folder: TSA-Seq-2020/Changed_Region)
