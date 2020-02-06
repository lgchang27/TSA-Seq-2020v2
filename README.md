# TSA-Seq-2.0

## TSA-Seq normalization

TSA-Seq normalization pipeline and softwares are at https://github.com/zocean/Norma. This pipeline will take raw sequencing fastq files (pulldown and input), align reads to genome, remove PCR duplicates, generate TSA-Seq enrichment scores (20kb bin) as .wig and .bigwig files.

Figures 1E (top), 2A (top), Supplementary Figures 3C were generated from the output bigwig file.

## TSA-Seq data smoothing

This code is used to smooth 20-kb binned TSA-Seq enrichment scores. (TSA-Seq-2020/Smooth)

## Distance and residuals

This set of codes is used for distance prediction from smoothed TSA-Seq enrichment scores and calculation of distance residuals between different TSA-Seq conditions. (TSA-Seq-2020/Distance)

## SPAD

This set of codes is used to call Speckle Associated Domains (SPADs), compare SPADs in different cell lines, and correlate SPADs with gene expression. (TSA-Seq-2020/SPAD)

## Changed regions

This set of codes is used for cell type pair-wise comparison to identify changed regions in SON TSA-Seq mapping and to correlate these regions with gene expression. (TSA-Seq-2020/Changed_Region)
