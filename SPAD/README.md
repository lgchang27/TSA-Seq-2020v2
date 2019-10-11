# SPAD
This set of codes is used to call Speckle Associated Domains (SPADs), compare SPADs from different cell lines, and correlate with gene expression.

## Percentile normalization
SPADs are defined as genomic regions with top 5 percentile SON TSA-Seq scores. So we first convert TSA-Seq scores into percentiles.

```shell
python TSA_percentile_norm_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -q 100 -g utilities/hg38_Gap.bed -o TSA-Seq_hanning_20kbx21_percentile -gg utilities/hg38F.genome
#Genome size file hg38F.genome was for female cell line (K562), hg38M.genome was for male cell lines (H1, HCT116, HFFc6).
```

Figure 2a (middle) was generated from the bigwig file (TSA-Seq_hanning_20kbx21_percentile.bw) by this code.

## SPADs calling
Use the percentile wig files from last step, identify bins above 95 percentile and merge adjacent bins to call SPADs

```shell
python BigPercentiles_TSA2.0.py -w TSA-Seq_hanning_20kbx21_percentile.wig -o TSA-Seq_hanning_20kbx21_percentile -p 95 -g utilities/hg38F.genome -win 20000
#Genome size file hg38F.genome was for female cell line (K562), hg38M.genome was for male cell lines (H1, HCT116, HFFc6).
```
This code will also generate simple statistics of region size and number with histograms.

Figure 2a (bottom) was generated from bigbed file (TSA-Seq_hanning_20kbx21_percentile_above_95.0_mergeAdjacent.bb) by this code.
