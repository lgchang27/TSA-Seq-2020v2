# Changed Region
This set of codes is for cell type pair-wise comparison to identify changed regions in SON TSA-Seq mapping and to correlate these regions with gene expression.

## Rescale TSA-Seq scores
Rescale TSA-Seq enrichment scores (20kb bin) linearly between their min and max values to a new 1-100 and round up to integers. 

The rescaling funciton is:

Scaled enrichment score (bin i) = (TSA-Seq enrichment score (bin i) - min) / (max - min) * 100

(min assigned to 1 instead of 0)

```shell
python plot_TSA_value_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -g utilities/hg38_Gap.bed -u 99.95 -l 0.05
```
This code will return a large and a small percentile of all ranked TSA-Seq enrichment scores (20kb bin) that will be used as the max (xx) and min (yy) scores for rescaling (e.g. 99.95% and 0.05% for HFFc6 and H1 comparison).

```shell
python TSA_max_min_norm_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -q 100 -g utilities/hg38_Gap.bed -o TSA-Seq_hanning_20kbx21_maxmin -gg utilities/hg38M.genome -max xx -min -yy
#Genome size file hg38F.genome was for female cell line (K562), hg38M.genome was for male cell lines (H1, HCT116, HFFc6).
```
This code will rescale TSA-Seq enrichmen scores into a 1 to 100 scale (showing -49 to 50 in generated wig files).

## Statistics
Statistical analysis based on biological replicates for two cell lines and generate thresholds to call changed 20kb-bins based on a P-value of 0.01 (see paper Methods)

```shell
python residual_4datasets_stat_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell1AndCell2 -P 0.01
```
This code will generate an upper (aa) and a lower (bb) shreshold to call significantly changed 20kb bins. 

The upper threshold (aa) is a positive number: residuals larger than it mean TSA-Seq signals in cell type 2 are significantly bigger than that in cell type 1.

The lower threshold (bb, bb = (-1) * aa) is a negative number: residuals smaller than it mean TSA-Seq signals in cell type 1 are significantly bigger than that in cell type 2.

This code will also output the mean (MEAN) and standard deviation (STD) for the fitted gaussian distribution (statistics). 

## Identify changed domains
Compare two cell lines with two biological replicates each, identify all 20 kb bins above the threshold generated from last step (always cell type 2 - cell type 1).

Merge adjacent bins to segment regions. Take a second threshold of domain size (100kb) and return the segments as changed domains when the region above the size threshold.

```shell
python residual_4datasets_compare_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell2-cell1_maxmin -c aa -g utilities/hg38M.genome -w 20000 -s 100000 -cell1 cell1_maxmin_ReplicateMean -cell2 cell2_maxmin_ReplicateMean
```

For regions with signals in cell type 2 significantly bigger than cell type 1: this code will output a wig file (cell2-cell1_maxmin.wig) and a bed file (cell2-cell1_maxmin.bed) for all 20 kb bins with residuals between the two cell lines, and a bed file for segmented domains (cell2-cell1_maxmin_mergeAdjacent.bed) with mean 20kb-bin residuals for each domain. And it will also output corresponding bigwig and bigbed files.


For regions with signals in cell type 1 significantly bigger than cell type 2: this code will output a wig file (otherwaycell2-cell1_maxmin.wig) and a bed file (otherwaycell2-cell1_maxmin.bed) for all 20kb bins with residuals between the two cell lines, and a bed file for segmented domains (otherwaycell2-cell1_maxmin_mergeAdjacent.bed) with mean 20kb-bin residuals for each domain. And it will also output corresponding bigwig and bigbed files.

This code will also output a replicate-mean signal bw file for each cell line.

This code will also output a bin residue distribution and fit of gaussian distribution.

Figures 2e (middle), 2f (region bars), Supplementary Figures 10a,c (middle), 10b,d (region bars), 11a (middle), 11b (region bars) were generated from the mergeAdjacent.bed (.bb) files.

## Calculate genome-wide P-values
Compare two cell lines with two replicates each, take mean (MEAN) and standard diviation (STD) generated from the "Statistics" step, and return p-values for each 20kb bins


```shell
python residual_4datasets_compare_Pvalue_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell2-cell1_maxmin -o cell2-cell1_maxmin_Pvalue -t aa -g utilities/hg38M.genome -w 20000 -mean MEAN -std STD
```

This code will generate a wig file (cell2-cell1_maxmin_Pvalue.wig) and a corresponding bigwig file to show P-values 20kb bins genome-wide.

Figures 2e (bottom), 2f (P-value track), Supplementary Figures 10a,c (bottom), 10b,d (P-value track), 11a (bottom), 11b (P-value track) were generated by this code.

## Correlate changed regions with gene expression in the two cell lines

Take the bed files (cell2-cell1_maxmin_mergeAdjacent.bed and otherwaycell2-cell1_maxmin_mergeAdjacent.bed) for the changed domains between the two cell lines generated from the "Identify changed domains" step. Also take the gene expression results (gencode_expr_cellX.txt) from "TSA-Seq-2.0-Analysis/SPAD/Expression/Report" for the two cell lines to be compared.

Take H1 vs HFF comparison as example:


For genes within domains with significantly higher SON TSA-Seq signals in HFF than in H1:

```shell
python TX_change_v2_TSA2.0.py -b HFF-H1_maxmin_mergeAdjacent.bed -g1 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_HFF_GSE100576.txt -g2 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_H1.txt -o HFF-H1 -geneID HFF-H1_geneID -y HFF/H1

# This code will generate a scatter plot showing log2-fold change of HFF/H1 against domain mean rescaled TSA-Seq scores for all protein coding genes within the domains with significantly higher TSA-Seq signal in HFF than in H1.

# This code will also report a number for all the genes within the domains and a number for genes with biased expression (log2-fold change of HFF/H1 > 0) comparing the two cell lines.

# This code will also generate a gene list for the genes with the biased expression (log2-fold change of HFF/H1 > 0).
```


For genes within domains with significantly higher SON TSA-Seq signals in H1 than in HFF:

```shell
python TX_change_otherway_v2_TSA2.0.py -b otherwayHFF-H1_maxmin_mergeAdjacent.bed -g1 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_HFF_GSE100576.txt -g2 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_H1.txt -o otherwayHFF-H1 -geneID otherwayHFF-H1_geneID -y HFF/H1

# This code will generate a scatter plot showing log2-fold changes of HFF/H1 against domain mean rescaled TSA-Seq scores for all protein coding genes within the domains with significantly higher TSA-Seq signal in H1 than in HFF.

# This code will also report a number for all the genes within the domains and a number for genes with biased expression (log2-fold change of HFF/H1 < 0) comparing the two cell lines.

# This code will also generate a gene list for the genes with the biased expression (log2-fold change of HFF/H1 < 0)
```


Plot kernel density of the two groups of genes:

```shell
python TX_change_v2_geneHist_TSA2.0.py -b HFF-H1_maxmin_mergeAdjacent.bed -b2 otherwayHFF-H1_maxmin_mergeAdjacent.bed -g1 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_HFF_GSE100576.txt -g2 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_H1.txt -o HFF-H1_gene -y HFF/H1
```

Figure 2g and Supplementary Figure 11c were generated by this set of codes
