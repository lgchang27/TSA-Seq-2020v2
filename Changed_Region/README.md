# Changed region
This set of codes was used for cell type pair-wise comparison to identify changed regions in SON TSA-Seq mapping and to correlate these regions with gene expression.

This set of codes was also used to compare K562 TSA-Seq data without and with heat shock (37C vs 42C 30min).

Examples below are shown using cell type paire-wise comparison.

## Rescale TSA-Seq scores
Rescale TSA-Seq enrichment scores (20kb bin) linearly between their min and max values to a new 1-100 scale and round up to integers. 

The rescaling funciton is:

Scaled enrichment score (bin i) = (TSA-Seq enrichment score (bin i) - min) / (max - min) * 100

(min assigned to 1 instead of 0)

```shell
python plot_TSA_value_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -g utilities/hg38_Gap.bed -u 99.95 -l 0.05
```
This code takes the .wig file for 20kb-binned smoothed TSA-Seq enrichment score (TSA-Seq_hanning_20kbx21.wig) as input. The code will return a large (-u) and a small (-l) percentile of all ranked TSA-Seq enrichment scores that will be used as the max (xx) and min (yy) scores for rescaling (e.g. 99.95% and 0.05% for HFFc6 and H1 comparison, 99.92% and 0.08% for HCT116 and H1 comparison, 99.95% and 0.05% for K562 and H1 comparison).

```shell
python TSA_max_min_norm_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -q 100 -g utilities/hg38_Gap.bed -o TSA-Seq_hanning_20kbx21_maxmin -gg utilities/hg38M.genome -max xx -min -yy

#Genome size file hg38F.genome was for female cell line (K562), hg38M.genome was for male cell lines (H1, HCT116, HFFc6).
```

This code takes the .wig file for 20kb-binned smoothed TSA-Seq enrichment score (TSA-Seq_hanning_20kbx21.wig) as input and takes a large and a small percentile of all ranked TSA-Seq enrichment scores generated from last step as max (-max) and min (-min) for the rescaling. This code will rescale TSA-Seq enrichment scores into a 1 to 100 scale (reporting as -49 to 50) and generate .wig and .bw files (TSA-Seq_hanning_20kbx21_maxmin.wig, TSA-Seq_hanning_20kbx21_maxmin.bw).

## Statistics
Statistical analysis based on biological replicates for two cell lines and generate thresholds to call changed 20kb-bins based on a P-value of 0.01 (see paper Methods)

```shell
python residual_4datasets_stat_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell1AndCell2 -P 0.01
```
This code takes the .wig files (TSA-Seq_hanning_20kbx21_maxmin.wig) generated from last step to compare two cell lines (cell1 and cell2), each with two biological replicates (rep1 and rep2). This code will generate an upper (aa) and a lower (bb) shreshold to call significantly changed 20kb bins. 

The upper threshold (aa) is a positive number: residuals larger than it mean TSA-Seq signals in cell type 2 (cell2) are significantly bigger than that in cell type 1 (cell1).

The lower threshold (bb, bb = (-1) * aa) is a negative number: residuals smaller than it mean TSA-Seq signals in cell type 1 are significantly bigger than that in cell type 2.

This code will also output the mean (MEAN) and standard deviation (STD) for the fitted gaussian distribution (statistics). 

## Identify changed domains
Compare two cell lines with two biological replicates each, identify all 20 kb bins above the threshold (-c: aa) generated from last step (always cell type 2 values - cell type 1 values).

Merge adjacent bins to segment regions. Take a second threshold (-s: 100000) of domain size (100kb) and return the segments as changed domains when the region size are above the size threshold (100kb).

```shell
python residual_4datasets_compare_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell2-cell1_maxmin -c aa -g utilities/hg38M.genome -w 20000 -s 100000 -cell1 cell1_maxmin_ReplicateMean -cell2 cell2_maxmin_ReplicateMean
```

For regions with signals in cell type 2 significantly bigger than cell type 1: this code will output a wig file (cell2-cell1_maxmin.wig) and a bed file (cell2-cell1_maxmin.bed) for all 20 kb bins with residuals between the two cell lines, and a bed file for segmented domains (cell2-cell1_maxmin_mergeAdjacent.bed) with mean 20kb-binned residuals for each domain. And it will also output corresponding .bw and .bb files for visualizaiton in genome browser.


For regions with signals in cell type 1 significantly bigger than cell type 2: this code will output a wig file (otherwaycell2-cell1_maxmin.wig) and a bed file (otherwaycell2-cell1_maxmin.bed) for all 20kb bins with residuals between the two cell lines, and a bed file for segmented domains (otherwaycell2-cell1_maxmin_mergeAdjacent.bed) with mean 20kb-binned residuals for each domain. And it will also output corresponding .bw and .bb files for visualizaiton in genome browser.

This code will also output a replicate-mean signal .bw file for each cell line.

This code will also output a 20kb-binned residual distribution and fit of gaussian distribution.

Figures 3A,B (middle bars), Supplementary Figures 6A,B,D,E (middle bars), 8A-D (middle bars) were generated from the _mergeAdjacent.bb files.

## Calculate genome-wide P-values
Compare two cell lines with two replicates each, take mean (MEAN) and standard diviation (STD) generated from the "Statistics" step, and return p-values for each 20kb bins


```shell
python residual_4datasets_compare_Pvalue_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell2-cell1_maxmin -o cell2-cell1_maxmin_Pvalue -t aa -g utilities/hg38M.genome -w 20000 -mean MEAN -std STD
```

This code will generate a wig file and a bw file (cell2-cell1_maxmin_Pvalue.wig, cell2-cell1_maxmin_Pvalue.wig) to show P-values for each 20kb bin genome-wide.














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


## For all differentially expressed genes between H1 and HFF, compare relocated vs non-relocated genes

### Differential expression (DE) analysis
DE analysis was done by R with DESeq2 package (see paper Methods). Thresholds to identify differentially expressed genes: adjusted P-value < 0.01, fold change > 2, gene type: protein-coding.

Results: 

Significantly higher expression in HFF: 
DE/20190421_HFFvsH1_2folds_padj0.01-DESeq2-results-with-normalized-counts-protein-coding.csv

Significantly higher expression in H1: 
DE/20190421_H1vsHFF_2folds_padj0.01-DESeq2-results-with-normalized-counts-protein-coding.csv

### Seperate relocated and non-relocated DE genes, and compare their log2-fold changes

Take the bed files (HFF-H1_maxmin_mergeAdjacent.bed and otherwayHFF-H1_maxmin_mergeAdjacent.bed) for the changed domains between the two cell lines generated from the "Identify changed domains" step. Also take the .csv files generated from DE analysis. Find genes that are entirely located within the changed domains as relocated genes. Find genes that are not entirely located within the changed domains as non-relocated genes.

#### For DE genes with significantly higher expression in HFF:

```shell
python DE_TSA_pos_v2_TSA2.0.py -b HFF-H1_maxmin_mergeAdjacent.bed -csv DE/20190421_HFFvsH1_2folds_padj0.01-DESeq2-results-with-normalized-counts-protein-coding.csv -o DE_HFF-H1 -overlap_geneID relocated -nonOverlap_geneID not_relocated -all_geneID all_gene -y HFF/H1
```
This code will 1) generate 3 gene lists with chromosome and positions (relocated, not_relocated, all_gene); 2) report gene number for the three lists; 3) plot a box plot to compare log2-fold changes (HFF/H1) for relocated vs not_relocated genes and calculate a P-value with two-tailed t-test.

Take the three gene lists (relocated, not_relocated, all_gene) and take the bigwig files for rescaled TSA-Seq scores (mean values of biological replicates, generated in the "Identify changed domains" step), generate scatter plots to compare their rescaled TSA-Seq scores in the two cell lines:

all_gene:
```shell
python gene_percentile_TSA2.0.py -g all_gene -p1 HFF_maxmin_ReplicateMean.bw -p2 H1_maxmin_ReplicateMean.bw -o all_DE_gene_maxmin -x HFFc6 -y H1
```

relocated:
```shell
python gene_percentile_TSA2.0.py -g relocated -p1 HFF_maxmin_ReplicateMean.bw -p2 H1_maxmin_ReplicateMean.bw -o relocated_DE_gene_maxmin -x HFFc6 -y H1
```

not_relocated:
```shell
python gene_percentile_TSA2.0.py -g not_relocated -p1 HFF_maxmin_ReplicateMean.bw -p2 H1_maxmin_ReplicateMean.bw -o not_relocated_DE_gene_maxmin -x HFFc6 -y H1
```

In the scatter plots, dashed lines show the threshold to call significantly changed domains identified in the "Statistics" step.

Also use the relocated and not_relocated lists to do gene ontology analysis (see paper Methods).

Supplementary Figures 12a,b were generated by this set of codes.

#### For DE genes with significantly higher expression in H1:

```shell
python DE_TSA_pos_v2_TSA2.0.py -b otherwayHFF-H1_maxmin_mergeAdjacent.bed -csv DE/20190421_H1vsHFF_2folds_padj0.01-DESeq2-results-with-normalized-counts-protein-coding.csv -o DE_otherwayHFF-H1 -overlap_geneID relocated -nonOverlap_geneID not_relocated -all_geneID all_gene -y H1/HFF
```
This code will 1) generate 3 gene lists with chromosome and positions (relocated, not_relocated, all_gene); 2) report gene number for the three lists; 3) plot a box plot to compare log2-fold changes (H1/HFF) for relocated vs not_relocated genes and calculate a P-value with two-tailed t-test.

Take the three gene lists (relocated, not_relocated, all_gene) and take the bigwig files for rescaled TSA-Seq scores (mean values of biological replicates, generated in the "Identify changed domains" step), generate scatter plots to compare their rescaled TSA-Seq scores in the two cell lines:

all_gene:
```shell
python gene_percentile_TSA2.0.py -g all_gene -p1 HFF_maxmin_ReplicateMean.bw -p2 H1_maxmin_ReplicateMean.bw -o all_DE_gene_maxmin -x HFFc6 -y H1
```

relocated:
```shell
python gene_percentile_TSA2.0.py -g relocated -p1 HFF_maxmin_ReplicateMean.bw -p2 H1_maxmin_ReplicateMean.bw -o relocated_DE_gene_maxmin -x HFFc6 -y H1
```

not_relocated:
```shell
python gene_percentile_TSA2.0.py -g not_relocated -p1 HFF_maxmin_ReplicateMean.bw -p2 H1_maxmin_ReplicateMean.bw -o not_relocated_DE_gene_maxmin -x HFFc6 -y H1
```

In the scatter plots, dashed lines show the threshold to call significantly changed domains identified in the "Statistics" step.

Also use the relocated and not_relocated lists to do gene ontology analysis (see paper Methods).

Supplementary Figures 12d,e were generated by this set of codes.
