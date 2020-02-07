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

This code will also output a replicate-mean signal .wig and a corresponding .bw file for each cell line (cell1_maxmin_ReplicateMean.wig, cell1_maxmin_ReplicateMean.bw, cell2_maxmin_ReplicateMean.wig, cell2_maxmin_ReplicateMean.bw).

This code will also output a 20kb-binned residual distribution and fit of gaussian distribution.

Figures 3A,B (middle bars), 5B (middle bar), Supplementary Figures 6A,B,D,E (middle bars), 8A-D (middle bars), 10 (panel middle bars) were generated from the _mergeAdjacent.bb files.

## Calculate genome-wide P-values
Compare two cell lines with two replicates each, take mean (MEAN) and standard diviation (STD) generated from the "Statistics" step, and return p-values for each 20kb bins


```shell
python residual_4datasets_compare_Pvalue_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell2-cell1_maxmin -o cell2-cell1_maxmin_Pvalue -t aa -g utilities/hg38M.genome -w 20000 -mean MEAN -std STD
```

This code will generate a wig file and a bw file (cell2-cell1_maxmin_Pvalue.wig, cell2-cell1_maxmin_Pvalue.wig) to show P-values for each 20kb bin genome-wide.

Figures 2D (middle P-value track), 3A (bottom), 3B (middle P-value track), 5B (middle P-value track), Supplementary Figures 6A,D (bottom), 6B,E (middle P-value track), 8A,C (bottom), 8B,D (middle P-value tracks), 10 (P-value track in each panel) were generated by this code.

## Correlate changed regions with gene expression in the two cell lines

Take the bed files (cell2-cell1_maxmin_mergeAdjacent.bed and otherwaycell2-cell1_maxmin_mergeAdjacent.bed) for the changed domains between the two cell lines generated from the "Identify changed domains" step. Also take the gene expression results (gencode_expr_cellX.txt) from "TSA-Seq-2.0-Analysis/SPAD/Expression/Report" for the two cell lines to be compared.

Take H1 vs HFF comparison as example:


For genes within domains that are closer to speckles (significantly higher SON TSA-Seq scores) in HFF than in H1:

```shell
python TX_change_v2_TSA2.0.py -b HFF-H1_maxmin_mergeAdjacent.bed -g1 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_HFF_GSE100576.txt -g2 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_H1.txt -o HFF-H1 -geneID HFF-H1_geneID -y HFF/H1

# This code will generate a scatter plot showing log2-fold change of HFF/H1 against domain mean rescaled TSA-Seq score changes for all protein coding genes within the domains with significantly higher SON TSA-Seq score in HFF than in H1.

# This code will also report a number for all the genes within the domains and a number for genes with biased expression (log2-fold change of HFF/H1 > 0) comparing the two cell lines.

# This code will also generate a gene list for the genes with the biased expression (log2-fold change of HFF/H1 > 0).
```


For genes within domains that are closer to speckles (significantly higher SON TSA-Seq scores) in H1 than in HFF:

```shell
python TX_change_otherway_v2_TSA2.0.py -b otherwayHFF-H1_maxmin_mergeAdjacent.bed -g1 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_HFF_GSE100576.txt -g2 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_H1.txt -o otherwayHFF-H1 -geneID otherwayHFF-H1_geneID -y HFF/H1

# This code will generate a scatter plot showing log2-fold changes of HFF/H1 against domain mean rescaled TSA-Seq score changes for all protein coding genes within the domains with significantly higher SON TSA-Seq score in H1 than in HFF.

# This code will also report a number for all the genes within the domains and a number for genes with biased expression (log2-fold change of HFF/H1 < 0) comparing the two cell lines.

# This code will also generate a gene list for the genes with the biased expression (log2-fold change of HFF/H1 < 0)
```


Plot kernel density of the two groups of genes:

```shell
python TX_change_v2_geneHist_TSA2.0.py -b HFF-H1_maxmin_mergeAdjacent.bed -b2 otherwayHFF-H1_maxmin_mergeAdjacent.bed -g1 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_HFF_GSE100576.txt -g2 ./TSA-Seq-2.0-Analysis/SPAD/Expression/Report/gencode_expr_H1.txt -o HFF-H1_gene -y HFF/H1
```

Figure 3C and Supplementary Figures 6C,F were generated by this set of codes.


## Profile local TSA-Seq signal and histone marks for changed domains between H1 and HFFc6 cells

Software: deepTools suite (Ramirez et al. 2016) (version 3.2.1): -computeMatrix -plotHeatmap

### Profile local TSA-Seq score percentiles for the changed domains in the two cell lines
Take the .bed files (HFF-H1_maxmin_mergeAdjacent.bed and otherwayHFF-H1_maxmin_mergeAdjacent.bed) for the changed domains between the two cell lines generated from the "Identify changed domains" step. Take the .bw files for TSA-Seq enrichment score percentiles in the two cell lines generated from the "SPAD/Percentile normalizaiton" step (H1_TSA-Seq_hanning_20kbx21_percentile.bw, HFF_TSA-Seq_hanning_20kbx21_percentile.bw).

Example: for regions that are closer to speckles in HFFc6, their local (2 Mbp downstream and 2 Mbp upstream of region center) percentile profiles in H1:

```shell
computeMatrix reference-point \
 -S H1_TSA-Seq_hanning_20kbx21_percentile.bw \
 -R HFF-H1_maxmin_mergeAdjacent.bed \
 --referencePoint center \
 -a 2000000 \
 -b 2000000 \
 -out matrix_HFF-H1_mergeAdjacent_localPercInH1.tab.gz

plotHeatmap \
 -m matrix_HFF-H1_mergeAdjacent_localPercInH1.tab.gz\
 -out matrix_HFF-H1_mergeAdjacent_localPercInH1.eps \
 --heatmapHeight 15  \
 --refPointLabel region_center \
 --regionsLabel HFF-H1_regions \
 --colorList 'blue,yellow,red' \
 --plotTitle 'HFF-H1_regions_localPercentileInH1' \
 --yMin -20 \
 --yMax 40
```

This set of codes will generate of a heatmap with each line as a region and a summary plot showing mean values of all regions per position.

Figures 3D,E (top) were generated by this set of codes.

### Profile local CUT&RUN-Seq singals for the changed domains in the two cell lines

CUT&RUN-Seq data processing codes see subfolder "CUT&RUN" (by Yuchuan Wang, Ma lab, CMU), resulting files that are used: H1-hESC__H3K4me1_FE.sort.bw, H1-hESC__H3K27ac_FE.sort.bw, H1-hESC__H3K27me3_FE.sort.bw, HFFc6__H3K4me1_FE.sort.bw, HFFc6__H3K27ac_FE.sort.bw, HFFc6__H3K27me3_FE.sort.bw

Take the .bed files (HFF-H1_maxmin_mergeAdjacent.bed and otherwayHFF-H1_maxmin_mergeAdjacent.bed) for the changed domains between the two cell lines generated from the "Identify changed domains" step. Take the .bw files of CUT&RUN-Seq. 

Example: for regions that are closer to speckles in HFFc6, their local (2 Mbp downstream and 2 Mbp upstream of region center) H3K4me1 CUT&RUN-Seq profiles in H1:

```shell
computeMatrix reference-point \
 -S H1-hESC__H3K4me1_FE.sort.bw \
 -R HFF-H1_maxmin_mergeAdjacent.bed \
 --referencePoint center \
 -a 2000000 \
 -b 2000000 \
 -out matrix_HFF-H1_mergeAdjacent_localH3K4me1InH1.tab.gz

plotHeatmap \
 -m matrix_HFF-H1_mergeAdjacent_localH3K4me1InH1.tab.gz\
 -out matrix_HFF-H1_mergeAdjacent_localH3K4me1InH1.eps \
 --heatmapHeight 15  \
 --refPointLabel region_center \
 --regionsLabel HFF-H1_regions \
 --colorList 'blue,yellow,red' \
 --plotTitle 'HFF-H1_regions_localH3K4me1InH1' \
 --yMin 0 \
 --yMax 2
```

This set of codes will generate of a heatmap with each line as a region and a summary plot showing mean values of all regions per position.

Figures 3D,E (bottom) were generated by this set of codes.


## Estimate changed distances relative to speckles in TSA-Seq H1 vs K562 comparison

### Correlate scaled SON TSA-Seq scores with predicted distances in K562 cells

Take the K562_TSA-Seq_hanning_20kbx21_maxmin.wig file generated in the "Changed_Region/Rescale TSA-Seq scores" step. Take the K562_TSA-Seq_hanning_20kbx21_distance.wig file generated in the "Distance/Distance conversion" step (condition E, same TSA-Seq dataset used as for generateing K562_TSA-Seq_hanning_20kbx21_maxmin.wig).

```shell
python TSAScore_correlation_v3_temp_TSA2.0.py -w1 K562_TSA-Seq_hanning_20kbx21_maxmin.wig -w2 K562_TSA-Seq_hanning_20kbx21_distance.wig -y distance -x maxmin_score -o K562_maxminVSdistance_20kb
```

This code will generate a scatter plot (K562_maxminVSdistance_20kb.eps) to show the correlation between predicted distances to speckles and the rescaled SON TSA-Seq scores (20kb bin) in K562 cells. 

Supplementary Figure 7A was generated by this code.

### Calculate mean distance for 20kb-bins with the same rescaled SON TSA-Seq scores (1-100) and generate a dictionary

Take the same two .wig files as used for last step (K562_TSA-Seq_hanning_20kbx21_maxmin.wig, K562_TSA-Seq_hanning_20kbx21_distance.wig)

```shell
python percentileDistance_TSA2.0.py -w1 K562_TSA-Seq_hanning_20kbx21_maxmin.wig -w2 K562_TSA-Seq_hanning_20kbx21_distance.wig -o K562_maxminVSdistance.dic
```

This code will generate a file (K562_maxminVSdistance.dic) with four columns. Column 1: rescaled SON TSA-Seq scores (1-100). Column 2: number of 20kb bins for each score in column 1. Column 3: mean predicted distance of all the 20kb bins for each score in column 1. Column 4: standard deviation of predicted distance of all the 20kb bins for each score in column 1.

### Estimate distance changes relative to speckles comparing H1 vs K562, based on distance calibration in K562 cells.

Take the bed files (K562-H1_maxmin_mergeAdjacent.bed and otherwayK562-H1_maxmin_mergeAdjacent.bed) for the changed domains between the two cell lines generated from the "Identify changed domains" step. 

Take the .wig files (TSA-Seq_hanning_20kbx21_maxmin.wig) generated from the "Rescale TSA-Seq scores" step for the two cell lines each with 2 biological replicates.

Take the score-distance dictionary file generated from last step (K562_maxminVSdistance.dic). 

```shell
#For regions closer to speckles in K562 cells:
python domainRsd_convert_TSA2.0.py -b K562-H1_maxmin_mergeAdjacent.bed -p1 H1_rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -p2 H1_rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -p3 K562_rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -p4 K562_rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -d K562_maxminVSdistance.dic -o K562-H1_maxmin_mergeAdjacent_k562DistanceRsd.bed

#For regions closer to speckles in H1 cells:
python domainRsd_convert_TSA2.0.py -b otherwayK562-H1_maxmin_mergeAdjacent.bed -p1 H1_rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -p2 H1_rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -p3 K562_rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -p4 K562_rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -d K562_maxminVSdistance.dic -o otherwayK562-H1_maxmin_mergeAdjacent_k562DistanceRsd.bed
```

This code will take all the changed domains from the "mergeAdjacent.bed" file, get the domain mean rescaled SON TSA-Seq score in both cell lines (replicate mean) for each region, and convert the rescaled SON TSA-Seq scores to distances based on the provided score-distance dictionary. This code will calculate the distance residuals between the converted distences in the two cell lines for each region, always cell line 2 distance (cell2) - cell line 1 distance (cell1). This code will output a .bed file with five columns: chrom, start, end, distance(cell1), distance(cell2), distance residual (cell2-cell1)






## For all differentially expressed genes between H1 and HFF, compare relocated vs non-relocated genes

### Differential expression (DE) analysis
DE analysis was done by R with DESeq2 package (see paper Methods). Thresholds to identify differentially expressed genes: adjusted P-value < 0.01, fold change > 2, gene type: protein-coding.

Results (see sub-folder "DE"): 

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
This code will 1) generate 3 gene lists with chromosome and positions (relocated, not_relocated, all_gene); 2) report gene number for the three lists; 3) plot a box plot to compare log2-fold changes (HFF/H1) for relocated vs not_relocated genes and calculate a P-value with Welch's t-test.

Take the three gene lists (relocated, not_relocated, all_gene) and take the .bw files for rescaled TSA-Seq scores (mean values of biological replicates, generated in the "Identify changed domains" step, HFF_maxmin_ReplicateMean.bw and H1_maxmin_ReplicateMean.bw), generate scatter plots to compare their rescaled TSA-Seq scores in the two cell lines:

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

Use the relocated and not_relocated lists to do gene ontology analysis (see paper Methods).

Figures 4A,B were generated by this set of codes.

#### For DE genes with significantly higher expression in H1:

```shell
python DE_TSA_pos_v2_TSA2.0.py -b otherwayHFF-H1_maxmin_mergeAdjacent.bed -csv DE/20190421_H1vsHFF_2folds_padj0.01-DESeq2-results-with-normalized-counts-protein-coding.csv -o DE_otherwayHFF-H1 -overlap_geneID relocated -nonOverlap_geneID not_relocated -all_geneID all_gene -y H1/HFF
```
This code will 1) generate 3 gene lists with chromosome and positions (relocated, not_relocated, all_gene); 2) report gene number for the three lists; 3) plot a box plot to compare log2-fold changes (H1/HFF) for relocated vs not_relocated genes and calculate a P-value with Welch's t-test.

Take the three gene lists (relocated, not_relocated, all_gene) and take the .bw files for rescaled TSA-Seq scores (mean values of biological replicates, generated in the "Identify changed domains" step, HFF_maxmin_ReplicateMean.bw and H1_maxmin_ReplicateMean.bw), generate scatter plots to compare their rescaled TSA-Seq scores in the two cell lines:

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

Use the relocated and not_relocated lists to do gene ontology analysis (see paper Methods).

Figures 4D,E were generated by this set of codes.
