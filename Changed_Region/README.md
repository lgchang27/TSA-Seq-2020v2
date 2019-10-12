# Changed Region
This set of codes is for cell type pair-wise comparison to identify changed regions in SON TSA-Seq mapping.

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
Statistical analysis based on biological replicates for two cell lines and generate thresholds to call changed 20kb-bins based on a P-value of 0.01

```shell
python residual_4datasets_stat_TSA2.0.py -c1r1 cell1Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c1r2 cell1Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r1 cell2Rep1_TSA-Seq_hanning_20kbx21_maxmin.wig -c2r2 cell2Rep2_TSA-Seq_hanning_20kbx21_maxmin.wig -o cell1AndCell2 -P 0.01
```
This code will generate an upper (aa) and a lower (bb) shreshold to call significantly changed 20kb bins. 

The upper threshold (aa) is a positive number, residuals larger than which mean TSA-Seq signals in cell type 2 are significantly bigger than that in cell type 1.

The lower threshold (bb) is a negative number, residuals smaller than which mean TSA-Seq signals in cell type 1 are significantly bigger than that in cell type 2.
