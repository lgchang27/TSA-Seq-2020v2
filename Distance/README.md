# Distance and residuals
This set of codes is used for distance prediction from smoothed TSA-Seq enrichment scores (TSA-Seq_hanning_20kbx21.wig) and calculation for distance residuals between different TSA-Seq conditions.

## Distance calculation
Cytological distances (um) (x) to the targetted nuclear body are predicted from TSA-Seq signals (y) by a fomula y = y0 + A * e^(R0 * x). (See paper Methods)

### Find y0, A and R0

```shell
python plot_TSA_value_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -g utilities/hg38_Gap.bed -u 100 -l 0
```
This code will find the max and min smoothed TSA-Seq enrichment scores based on which y0 and A will be obtained according to:

y0 = ymin = 2^(min)

A = ymax - ymin = 2^(max) - 2^(min)

R0 is obtained based on FISH calibration and fitting (see paper Methods and Sup Table 2)

### Distance conversion

```shell
python TSAtoDistance_v2_TSA2.0.py -i TSA-Seq_hanning_20kbx21.wig -o TSA-Seq_hanning_20kbx21_distance -y0 xx -A yy -R0 zz -g utilities/hg38F.genome -gap utilities/hg38_Gap.bed
```

This code will take the parameters y0, A, R0 to convert smoothed TSA-Seq enrichment scores (TSA-Seq_hanning_20kbx21.wig) to distances in um for each 20kb bin (TSA-Seq_hanning_20kbx21_distance.wig, TSA-Seq_hanning_20kbx21_distance.bw).

Figure 1E (middle) and Supplementary Figure 2A (top) were generated from the .bw files (TSA-Seq_hanning_20kbx21_distance.bw) by thie code.

## Residual calculation

```shell
python distance_residual_v2_TSA2.0.py -w1 TSA_Seq_hanning_20kbx21_distance_conditionX.wig -w2 TSA_Seq_hanning_20kbx21_distance_conditionY.wig -o distance_residual -gap utilities/hg38_Gap.bed -g utilities/hg38F.genome -w 20000
```
This code will compare two 20kb-binned distance .wig files (TSA_Seq_hanning_20kbx21_distance.wig) from different TSA-Seq conditions, return distance residuals for each 20kb bin (distance_residual.wig, distance_residual.bw). This code will also generate a histgram showing absolute residuals between the two conditions (distance_residual_hst.eps).

Figure 1E (bottom) and Supplementary Figure 2A (bottom) were generated from the .bw files (distance_residual.bw).

Figure 1F and Supplementary Figure 2B were generated from the residual histograms (distance_residual_hst.eps).

