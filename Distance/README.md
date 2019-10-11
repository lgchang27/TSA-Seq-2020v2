# Distance and residuals
This set of codes is used for distance prediction by TSA-Seq data (smoothed) and calculation of distance residuals between different TSA-Seq conditions.

## Distance calculation
Cytological distances (um) (x) to the targetted nuclear body are predicted from TSA-Seq signals (y) by a fomula y = y0 + A * e^(R0 * x). (See paper Method)

### Find y0, A and R0

```shell
python plot_TSA_value_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -g utilities/hg38_Gap.bed
```
This code will find the max and min TSA-Seq enrichment scores.

y0 = ymin = 2^(min)

A = ymax - ymin = 2^(max) - 2^(min)

R0 is based on FISH calibration and fitting (see paper Method and Sup Table 2)

### Distance conversion

```shell
python TSAtoDistance_v2_TSA2.0.py -i TSA-Seq_hanning_20kbx21.wig -o TSA-Seq_hanning_20kbx21_distance -y0 xx -A yy -R0 zz -g utilities/hg38F.genome -gap utilities/hg38_Gap.bed
```

This code will generate a ditance track (TSA-Seq_hanning_20kbx21_distance.bw) that was used to generate figure 1e (middle) and supplementary figure 4a (top).

## Residual calculation

```shell
python distance_residual_v2_TSA2.0.py -w1 TSA_Seq_hanning_20kbx21_distance_conditionX.wig -w2 TSA_Seq_hanning_20kbx21_distance_conditionY.wig -o distance_residual -gap utilities/hg38_Gap.bed -g utilities/hg38F.genome -w 20000
```
This code will compare two distance .wig files, retun distance residual wig and bigwig files, and a histgram of absolute residuals between the two files.

The generated residual tracks were used for figure 1e (bottom) and supplementary figure 4a (bottom).

The generated residual histograms were used for figure 1f and supplementary figure 4b.
