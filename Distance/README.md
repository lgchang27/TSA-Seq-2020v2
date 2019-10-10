# Distance and residuals
This set of codes is used for distance prediction by TSA-Seq data (smoothed) and calculation of distance residuals between different TSA-Seq conditions.

## Distance calculation
Cytological distances (um) (x) to a targetted nuclear body are predicted from TSA-Seq signals (y) by a fomula y = y0 + A*e^(R0*x). (See paper Method)

### Find R0 and A

```shell
python plot_TSA_value_TSA2.0.py -w TSA-Seq_hanning_20kbx21.wig -g utilities/hg19_Gap.bed
