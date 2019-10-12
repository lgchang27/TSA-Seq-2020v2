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
