# TSA-Seq data smoothing
This code is used to smooth the TSA-Seq enrichment scores (20kb bin) using the the TSA-Seq_20kb.wig file generated from the "TSA-Seq normalization" step.

This code will generate a TSA-Seq_hanning_20kbx21.wig file for the smoothed TSA-seq enrichment scores and a corresponding .bw file. This code will also generate a 200kb-binned TSA-Seq_hanning_20kbx21_agg_200kb.wig file averaging the smoothed scores for each 200kb bin and a corresponding .bw file.

```shell
python TSA_smooth_hanningFor20kbNonsliding_TSA2.0.py --wig TSA-Seq_20kb.wig -w 20000 -aggwin 200000 --smooth -n1 TSA-Seq_hanning_20kbx21 -n2 TSA-Seq_hanning_20kbx21_agg_200kb -g utilities/hg38F.genome


#Genome size file hg38F.genome was for female cell line (K562), hg38M.genome was for male cell lines (H1, HCT116, HFFc6).
```

Figures 2D (top), 3A (top), 3B (top), Supplementary Figures 6A (top), 6B (top), 6D (top), 6E (top), 8A (top), 8B (top), 8C (top), 8D (top), 9B, 10 (each panel top) were generated from the output TSA-Seq_hanning_20kbx21.bw file by this code.
