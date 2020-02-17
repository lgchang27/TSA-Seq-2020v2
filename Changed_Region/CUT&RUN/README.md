
CUT&RUN Data Process

CUT&RUN data in H1 and HFF is downloaded from 4DN data portal. First, we trimmed the paired end reads with trimmomatic (1). The used parameters are

ILLUMINACLIP:$adapter/TruSeq3-PE-2.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25.

Then we used bowtie 2 (2) to align the trimmed reads to reference genome hg38. We set the minimum fragment length to be 10 and maximum fragment length to be 700. Then we sort and remove duplicate reads with samtools (3). Finally, we used MACS2 (4) to generate fold enrichment (FE) tracks and perform peak calling. 

Reference:

1. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.

2. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357.

3. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The sequence alignment/map format and SAMtools. Bioinformatics, 25(16), 2078-2079.

4. Zhang, Y., Liu, T., Meyer, C. A., Eeckhoute, J., Johnson, D. S., Bernstein, B. E., ... & Liu, X. S. (2008). Model-based analysis of ChIP-Seq (MACS). Genome biology, 9(9), R137.
