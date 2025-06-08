# RNAseq-analysis
A regular script for upstream analysis of rnaseq and my differential analysis plan

# Desciption
RNA_seq.py：A script for upstream analysis of rnaseq. If you want to use it, you need to replace the software installation path I specified. At the same time, it has to be mentioned that for my script, my transcriptome assembly step did not actually work. You can delete this part.

RNAseq.sh：A script for executing RNAseq.py files on Linux. If you want to use RNAseq.py, you can input parameters similar to RNAseq.sh to execute RNAseq.py

RNAseq_differential_analysis.r:A R language script that uses deseq2 for differential analysis. It is worth noting that if you use this script, you will need to modify the merging method of my count files. Also, please note that this script does not include drawing
