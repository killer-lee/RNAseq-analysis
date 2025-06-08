rm(list=ls())
control1 <- read.table("M1_count.txt",sep ="\t",col.names = c("gene_id","control"))
control2 <- read.table("M2_count.txt",sep ="\t",col.names = c("gene_id","control"))
control3 <- read.table("M3_count.txt",sep ="\t",col.names = c("gene_id","control"))
control4 <- read.table("M4_count.txt",sep ="\t",col.names = c("gene_id","control"))
control5 <- read.table("M5_count.txt",sep ="\t",col.names = c("gene_id","control"))
treat1 <- read.table("M6_count.txt",sep ="\t",col.names = c("gene_id","treat"))
treat2 <- read.table("M7_count.txt",sep ="\t",col.names = c("gene_id","traet")
treat3 <- read.table("M8_count.txt",sep ="\t",col.names = c("gene_id","treat"))
treat4 <- read.table("M9_count.txt",sep ="\t",col.names = c("gene_id","traet"))
treat5 <- read.table("M10_count.txt",sep ="\t",col.names = c("gene_id","treat"))
merged_controls <- merge(control1, control2, by = "gene_id")
merged_controls <- merge(merged_controls, control3, by = "gene_id")
merged_controls <- merge(merged_controls, control4, by = "gene_id")
merged_controls <- merge(merged_controls, control5, by = "gene_id")
merged_treat <- merge(treat1, treat2, by = "gene_id")
merged_treat <- merge(merged_treat, treat3, by = "gene_id")
merged_traet <- merge(merged_treat, treat4, by = "gene_id")
merged_treat <- merge(merged_treat, treat5, by = "gene_id")
raw_count <- merge(merged_controls,merged_treat,by="gene_id")
head(raw_count)
raw_count_filt <- raw_count[-1:-5,]
head(raw_count_filt)
ENSEMBL <- gsub("\\.\\d*\\_\\d*", "", raw_count_filt$gene_id)
rownames(raw_count_filt) <- ENSEMBL
raw_count_filt1 <- cbind(ENSEMBL, raw_count_filt)
colnames(raw_count_filt1) <- c("ensembl_gene_id", "gene_id", "control1", "control2","control3", "control4","control5","treat1", "treat2","treat3","treat4","treat5")
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
attributes <- listAttributes(mart)
head(attributes)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <- row.names(raw_count_filt1)
options(timeout = 4000000)
mouse_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'mgi_symbol', "chromosome_name", "start_position", "end_position", "band"),
  filters = 'ensembl_gene_id',
  values = my_ensembl_gene_id,
  mart = mart
)
readcount1 <- merge(raw_count_filt1,mouse_symbols,by="ensembl_gene_id")
head(readcount1)
readcount <- raw_count_filt1[ ,-1:-2]
head(readcount)
write.csv(readcount, file='readcount.csv')
library(DESeq2)
mycounts <- read.csv("readcount.csv")
head(mycounts)
rownames(mycounts) <- mycounts[,1]
head(mycounts)
condition <- factor(c(rep("control",5), rep("treat", 5)), levels = c("control", "treat"))
condition
colData <- data.frame(row.names = colnames(mycounts)[-1], condition)
mycounts <- mycounts[ , -1]
head(mycounts)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~ condition)
dds_norm <- DESeq(dds)
dds_norm
normalized_counts <- counts(dds_norm, normalized=TRUE)
head(normalized_counts)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts, file="dds_normalized_counts.xls", quote=F, sep="\t", row.names=T, col.names=T)
rld <- rlog(dds_norm, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="dds_normalized_counts_rlog.xls", quote=F, sep="\t", row.names=T, col.names=T)
res = results(dds_norm, contrast=c("condition", "treat", "control"))
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res, file="All_results.csv")
table(res$padj < 0.05)
diff_gene_deseq2 <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2, file="DEG_treat_vs_control.csv")
up <- sum(diff_gene_deseq2$log2FoldChange > 1)
down <- sum(diff_gene_deseq2$log2FoldChange < -1)
total <- nrow(diff_gene_deseq2)
cat("Total DEGs:", total, "\nUpregulated:", up, "\nDownregulated:", down)
library(biomaRt) 
library(curl)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <- row.names(diff_gene_deseq2)
mouse_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
  filters = 'ensembl_gene_id',
  values = my_ensembl_gene_id,
  mart = mart
)
sum(mouse_symbols$external_gene_name == "")
mouse_symbols_clean <- mouse_symbols[!is.na(mouse_symbols$external_gene_name) & mouse_symbols$external_gene_name != "", ]
nrow(mouse_symbols_clean)
diff_gene_deseq2 <- cbind(my_ensembl_gene_id, diff_gene_deseq2)
colnames(diff_gene_deseq2)[1] <- c("ensembl_gene_id")
head(diff_gene_deseq2)
diff_name <- merge(diff_gene_deseq2, mouse_symbols_clean, by = "ensembl_gene_id")
head(diff_name)
diff_name$regulation <- ifelse(diff_name$log2FoldChange > 0, "up", "down")
table(diff_name$regulation)
library(openxlsx)
write.xlsx(diff_name, file = "treat_vs_control_annotated.xlsx", rowNames = FALSE)
