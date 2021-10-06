# Installing DESeq2 package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggrepel)

# Sample info
sample_info <- read_csv("/projectnb/bf528/users/tinman/Project3/programmer/deseq/group_4_rna_info.csv")
sample_info_control <- sample_info %>% filter(mode_of_action == "Control")

# extract the treatment samples
treatment_sample = sample_info %>% filter(!mode_of_action %in% "Control")

# extract the control samples
control_sample = sample_info %>% filter(mode_of_action %in% "Control")


# Reading individual count tables
SRR1177960 <- fread("featureCounts/count_files/SRR1177960_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177967 <- fread("featureCounts/count_files/SRR1177967_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177968 <- fread("featureCounts/count_files/SRR1177968_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177971 <- fread("featureCounts/count_files/SRR1177971_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177984 <- fread("featureCounts/count_files/SRR1177984_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177985 <- fread("featureCounts/count_files/SRR1177985_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177986 <- fread("featureCounts/count_files/SRR1177986_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178049 <- fread("featureCounts/count_files/SRR1178049_star_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)

# Using only the geneid and count columns
SRR1177960 <- SRR1177960[,c(1,7)] %>% dplyr::rename(SRR1177960 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1177960_star_Aligned.sortedByCoord.out.bam")
SRR1177967 <- SRR1177967[,c(1,7)] %>% dplyr::rename(SRR1177967 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1177967_star_Aligned.sortedByCoord.out.bam")
SRR1177968 <- SRR1177968[,c(1,7)] %>% dplyr::rename(SRR1177968 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1177968_star_Aligned.sortedByCoord.out.bam")
SRR1177971 <- SRR1177971[,c(1,7)] %>% dplyr::rename(SRR1177971 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1177971_star_Aligned.sortedByCoord.out.bam")
SRR1177984 <- SRR1177984[,c(1,7)] %>% dplyr::rename(SRR1177984 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1177984_star_Aligned.sortedByCoord.out.bam")
SRR1177985 <- SRR1177985[,c(1,7)] %>% dplyr::rename(SRR1177985 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1177985_star_Aligned.sortedByCoord.out.bam")
SRR1177986 <- SRR1177986[,c(1,7)] %>% dplyr::rename(SRR1177986 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1177986_star_Aligned.sortedByCoord.out.bam")
SRR1178049 <- SRR1178049[,c(1,7)] %>% dplyr::rename(SRR1178049 = "/projectnb/bf528/users/tinman/Project3/datacurator/star_results2/SRR1178049_star_Aligned.sortedByCoord.out.bam")

# Merging by geneid
merged_gid <- Reduce(function(...) merge(..., by = "Geneid", all=TRUE), list(SRR1177960, SRR1177967, SRR1177968, SRR1177971, SRR1177984, SRR1177985, SRR1177986, SRR1178049))
rm(SRR1177960, SRR1177967, SRR1177968, SRR1177971, SRR1177984, SRR1177985, SRR1177986, SRR1178049)

# Writing CSV file
write.csv(merged_gid, file="featureCounts_concatenated.csv", row.names=FALSE)

# Box plot of distribution
merged_box <- merged_gid %>% pivot_longer(!Geneid, names_to = "Sample", values_to = "Count")
png(filename = "boxplots/box_count_distribution.png")
merged_box %>% ggplot(mapping = aes(x = Sample, y = Count)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Boxplot of Count Distributions") +
  ylim(0, 1000)
dev.off()

### Part 3: RNA-Seq Differential Expression using DESeq2

# Installing DESeq2 package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")

library(apeglm)

# Load control counts 
control_counts <- read.csv("/project/bf528/project_3/samples/control_counts.csv") %>%
  select(Geneid, as.character(sample_info_control[[1]]))

# Merge controls with count matrices
merged_controls <- merge(merged_gid, control_counts, by = "Geneid") %>% column_to_rownames(var="Geneid")

# Writing CSV for merged controls
#write.csv(merged_controls, file="featureCountsandcontrol_combined.csv", row.names=FALSE)

# filter out rows that have any zeros for funzies
merged_filtered <- subset(merged_controls,rowSums(merged_controls==0)==0)


# Splitting up per MOA in sample_info
dna <- sample_info %>% filter(mode_of_action == "DNA_Damage" | mode_of_action == "Control" & vehicle == "SALINE_100_%") %>% select(Run)
er <- sample_info %>% filter(mode_of_action == "ER" | mode_of_action == "Control" & vehicle == "CORN_OIL_100_%") %>% select(Run)
ppara <- sample_info %>% filter(mode_of_action == "PPARA" | mode_of_action == "Control" & vehicle == "CORN_OIL_100_%") %>% select(Run)

# Creating MOA DESeq objects
ds_dna <- DESeqDataSetFromMatrix(
  countData = merged_filtered %>% select(dna[[1]]),
  colData = sample_info[c(1:3, 9:11),],
  design = ~ mode_of_action
)

# MOA: ER, adjust for unattainable counts file SRR1178023
ds_er <- DESeqDataSetFromMatrix(
  countData = merged_filtered %>% select(er[[1]]),
  colData = sample_info[c(4:5, 12:14),],
  design = ~ mode_of_action
)

ds_ppara <- DESeqDataSetFromMatrix(
  countData = merged_filtered %>% select(ppara[[1]]),
  colData = sample_info[c(6:8, 12:14),],
  design = ~ mode_of_action
)

# Relevel mode_of_action as factor for each DESeq object
ds_dna$mode_of_action <- relevel(ds_dna$mode_of_action, ref="Control")
ds_er$mode_of_action <- relevel(ds_er$mode_of_action, ref="Control")
ds_ppara$mode_of_action <- relevel(ds_ppara$mode_of_action, ref="Control")

# Running DESeq on objects
ds_dna <- DESeq(ds_dna)
ds_er <- DESeq(ds_er)
ds_ppara <- DESeq(ds_ppara)

# Storing DEseq results
dna_results <- results(ds_dna, contrast=c('mode_of_action','DNA_Damage','Control'))
dna_results <- lfcShrink(ds_dna, coef=2)

er_results <- results(ds_er, contrast=c('mode_of_action','ER','Control'))
er_results <- lfcShrink(ds_er, coef=2)

ppara_results <- results(ds_ppara, contrast=c('mode_of_action','PPARA','Control'))
ppara_results <- lfcShrink(ds_ppara, coef=2)

# Writing results to CSV files
write.csv(dna_results,"deseq/dna_deseq.csv")
write.csv(er_results,"deseq/er_deseq.csv")
write.csv(ppara_results,"deseq/ppara_deseq.csv")

# Storing counts for each file
dna_counts <- read.csv("/projectnb/bf528/users/tinman/Project3/programmer/deseq/dna_deseq.csv", header = TRUE, row.names = 1)
er_counts <- read.csv("/projectnb/bf528/users/tinman/Project3/programmer/deseq/er_deseq.csv", header = TRUE, row.names = 1)
ppara_counts <- read.csv("/projectnb/bf528/users/tinman/Project3/programmer/deseq/ppara_deseq.csv", header = TRUE, row.names = 1)

# Number of genes significant at p-adj < 0.05
#DNA = 4384 genes
dna_counts %>% filter(padj < 0.05) %>% summarize(count = n())

#ER = 2627 genes
er_counts %>% filter(padj < 0.05) %>% summarize(count = n())

#PPARA = 2838 genes
ppara_counts %>% filter(padj < 0.05) %>% summarize(count = n())

# Creating data frame of significant gene counts
sig_genes <- data.frame(mode_of_action=c('DNA Damage','Estrogen Receptor (ER)', 'Peroxisome Proliferator-Activated Receptor Alpha (PPARA)'),
  gene_count=c('4,384','2,627','2,838'))

# Writing gene counts to CSV
write.csv(sig_genes,"deliverables_part4/number_significant_genes.csv")

# Top 10 diff expressed genes by p-value
dna_counts %>% top_n(-10, pvalue) %>% write.csv("deseq/top10DE/dna.csv", row.names = TRUE)
er_counts %>% top_n(-10, pvalue) %>% write.csv("deseq/top10DE/er.csv", row.names = TRUE)
ppara_counts %>% top_n(-10, pvalue) %>% write.csv("deseq/top10DE/ppara.csv", row.names = TRUE)

# Histograms of FC Values
#DNA
png("deseq/figures/histo_dna.png")
dna_counts %>% filter(padj < 0.05) %>%
  ggplot(mapping = aes(x = log2FoldChange)) +
  geom_histogram() +
  labs(title = "DNA Damage Fold Change Histogram", x = "Log2 Fold Change", y = "Frequency")
dev.off()

#ER
png("deseq/figures/histo_er.png")
er_counts %>% filter(padj < 0.05) %>%
  ggplot(mapping = aes(x = log2FoldChange)) +
  geom_histogram() +
  labs(title = "Estrogen Receptor (ER) Fold Change Histogram", x = "Log2 Fold Change", y = "Frequency")
dev.off()

#PPARA
png("deseq/figures/histo_ppara.png")
ppara_counts %>% filter(padj < 0.05) %>%
  ggplot(mapping = aes(x = log2FoldChange)) +
  geom_histogram() +
  labs(title = "Peroxisome Proliferator-Activated Receptor Alpha (PPARA) Fold Change Histogram", x = "Log2 Fold Change", y = "Frequency")
dev.off()

# Reloading counts into data frames
df_dna <- as.data.frame(dna_counts)
row.names(df_dna) <- df_dna[, 1]
df_dna <- df_dna[,-1]
ds = df_dna[order(df_dna$log2FoldChange),]

df_er <- as.data.frame(er_counts)
row.names(df_er) <- df_dna[, 1]
df_er <- df_er[,-1]
es = df_er[order(df_er$log2FoldChange),]

df_ppara <- as.data.frame(ppara_counts)
row.names(df_ppara) <- df_ppara[, 1]
df_ppara <- df_ppara[,-1]
ps = df_ppara[order(df_ppara$log2FoldChange),]

# Scatter plots
#DNA
ggplot(df_dna, aes(log2FoldChange, -log(pvalue))) + 
  geom_point() +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=tail(ds, 10), col="red") +
  geom_point(data=head(ds, 10), col="dodgerblue2") +
  ggtitle("DNA Damage Fold Change Scatter Plot")
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))

#ER
ggplot(df_er, aes(log2FoldChange, -log(pvalue))) + 
  geom_point() +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=tail(es, 10), col="red") +
  geom_point(data=head(es, 10), col="dodgerblue2") +
  ggtitle("Estrogen Receptor (ER) Fold Change Scatter Plot")
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))
  
#PPARA
ggplot(df_ppara, aes(log2FoldChange, -log(pvalue))) + 
  geom_point() +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=tail(ps, 10), col="red") +
  geom_point(data=head(ps, 10), col="dodgerblue2") +
  ggtitle("Peroxisome Proliferator-Activated Receptor Alpha (PPARA) Fold Change Scatter Plot")
  theme_classic() +
  theme(axis.title.x = element_text(size = 13, face="bold"),
        axis.title.y = element_text(size = 13, face="bold"))

# write out matrix of normalized counts
write.csv(counts(ds_dna,normalized=TRUE),'dna_norm_deseq.csv')
write.csv(counts(ds_er,normalized=TRUE),'er_norm_deseq.csv')
write.csv(counts(ds_ppara,normalized=TRUE),'ppara_norm_deseq.csv')
