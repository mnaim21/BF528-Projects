# Loading libraries
library('ggplot2')
library('ggpubr')

# Setting working directory
setwd("/projectnb2/bf528/users/tinman/Project4/project5/maha")

# Loading in expression file and selecting 10 most significant genes
gene_diff <- read.delim("gene_exp.diff")
gene_diff <- gene_diff[order(gene_diff$q_value,decreasing = FALSE),]
top_ten_genes <- gene_diff[1:10,c(3,8,9,10,12,13)]
write.table(top_ten_genes,'top_ten_diffexp_genes.txt',sep = '\t')

# Histogram of log2 fold change
ggplot(gene_diff, aes(x = log2.fold_change.)) + geom_histogram(binwidth = 0.3, fill = 'dodgerblue3', color = 'black') + 
  theme_light() + xlab('Log2 Fold Change') + ylab('Count') + scale_y_log10() +
  ggtitle('All Genes') + theme(plot.title = element_text(size=11)) 


# New dataframe for significant genes
sig_genes <- subset(gene_diff,significant == 'yes')

# Histogram for significant genes
ggplot(sig_genes, aes(x = log2.fold_change.)) + geom_histogram(binwidth = 0.3, fill = 'darkorange1', color = 'black') + 
  theme_light() + xlab('Log2 Fold Change') + ylab('Count') + scale_y_log10() +
  ggtitle('Significant Genes')+ theme(plot.title = element_text(size=11))

# PDF containing combined histograms
pdf('DiffExp_Hist.pdf',height = 4,width = 6)
ggplot() + 
  geom_histogram(data = gene_diff, aes(x = log2.fold_change., color = 'Non-Significant'), binwidth = 0.3, fill = 'grey85') +
  geom_histogram(data = sig_genes, aes(x = log2.fold_change., color = 'Significant'),binwidth = 0.3, fill = 'grey85', size = 0.75) + 
  scale_y_log10() + xlab('Log2 Fold Change') + ylab('Count') + theme_light() + 
  ggtitle('Distribution of Fold Change in Differentailly Expressed Genes\nfor Natal Versus Adult Cardiac Myocytes') +
  scale_color_manual(values = c('Non-Significant' = 'blue', 'Significant' = 'red'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = c(0.85,0.85),
        plot.title = element_text(hjust = 0.5))
dev.off()

# Subsetting list into upregulated and downregulated
up_reg <- sig_genes[sig_genes$log2.fold_change. > 0,] #242 up regulated genes
down_reg <- sig_genes[sig_genes$log2.fold_change. < 0,] #215 down regulated genes

# Writing gene list to file for DAVID analysis
write(up_reg$gene, file = 'UpReg_genes.txt')
write(down_reg$gene, file = 'DownReg_genes.txt')