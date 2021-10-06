# Installed hgu133plus2.db 
# Installed affy for exploration and analysis of oligonucleotide
# Installed sva to remove batch effects and unwanted variation
# Installed GSEAbase to load gene sets
# Installed AnnotationDbi to query data packages

# Load packages
library(affy)
library(hgu133plus2.db)
library(sva)
library(AnnotationDbi)
library(GSEABase)
library(dplyr)
library(tidyverse)

# Checking working directory
getwd()

# Read differential expression file and organize by descending t-stat trend
de_data <- read.csv("/projectnb/bf528/project_1/data/differential_expression_results.csv", row.names = 1, header = TRUE)
View(de_data)
de_data <- de_data %>% arrange(desc(t))

# Map probeset IDs to appropriate gene symbol
de_matches <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(row.names(de_data)), columns = ("SYMBOL"))

# Remove duplicated matches
de_duplicatedmatches <- de_matches[!duplicated(de_matches[1]),]

# Combine symbols with initial de_data to add symbol column
de_data <- cbind(de_duplicatedmatches, de_data)

# Keeping only the probes with significant padj (adjust p-values)
de_data <- de_data %>%
  group_by(SYMBOL) %>%
  filter(padj == min(padj)) %>%
  ungroup(SYMBOL)

# KEGG, GO, and Hallmark genesets downloaded
# Loading gene sets
hallmarks <- getGmt("/projectnb/bf528/users/Tinman/Project_1/Biologist Files/h.all.v7.2.symbols.gmt")
kegg <- getGmt("/projectnb/bf528/users/Tinman/Project_1/Biologist Files/c2.cp.kegg.v7.2.symbols.gmt")
go <- getGmt("/projectnb/bf528/users/Tinman/Project_1/Biologist Files/c5.go.v7.2.symbols.gmt")

# Number of gene sets per collection
cat("Number of gene sets in Hallmark: ", length(names(hallmarks)))
cat("Number of gene sets in GO: ", length(names(go)))
cat("Number of gene sets in KEGG: ", length(names(kegg)))

# Remove missing values
de_data <- de_data[!is.na(de_data$SYMBOL),]

# Gather top 1000 up-regulated and down-regulated genes 
top1000_upreg <- head(de_data, 1000)
top1000_downreg <- tail(de_data, 1000)

# Gather top 10 up-regulated and down-regulated genes from top 1000 list
top10_upreg <- head(top1000_upreg, 10)
top10_downreg <- tail(top1000_downreg, 10)

# Install and load packages to create a table
library(knitr)
library(kableExtra)

# Create tables for both top 10 up- and down-regulated genes
kable(top10_upreg, digits = 50)
kable(top10_downreg, digits = 50)

# Save as csv
write.csv(top10_upreg, file = "top10upreg.csv", row.names = FALSE)
write.csv(top10_downreg, file = "top10downreg.csv", row.names = FALSE)

# Obtain genes which were not expressed
up_nde <- subset(de_data, !de_data$SYMBOL %in% top1000_upreg$SYMBOL)
down_nde <- subset(de_data, !de_data$SYMBOL %in% top1000_downreg$SYMBOL)

# Defining fisher test function with gl = genelist, gs= geneset, nde= not differentially expressed
fishertest <- function(gl, gs, nde)           
{ diffexp_ings <- length(intersect(gl,gs))    #diffexp_ings = differentially expressed genes that are present in geneset
diffexp_notgs <- length(gl) - diffexp_ings    #diffexp_notgs = differentially expressed genes that are not present in geneset 
notde_ings <- length(intersect(nde,gs))       #notde_ings = not expressed but present in geneset
notde_notgs <- length(nde) - notde_ings       #notde_notgs = not differentially expressed and not in geneset
return(c(diffexp_ings,diffexp_notgs,notde_ings,notde_notgs))}   #returns the fishertest values
values


# Initialize data frame to store fisher resultsf for each gene set
kegg_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)
go_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)
hallmark_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

# Using a for loop to obtain fisher test values and input into dataframe
# Kegg
for (i in 1:length(kegg))
{
  geneid <- geneIds(kegg[i])
  fisher_up <- fishertest(top1000_upreg$SYMBOL, geneid[[names(geneid)]], up_nde$SYMBOL)
  fisher_down <- fishertest(top1000_downreg$SYMBOL, geneid[[names(geneid)]], down_nde$SYMBOL)
  up_reg <- fisher.test(matrix(fisher_up,nrow=2))
  down_reg <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid), up_reg$p.value, up_reg$estimate, 'Up')
  kegg_results[nrow(kegg_results) +1, ] <- c(names(geneid), down_reg$p.value, down_reg$estimate, 'Down')}

kegg_results <- kegg_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))

# GO
for (i in 1:length(go))
{
  geneid <- geneIds(go[i])
  fisher_up <- fishertest(top1000_upreg$SYMBOL, geneid[[names(geneid)]], up_nde$SYMBOL)
  fisher_down <- fishertest(top1000_downreg$SYMBOL, geneid[[names(geneid)]], down_nde$SYMBOL)
  up_reg <- fisher.test(matrix(fisher_up,nrow=2))
  down_reg <- fisher.test(matrix(fisher_down, nrow=2))
  go_results[nrow(go_results) +1, ] <- c(names(geneid), up_reg$p.value, up_reg$estimate, 'Up')
  go_results[nrow(go_results) +1, ] <- c(names(geneid), down_reg$p.value, down_reg$estimate, 'Down')}

go_results <- go_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))

# Hallmark
for (i in 1:length(hallmarks))
{
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(top1000_upreg$SYMBOL, geneid[[names(geneid)]], up_nde$SYMBOL)
  fisher_down <- fishertest(top1000_downreg$SYMBOL, geneid[[names(geneid)]], down_nde$SYMBOL)
  up_reg <- fisher.test(matrix(fisher_up,nrow=2))
  down_reg <- fisher.test(matrix(fisher_down, nrow=2))
  hallmark_results[nrow(hallmark_results) +1, ] <- c(names(geneid), up_reg$p.value, up_reg$estimate, 'Up')
  hallmark_results[nrow(hallmark_results) +1, ] <- c(names(geneid), down_reg$p.value, down_reg$estimate, 'Down')}

hallmark_results <- hallmark_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))

# Adjusting p-value to FDR method
kegg_results$BH <- p.adjust(kegg_results$pvalue, method = "BH", n = length(kegg_results$pvalue))
write.csv(kegg_results, "kegg.csv")

go_results$BH <- p.adjust(go_results$pvalue, method = "BH", n = length(go_results$pvalue))
write.csv(go_results, "go.csv")     

hallmark_results$BH <- p.adjust(hallmark_results$pvalue, method = "BH", n = length(hallmark_results$pvalue))
write.csv(hallmark_results, "hallmarks.csv")

# Statistically significant and enriched genesets
# Kegg
kegg_enr <- kegg_results[kegg_results$pvalue<0.05,]
kegg_enr_num <- length(kegg_enr$gene_set)
cat("Statistically enriched genesets in KEGG: ", kegg_enr_num, "\n")

# GO:
go_enr <- go_results[go_results$pvalue<0.05,]
go_enr_num <- length(go_enr$gene_set)
cat("Statistically enriched genesets in GO: ", go_enr_num, "\n")

# Hallmark:
hm_enr <- hallmark_results[hallmark_results$pvalue<0.05,]
hm_enr_num <- length(hm_enr$gene_set)
cat("Statistically enriched genesets in Hallmark: ", hm_enr_num, "\n")

# Top 3 genesets for each:
top3_kegg <- slice_min(kegg_results, order_by=pvalue, n=3)
top3_go <- slice_min(go_results, order_by=pvalue, n=3)
top3_hm <- slice_min(hallmark_results, order_by=pvalue, n=3)
top3s <- rbind(top3_kegg, top3_go, top3_hm)
print(top3s)

write.csv(top3s, file="geneset_top_results.csv")


#############################################
# Creating Fisher test function for p-values
fisher_pvals <- function(gmt){
  pvalues <- c()
  df <- list()
  for(geneset in gmt){
    setname <- setName(geneset)
    geneids <- geneIds(geneset)
    differentially_expressed <- length(de_data$SYMBOL)
    in_set <- length(geneids)
    in_set_differential <- sum(geneids %in% de_data$SYMBOL)
    in_set_not_differential <- in_set - in_set_differential
    not_in_set_differential <- differentially_expressed - in_set_differential
    not_in_set_not_differential <- 0
    fishervals <- fisher.test(matrix(c(in_set_differential, in_set_not_differential, not_in_set_differential, not_in_set_not_differential), nrow = 2))
    pval <- fishervals$p.value
    pvalues[setname] <- pval
  }
  return(pvalues)
}

# Creating function to generate dataframe of fisher test values
fisher_df <- function(gmt){
  df <- list()
  for(geneset in gmt){
    setname <- setName(geneset)
    geneids <- geneIds(geneset)
    differentially_expressed <- length(de_data$SYMBOL)
    in_set <- length(geneids)
    in_set_differential <- sum(geneids %in% de_data$SYMBOL)
    in_set_not_differential <- in_set - in_set_differential
    not_in_set_differential <- differentially_expressed - in_set_differential
    not_in_set_not_differential <- 0
    fishervals <- fisher.test(matrix(c(in_set_differential, in_set_not_differential, not_in_set_differential, not_in_set_not_differential), nrow = 2))
    pval <- fishervals$p.value
    est <- fishervals$estimate
    padj <- p.adjust(pval, method="fdr")
    df[[setname]] <- data.frame(geneset = setname, statistic = est, pvalue = pval, p.adj = padj)
  }
  return(df)
}

# KEGG Pathways
pvalues_kegg <- fisher_pvals(kegg)
df_kegg <- fisher_df(kegg)

# GO Pathways
pvalues_go <- fisher_pvals(go)
df_go <- fisher_df(go)

# Hallmark Pathways
pvalues_hallmarks <- fisher_pvals(hallmarks)
df_hallmarks <- fisher_df(hallmarks)

#to do - Enriched gene sets at p<0.05 and top 3 enriched sets

hallmarks_enr <- df_hallmarks[df_hallmarks$padj<0.05, 3]
hallmarks_enr_total <- length(hallmarks_enr$geneset)
cat("Statistically enriched genesets in KEGG: ", sig_kegg_num, "\n")

kegg_top3 <- names(head(sort(pvalues_kegg), 3))
keggs <- rbind(df_kegg[[kegg_top3[[1]]]], df_kegg[[kegg_top3[[2]]]], df_kegg[[kegg_top3[[3]]]])

# write.csv(keggs, "top_kegg.csv")
