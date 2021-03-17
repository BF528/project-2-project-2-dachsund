
library(dplyr)
setwd("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/")
fpkm <- read.csv('fpkm_matrix.csv', header=TRUE, sep="\t")
genediff <- read.table("/projectnb2/bf528/project_2/data/P4_vs_P7_cuffdiff_out/gene_exp.diff", header=TRUE, sep= "\t")

#Reading data from fpkm tracking tables
P0_1 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/programmer/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE, sep="\t")

P0_2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P0_2_genes.fpkm_tracking", header=TRUE, sep="\t")

P4_1 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P4_1_genes.fpkm_tracking", header=TRUE, sep="\t")

P4_2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P4_2_genes.fpkm_tracking", header=TRUE, sep="\t")

P7_1 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P7_1_genes.fpkm_tracking", header=TRUE, sep="\t")

P7_2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P7_2_genes.fpkm_tracking", header=TRUE, sep="\t")

Ad1 <-read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/Ad_1_genes.fpkm_tracking", header=TRUE, sep="\t")

Ad2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/Ad_2_genes.fpkm_tracking", header=TRUE, sep="\t")





#Clustered heatmap of at MOST 1000 DE genes dound in P0 vs Adult
#I will be only looking at top 100
Gene_diff_2 <-genediff %>% filter(genediff$significant == 'yes' & genediff$significant == 'yes' & genediff$log2.fold_change. != 'Inf' & genediff$log2.fold_change. != '-Inf') %>% arrange(log2.fold_change.)

head50 <- Gene_diff_2 %>% head(Gene_diff_2[unique(Gene_diff_2$gene),]$log2.fold_change.,n=50)
tail50 <- Gene_diff_2 %>% head(Gene_diff_2[unique(genediff$gene),]$log2.fold_change.,n=50)
diff100 <- rbind(head50,tail50)
diff100_data <- data.frame('Gene'=diff100$gene,'log fold change' = diff100$log2.fold_change.)


#FPKM valuess


library(tidyverse)
FPKM_P0_1 <- P0_1 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM)

FPKM_P0_2 <- P0_2 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM) 

FPKM_P4_1 <- P4_1 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM)

FPKM_P4_2 <- P4_2 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM) 

FPKM_P7_1 <- P7_1 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM) 

FPKM_P7_2 <- P7_2 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM)

FPKM_Adult_1 <- Ad1 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM)

FPKM_Adult_2 <- Ad2 %>% filter(gene_short_name %in% diff100_data$Gene) %>% select(gene_short_name,FPKM) 

P01_FPKM_2 <- FPKM_P0_1[!duplicated(FPKM_P0_1$gene_short_name), ] %>% rename(P0_1 = FPKM)
P02_FPKM_2 <- FPKM_P0_2[!duplicated(FPKM_P0_2$gene_short_name), ] %>% rename(P0_2 = FPKM)
P41_FPKM_2 <- FPKM_P4_1[!duplicated(FPKM_P4_1$gene_short_name), ] %>% rename(P4_1 = FPKM)
P42_FPKM_2 <- FPKM_P4_2[!duplicated(FPKM_P4_2$gene_short_name), ] %>% rename(P4_2 = FPKM)
P71_FPKM_2 <- FPKM_P7_1[!duplicated(FPKM_P7_1$gene_short_name), ] %>% rename(P7_1 = FPKM)
P72_FPKM_2 <- FPKM_P7_2[!duplicated(FPKM_P7_2$gene_short_name), ] %>% rename(P7_2 = FPKM)
Adult1_FPKM_2 <- FPKM_Adult_1[!duplicated(FPKM_Adult_1$gene_short_name), ] %>% rename(Ad1 = FPKM)
Adult2_FPKM_2 <- FPKM_Adult_2[!duplicated(FPKM_Adult_2$gene_short_name), ] %>% rename(Ad2 = FPKM)

concat <- inner_join(P01_FPKM_2,P02_FPKM_2, by='gene_short_name') %>% inner_join(., P41_FPKM_2, by='gene_short_name') %>% inner_join(., P42_FPKM_2, by='gene_short_name') %>% inner_join(., P71_FPKM_2, by='gene_short_name') %>% inner_join(., P72_FPKM_2, by='gene_short_name') %>% inner_join(., Adult1_FPKM_2, by='gene_short_name') %>% inner_join(., Adult2_FPKM_2, by='gene_short_name')

concat2 <- concat[!duplicated(concat$gene_short_name),]
rownames(concat2) <- NULL
finished_concat <- concat2 %>% column_to_rownames('gene_short_name')

heatmap(as.matrix(finished_concat))






