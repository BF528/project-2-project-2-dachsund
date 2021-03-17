

library(dplyr)
library(tibble)
library(ggplot2)


setwd("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/")
fpkm <- read.csv('fpkm_matrix.csv', header=TRUE, sep="\t")
genediff <- read.table("/projectnb2/bf528/project_2/data/P4_vs_P7_cuffdiff_out/gene_exp.diff", header=TRUE, sep= "\t")
david_down5 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/down_regulated_genes.csv")
david_top5 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/up_regulated_genes.csv")


#Reading data from fpkm tracking tables
P0_1 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/programmer/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE, sep="\t")

P0_2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P0_2_genes.fpkm_tracking", header=TRUE, sep="\t")

P4_1 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P4_1_genes.fpkm_tracking", header=TRUE, sep="\t")

P4_2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P4_2_genes.fpkm_tracking", header=TRUE, sep="\t")

P7_1 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P7_1_genes.fpkm_tracking", header=TRUE, sep="\t")

P7_2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/P7_2_genes.fpkm_tracking", header=TRUE, sep="\t")

Ad1 <-read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/Ad_1_genes.fpkm_tracking", header=TRUE, sep="\t")

Ad2 <- read.csv("/projectnb2/bf528/users/dachshund/project_2/project-2-project-2-dachsund/biologist/Ad_2_genes.fpkm_tracking", header=TRUE, sep="\t")

#Filtering data

#Gene names identified from Figure 1D for , Mitochondria, Cell Cycle 
sarcomere_genes <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab")
mitochondria_genes <- c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
cell_cycle_genes <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")

#Filtering step genes -- SACROMERE 
#Adult filter
sarcomere_adult1 <- Ad1 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
sarcomere_adult2 <- Ad2 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

s_adult_combined <- cbind(sarcomere_adult1, sarcomere_adult1$FPKM)
s_adult_combined$mean <-rowMeans(s_adult_combined)
s_adult_mean <- s_adult_combined %>% rownames_to_column('Genes')

adult_sarcomere_genes <- data.frame('Genes' = s_adult_mean$Genes, 'FPKM' = s_adult_mean$mean, 'Status' = c('Adult', 'Adult', 'Adult', 'Adult', 'Adult', 'Adult', 'Adult'))

#P0 filter
sarcomere_P0_1 <- P0_1 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
sarcomere_P0_2 <- P0_2 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

s_P0_combined <- cbind(sarcomere_P0_1, sarcomere_P0_2$FPKM)
s_P0_combined$mean <-rowMeans(s_P0_combined)
s_P0_mean <- s_P0_combined %>% rownames_to_column('Genes')

P0_sarcomere_genes <- data.frame('Genes' = s_P0_mean$Genes, 'FPKM' = s_P0_mean$mean, 'Status' = c('P0', 'P0', 'P0', 'P0', 'P0', 'P0', 'P0'))

#P4 filter
sarcomere_P4_1 <- P4_1 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
sarcomere_P4_2 <- P4_2 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

s_P4_combined <- cbind(sarcomere_P4_1, sarcomere_P4_2$FPKM)
s_P4_combined$mean <-rowMeans(s_P4_combined)
s_P4_mean <- s_P4_combined %>% rownames_to_column('Genes')

P4_sarcomere_genes <- data.frame('Genes' = s_P4_mean$Genes, 'FPKM' = s_P4_mean$mean, 'Status' = c('P4', 'P4', 'P4', 'P4', 'P4', 'P4', 'P4'))

#P7 filter
sarcomere_P7_1 <- P7_1 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
sarcomere_P7_2 <- P7_2 %>% filter(gene_short_name %in% sarcomere_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

s_P7_combined <- cbind(sarcomere_P7_1, sarcomere_P7_2$FPKM)
s_P7_combined$mean <-rowMeans(s_P7_combined)
s_P7_mean <- s_P7_combined %>% rownames_to_column('Genes')

P7_sarcomere_genes <- data.frame('Genes' = s_P7_mean$Genes, 'FPKM' = s_P7_mean$mean, 'Status' = c('P7', 'P7', 'P7', 'P7', 'P7', 'P7', 'P7'))

finished_sarcomere_results <- bind_rows(P0_sarcomere_genes,P4_sarcomere_genes,P7_sarcomere_genes, adult_sarcomere_genes)

#Filtering step genes -- MITOCHONDRIA

#Adult filter
mitochondria_adult1 <- Ad1 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
mitochondria_adult2 <- Ad2 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

m_adult_combined <- cbind(mitochondria_adult1, mitochondria_adult1$FPKM)
m_adult_combined$mean <-rowMeans(m_adult_combined)
m_adult_mean <- m_adult_combined %>% rownames_to_column('Genes')

adult_mitochondria_genes <- data.frame('Genes' = m_adult_mean$Genes, 'FPKM' = m_adult_mean$mean, 'Status' = c('Adult', 'Adult', 'Adult', 'Adult', 'Adult'))

#P0 filter
mitochondria_P0_1 <- P0_1 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
mitochondria_P0_2 <- P0_2 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

m_P0_combined <- cbind(mitochondria_P0_1, mitochondria_P0_2$FPKM)
m_P0_combined$mean <-rowMeans(m_P0_combined)
m_P0_mean <- m_P0_combined %>% rownames_to_column('Genes')

P0_mitochondria_genes <- data.frame('Genes' = m_P0_mean$Genes, 'FPKM' = m_P0_mean$mean, 'Status' = c('P0', 'P0', 'P0', 'P0', 'P0'))

#P4 filter
mitochondria_P4_1 <- P4_1 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
mitochondria_P4_2 <- P4_2 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

m_P4_combined <- cbind(mitochondria_P4_1, mitochondria_P4_2$FPKM)
m_P4_combined$mean <-rowMeans(m_P4_combined)
m_P4_mean <- m_P4_combined %>% rownames_to_column('Genes')

P4_mitochondria_genes <- data.frame('Genes' = m_P4_mean$Genes, 'FPKM' = m_P4_mean$mean, 'Status' = c('P4', 'P4', 'P4', 'P4', 'P4'))

#P7 filter
mitochondria_P7_1 <- P7_1 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
mitochondria_P7_2 <- P7_2 %>% filter(gene_short_name %in% mitochondria_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

m_P7_combined <- cbind(mitochondria_P7_1, mitochondria_P7_2$FPKM)
m_P7_combined$mean <-rowMeans(m_P7_combined)
m_P7_mean <- m_P7_combined %>% rownames_to_column('Genes')

P7_mitochondria_genes <- data.frame('Genes' = m_P7_mean$Genes, 'FPKM' = m_P7_mean$mean, 'Status' = c('P7', 'P7', 'P7', 'P7', 'P7'))

finished_mitochondria_results <- bind_rows(P0_mitochondria_genes,P4_mitochondria_genes,P7_mitochondria_genes, adult_mitochondria_genes)

#Filtering step genes -- CELL CYCLE

#Adult filter
cellcycle_adult1 <- Ad1 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
cellcycle_adult2 <- Ad2 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

cc_adult_combined <- cbind(cellcycle_adult1, cellcycle_adult1$FPKM)
cc_adult_combined$mean <-rowMeans(cc_adult_combined)
cc_adult_mean <- cc_adult_combined %>% rownames_to_column('Genes')

adult_cellcycle_genes <- data.frame('Genes' = cc_adult_mean$Genes, 'FPKM' = cc_adult_mean$mean, 'Status' = c('Adult', 'Adult', 'Adult', 'Adult', 'Adult','Adult', 'Adult', 'Adult', 'Adult', 'Adult'))

#P0 filter
cellcycle_P0_1 <- P0_1 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
cellcycle_P0_2 <- P0_2 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

cc_P0_combined <- cbind(cellcycle_P0_1, cellcycle_P0_2$FPKM)
cc_P0_combined$mean <-rowMeans(cc_P0_combined)
cc_P0_mean <- cc_P0_combined %>% rownames_to_column('Genes')

P0_cellcycle_genes <- data.frame('Genes' = cc_P0_mean$Genes, 'FPKM' = cc_P0_mean$mean, 'Status' = c('P0', 'P0', 'P0', 'P0', 'P0','P0', 'P0', 'P0', 'P0', 'P0'))

#P4 filter
cellcycle_P4_1 <- P4_1 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
cellcycle_P4_2 <- P4_2 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

cc_P4_combined <- cbind(cellcycle_P4_1, cellcycle_P4_2$FPKM)
cc_P4_combined$mean <-rowMeans(cc_P4_combined)
cc_P4_mean <- cc_P4_combined %>% rownames_to_column('Genes')

P4_cellcycle_genes <- data.frame('Genes' = cc_P4_mean$Genes, 'FPKM' = cc_P4_mean$mean, 'Status' = c('P4', 'P4', 'P4', 'P4', 'P4','P4', 'P4', 'P4', 'P4', 'P4'))

#P7 filter
cellcycle_P7_1 <- P7_1 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 
cellcycle_P7_2 <- P7_2 %>% filter(gene_short_name %in% cell_cycle_genes) %>% select(gene_short_name, FPKM) %>% column_to_rownames('gene_short_name') 

cc_P7_combined <- cbind(cellcycle_P7_1, cellcycle_P7_2$FPKM)
cc_P7_combined$mean <-rowMeans(cc_P7_combined)
cc_P7_mean <- cc_P7_combined %>% rownames_to_column('Genes')

P7_cellcycle_genes <- data.frame('Genes' = cc_P7_mean$Genes, 'FPKM' = cc_P7_mean$mean, 'Status' = c('P7', 'P7', 'P7', 'P7', 'P7','P7', 'P7', 'P7', 'P7', 'P7'))

finished_cellcycle_results <- bind_rows(P0_cellcycle_genes,P4_cellcycle_genes,P7_cellcycle_genes, adult_cellcycle_genes)

#Plotting data

#Sacromere plot
finished_sarcomere_results$Status = as.factor(finished_sarcomere_results$Status)
finished_sarcomere_results$Status <- factor(finished_sarcomere_results$Status, levels=c('P0','P4','P7','Adult'))
ggplot(finished_sarcomere_results, aes(x=Status, y=FPKM, col=Genes, group=Genes)) +geom_point() + geom_path()+labs(x='', y="FPKM")+scale_y_continuous(breaks=seq(100,1300,100))+ggtitle('Sarcomere')+theme(plot.title = element_text(hjust=0.5))

#Mitochondria plot
finished_mitochondria_results$Status = as.factor(finished_mitochondria_results$Status)
finished_mitochondria_results$Status <- factor(finished_mitochondria_results$Status, levels=c('P0','P4','P7','Adult'))
ggplot(finished_mitochondria_results, aes(x=Status, y=FPKM, col=Genes, group=Genes)) +geom_point() + geom_path()+scale_y_continuous(breaks=seq(0,280,50))+ggtitle('Mitochondria')+theme(plot.title = element_text(hjust=0.5))

#Cell Cycle plot
finished_cellcycle_results$Status = as.factor(finished_cellcycle_results$Status)
finished_cellcycle_results$Status <- factor(finished_cellcycle_results$Status, levels=c('P0','P4','P7','Adult'))
ggplot(finished_cellcycle_results, aes(x=Status, y=FPKM, col=Genes, group=Genes)) +geom_point() + geom_path()+labs(x='', y="FPKM")+scale_y_continuous(breaks=seq(0,80,10))+ggtitle('Cell Cycle')+theme(plot.title = element_text(hjust=0.75))