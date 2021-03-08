# title: project 2 analyst
# author: Sheila Yee

# load differential expression file analysis
diff_expr_file <- read.table("/projectnb/bf528/users/dachshund/project_2/project-2-project-2-dachsund/programmer/cuffdiff_out/gene_exp.diff",
                             header = TRUE)

# sort the above data table so that the smallest q_values are at the top 
dif_expr_file_sorted <- diff_expr_file[order(diff_expr_file$q_value),]

# create a table with top ten differentially expressed genes, with their names, 
# FPKM values, log fold change, p-value, and q-value
dif_expr_top_10 <- dif_expr_file_sorted[1:10,]  
dif_expr_top_10_final <- dif_expr_top_10[c("gene", "value_1", "value_2", "log2.fold_change.", "p_value", "q_value")]

# produce a histogram of the log2.foldchange column for all genes
log2change_hist <- hist(diff_expr_file$log2.fold_change., 
                        breaks = 70, 
                        main = NULL, 
                        ylim = c(0, 20000),
                        xlim = c(-10, 10),
                        xlab = "Log2 fold change", 
                        col = c("lightblue")
                        )

# create a new data frame that contains only the genes where the last column named significant is equal to yes
dif_expr_sig <- subset(diff_expr_file, significant == "yes")

# create a second histogram of the log2 fold change values only for significant genes
dif_expr_sig_hist <- hist(dif_expr_sig$log2.fold_change.,
                          breaks = 70, 
                          main = NULL, 
                          ylim = c(0, 500),
                          xlim = c(-10, 10),
                          xlab = "Log2 fold change", 
                          col = c("darkred")
                          )

# create separate dataframes for only significant genes that had positive and negative log fold change
dif_expr_sig_pos <-subset(dif_expr_sig,log2.fold_change. > 0 )
dif_expr_sig_neg <-subset(dif_expr_sig,log2.fold_change. < 0 )
dim(dif_expr_sig_pos) # find out how many genes with positive log fold change
dim(dif_expr_sig_neg) # find out how many genes with negative log fold change

# write the up and down regulated gene names to separate different files
up_regulated_genes <- write.csv(dif_expr_sig_pos$gene, "up_regulated_genes.csv")
down_regulated_genes <- write.csv(dif_expr_sig_neg$gene, "down_regulated_genes.csv")
 


  

