#Sets working directory; Reads table
setwd("/projectnb/bf528/users/dachshund/project_2/project-2-project-2-dachsund/programmer/")
cd <- read.table("P0_1_cufflinks/genes.fpkm_tracking")

colnames(cd) <- cd[1,]#Set column names 
cd <- cd [-1,] #delete row of column names

#Prepare the data
cd$FPKM <- as.numeric(cd$FPKM) #Change to numeric
#Create numeric lists; to.eX means to 10*X power.
to.e1 <-cd[cd$FPKM > 0 & cd$FPKM <= 10,]$FPKM
to.e2<-cd[cd$FPKM > 10 & cd$FPKM <= 100,]$FPKM
to.e3 <-cd[cd$FPKM > 100 & cd$FPKM <= 1000,]$FPKM
to.e4 <-cd[cd$FPKM > 1000,]$FPKM
#create pdf with histograms
pdf('cufflink.density.pdf')
hist(cd$FPKM, main = "ALL FPKM", xlab="FPKM", col = 'green')
par(mfrow=c(2,2))
hist(to.e1, main = "A) FPKM less than or equal to 10", xlab="FPKM", col = 'magenta')
hist(to.e2, main = "B) FPKM between 10 and 100", xlab = "FPKM", col = 'skyblue2')
hist(to.e3, main = "C) FPKM between 100 and 1000", xlab = "FPKM", col = 'pink')
hist(to.e4, main = "D) FPKM greater than 1000", xlab = "FPKM", col = 'lavender')
dev.off()
