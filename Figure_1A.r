### Read data
CEM_A <- read.table("CEM_1A_3A_merged_sort_collapse_average.bed", sep = "\t", header = F, as.is = TRUE)
CEM_B <- read.table("CEM_1B_3B_merged_sort_collapse_average.bed", sep = "\t", header = F, as.is = TRUE) #
CEM_C <- read.table("CEM_1C_3C_merged_sort_collapse_average.bed", sep = "\t", header = F, as.is = TRUE) #

##### density

CEM_A$length = (CEM_A[,3] - CEM_A[,2])/100
CEM_B$length = (CEM_B[,3] - CEM_B[,2])/100
CEM_C$length = (CEM_C[,3] - CEM_C[,2])/100

tiff("Figure_1A.tif", units= "cm", width =5, height = 5, res = 300, family = "Arial")
par(cex = 0.6, lwd = 0.8, mar = c(4, 4, 1, 1) + 0.1, mgp = c(2.6, 0.8, 0))

plot(1,type = "n",xlab ="CCR length (100 bp)", ylab = "Density",xlim = c(0,10), ylim = c(0,1))

lines(density((CEM_A$length),adjust=1),col = "red", lwd = 0.8)
lines(density((CEM_B$length),adjust=1),col = "green", lwd = 0.8)
lines(density((CEM_C$length),adjust=1),col = "blue", lwd = 0.8)

legend("topright",legend = c("DMSO","VP16","pBQ"), cex = 0.6, lwd = 0.8,col=c("red","green","blue"))
dev.off()



