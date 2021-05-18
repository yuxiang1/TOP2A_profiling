

### Figure_1D
CEM_A <- read.table("CEM_DMSO_merged_sort_collapse_average_1MB.bed", sep = "\t", header = T, as.is = TRUE)
CEM_B <- read.table("CEM_VP16_merged_sort_collapse_average_1MB.bed", sep = "\t", header = T, as.is = TRUE) #
CEM_C <- read.table("CEM_pBQ_merged_sort_collapse_average_1MB.bed", sep = "\t", header = T, as.is = TRUE) #

chr_length <- read.table("hg38_chro_length.txt", sep = "\t", header = F, as.is = TRUE)
chr_length[,2] = chr_length[,2]/1000000

CEM_A_chr_density = matrix(0,23,2)
CEM_A_chr_density[,1] = c(1:22,"X")

CEM_B_chr_density = matrix(0,23,2)
CEM_B_chr_density[,1] = c(1:22,"X")

CEM_C_chr_density = matrix(0,23,2)
CEM_C_chr_density[,1] = c(1:22,"X")

for(i in 1:23){
    CEM_A_chr_density[i,2] = mean(subset(CEM_A,chrom==i)$clone)
    CEM_B_chr_density[i,2] = mean(subset(CEM_B,chrom==i)$clone)
    CEM_C_chr_density[i,2] = mean(subset(CEM_C,chrom==i)$clone)
}


cor.test(chr_length[1:23,2],as.numeric(CEM_A_chr_density[,2])) # 0.03129
cor.test(chr_length[1:23,2],as.numeric(CEM_B_chr_density[,2])) # 0.009806
cor.test(chr_length[1:23,2],as.numeric(CEM_C_chr_density[,2])) # 0.01405


##### density
pdf("Figure_1D.pdf",height = 6/2.54, width = 7/2.54,pointsize = 10)
#par(mar = c(5,4,1,1) + 0.1, mgp = c(2.6, 0.8, 0))
#tiff("Fig_1B_CEM_DMSO_CCR_along_chr.tif", units= "cm", width =5, height = 5, res = 300, family = "Arial")
par(cex =1, lwd = 1, mar = c(4,5, 1, 1) + 0.1, mgp = c(2.6, 0.8, 0))
plot(chr_length[1:23,2],as.numeric(CEM_A_chr_density[,2]),type = "p",xlab ="Chromosome Length (Mbp)", ylab = "CCR Signal Density\n(HP10M/Mb)",pch = 16,col = "brown3",las=1,xlim = c(40,260),ylim = c(400,1050))

average = mean(as.numeric(CEM_A_chr_density[,2]))
pos <- rep(1, (nrow(chr_length)-1))
pos[c(10,12,13,16,17,19,21,22)] <- c(2,3,3,3,2,4,3,3)
text(chr_length[1:23,2],as.numeric(CEM_A_chr_density[,2]), labels = CEM_A_chr_density[,1], pos = pos)#, cex = 0.6)
# Add line
abline(h = average, lwd = 1, lty = 2)
cor <- cor(chr_length[1:23,2],as.numeric(CEM_A_chr_density[,2]))
cor <- as.numeric(sprintf("%.3f", cor))
legend("bottomright",legend = bquote(italic(r) == .(format(cor,digits=3))), pch = 16,col=c("white"),bty = "n",plot=T)
legend("topright",legend = c("DMSO"), pch = 16,col=c("white"),bty = "n",plot=T)
dev.off()
