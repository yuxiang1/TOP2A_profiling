
Sample1 <- read.csv("CEM_DMSO_CCR_anno_type.txt",sep="\t",header=F,check.names=F,row.names=1)
Sample2 <- read.csv("CEM_VP16_CCR_anno_type.txt",sep="\t",header=F,check.names=F,row.names=1)
Sample3 <- read.csv("CEM_pBQ_CCR_anno_type.txt",sep="\t",header=F,check.names=F,row.names=1)
Sample4 <- read.csv("random-CCR_anno_type.txt",sep="\t",header=F,check.names=F,row.names=1)
type = c("promoter-TSS ","5' UTR ","exon ","intron ","3' UTR ","TTS ","non-coding ","Intergenic")

Sample1[,2] = Sample1[,1]
Sample2[,2] = Sample2[,1]
Sample3[,2] = Sample3[,1]
Sample4[,2] = Sample4[,1]
Sample1 = Sample1[type,]
Sample2 = Sample2[type,]
Sample3 = Sample3[type,]
Sample4 = Sample4[type,]
## data from annotated CCRs overlapped with lincRNAs
lincRNAs  =c(35703,57312,24738,6188)
Sample1["Intergenic",1] = Sample1["Intergenic",1] - (lincRNAs[1] -Sample1["non-coding ",1])
Sample2["Intergenic",1] = Sample2["Intergenic",1] - (lincRNAs[2] -Sample2["non-coding ",1])
Sample3["Intergenic",1] = Sample3["Intergenic",1] - (lincRNAs[3] -Sample3["non-coding ",1])
Sample4["Intergenic",1] = Sample4["Intergenic",1] - (lincRNAs[4] -Sample4["non-coding ",1])


Sample1["non-coding ",1] = lincRNAs[1]
Sample2["non-coding ",1] = lincRNAs[2]
Sample3["non-coding ",1] = lincRNAs[3]
Sample4["non-coding ",1] = lincRNAs[4]

type = rownames(Sample4)
Sample1[type,1] =log2(Sample1[type,1]/sum(Sample1[type,1]) /(Sample4[type,1]/sum(Sample4[type,1])))
Sample2[type,1] =log2(Sample2[type,1]/sum(Sample2[type,1]) /(Sample4[type,1]/sum(Sample4[type,1])))
Sample3[type,1] =log2(Sample3[type,1]/sum(Sample3[type,1]) /(Sample4[type,1]/sum(Sample4[type,1])))

enrichment = cbind(Sample1[type,1],Sample2[type,1],Sample3[type,1])
#Sample1 = Sample1[-8,]
#enrichment = enrichment[-8,]
enrichment = t(enrichment)
type = c("Promoter","5' UTR","Exon","Intron","3' UTR","TTS ","LincRNAs","Intergenic")
#tiff("CEM_peak_enrichment.tif", units= "cm", width = 5, height = 5, res = 300, family = "Arial")
#par(cex = 0.5, lwd = 0.5, mar = c(6,4,2,2) + 0.1, mgp = c(2.6, 0.8, 0))
pdf("Fig_1B.pdf",height = 7/2.54, width = 7/2.54,pointsize = 10)
par(cex =1, lwd = 1, mar = c(6,4,1,1) + 0.1, mgp = c(2.6, 0.8, 0))

barplot(enrichment,col = c("brown3","darkolivegreen4","royalblue3"), name = as.character(type), beside = TRUE, las = 2, ylab = expression(paste(Log[2],"(Enrichment)",sep="")),ylim = c(-1,1))
# Add legend
par(xpd = TRUE)
legend("topleft", legend = c("DMSO","VP16","pBQ"), fill = c("brown3","darkolivegreen4","royalblue3"),bty = "n")
par(xpd = FALSE)
dev.off()
