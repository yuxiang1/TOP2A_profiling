
############### Figure 4B

CEM_A = read.table("CEM_1A_3A_merged_CCR_CEL.txt",sep= "\t",row.names=1)
CEM_B = read.table("CEM_1B_3B_merged_CCR_CEL.txt",sep= "\t",row.names=1)
CEM_C = read.table("CEM_1C_3C_merged_CCR_CEL.txt",sep= "\t",row.names=1)

CCR_A = subset(CEM_A,CEM_A[,1]>0)
noCCR_A = subset(CEM_A,CEM_A[,1]==0)

CCR_B = subset(CEM_B,CEM_B[,1]>0)
noCCR_B = subset(CEM_B,CEM_B[,1]==0)


CCR_C = subset(CEM_C,CEM_C[,1]>0)
noCCR_C = subset(CEM_C,CEM_C[,1]==0)


type <- c(rep("DMSO_Y",nrow(CCR_A)),rep("DMSO_N",nrow(noCCR_A)),rep("VP16_Y",nrow(CCR_B)),rep("VP16_N",nrow(noCCR_B)),rep("pBQ_Y",nrow(CCR_C)),rep("pBQ_N",nrow(noCCR_C)))

sample <-c(CCR_A[,2],noCCR_A[,2],CCR_B[,2],noCCR_B[,2],CCR_C[,2],noCCR_C[,2])


data <- data.frame(sample,type)
names(data)=c("sample","type")
data$type <- factor(data$type, levels=unique(data$type))

pdf("Figure_4B.pdf",height = 8/2.54, width = 8/2.54,pointsize = 10)
par(cex =1, lwd = 1, mar = c(4,4,1,1) + 0.1, mgp = c(2.6, 0.8, 0),las=1)


plot(1,type = "n",xlab ="Gene expression (RMA)", ylab = "Probability Density",xlim = c(0,15), ylim = c(0,0.3))
ADJ = 1
lines(density((noCCR_A[,2]),adjust=ADJ),col = "gray20", lwd = 1)
lines(density((noCCR_B[,2]),adjust=ADJ),col = "gray40", lwd = 1)
lines(density((noCCR_C[,2]),adjust=ADJ),col = "gray60", lwd = 1)

lines(density((CCR_A[,2]),adjust=ADJ),col = "brown3", lwd = 1)
lines(density((CCR_B[,2]),adjust=ADJ),col = "darkolivegreen4", lwd = 1)
lines(density((CCR_C[,2]),adjust=ADJ),col = "royalblue3", lwd = 1)

legend("topright",legend = c("CCR-, DMSO","CCR-, VP16","CCR-, pBQ","CCR+, DMSO","CCR+, VP16","CCR+, pBQ"), lwd = 1,col=c("gray20","gray40","gray60","brown3","darkolivegreen4","royalblue3"),bty="n")
dev.off()
