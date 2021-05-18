############### Figure 1E


library("ggbio")
library(GenomicRanges)
genome = read.csv("hg38_chro_length.txt",header=F,sep="\t")

DMSO = read.table("K562_DMSO_merged_sort_collapse_average_1MB.bed",header=T)
VP16 = read.table("K562_VP16_merged_sort_collapse_average_1MB.bed",header=T)
pBQ = read.table("K562_pBQ_merged_sort_collapse_average_1MB.bed",header=T)

DMSO = DMSO[order(-DMSO[,5]),]
VP16 = VP16[order(-VP16[,5]),]
pBQ = pBQ[order(-pBQ[,5]),]

DMSO[1,5] = 1500
VP16[1,5] = 900
pBQ[1,5] = 900

DMSO = DMSO[order(DMSO[,1]),]
VP16 = VP16[order(VP16[,1]),]
pBQ = pBQ[order(pBQ[,1]),]


Mesculenta <- GRanges(seqnames = c(1:24), ranges = IRanges(start=1, width = as.numeric(genome[1:24,2])))

seqlevels(Mesculenta)<-as.character(c(1:22,"X","Y"))
seqlengths(Mesculenta)<-as.numeric(genome[1:24,2])



DMSO_bed <- GRanges(seqnames= DMSO[,1], ranges = IRanges(start = DMSO[,2]+1,width = DMSO[,3] - DMSO[,2]))
VP16_bed <- GRanges(seqnames= VP16[,1], ranges = IRanges(start = VP16[,2]+1,width = VP16[,3] - VP16[,2]))
pBQ_bed <- GRanges(seqnames= pBQ[,1], ranges = IRanges(start = pBQ[,2]+1,width = pBQ[,3] - pBQ[,2]))

DMSO_bed$score = DMSO[,5]
VP16_bed$score = VP16[,5]
pBQ_bed$score = pBQ[,5]

#png(width=600,height=600,units='px','circos.png')

#pdf("peaks_circos.pdf")
#tiff("Fig_1C_CEM_peaks_circos.tif", units= "cm", width = 30, height = 30, res = 300, family = "Arial")
#par(cex = 0.5, lwd = 1, mar = c(2,2,2,1) + 0.1, mgp = c(2.6, 0.8, 0))

pdf("Fig_1C_CEM_peaks_circos.pdf")#,height = 6/2.54, width = 6/2.54)#,pointsize = 10)
#par(cex =1, lwd = 1, mar = c(4,5, 1, 1) + 0.1, mgp = c(2.6, 0.8, 0))

ggplot() + layout_circle(Mesculenta, geom = "ideo", fill = "orange", radius = 39, trackWidth = 1) + layout_circle(Mesculenta, geom = "text", aes(label = seqnames), vjust = 0, radius = 40, trackWidth = 10,size = 15) + layout_circle(DMSO_bed, geom= 'bar', color = "brown3", aes(x=start,y=score,fill='a'),radius = 9, trackWidth = 10,size =0.1) + layout_circle(VP16_bed, geom= 'bar', color = "darkolivegreen4", aes(x=start,y=score,fill='b'),radius = 19, trackWidth = 10,size =0.1) + layout_circle(pBQ_bed, geom= 'bar', color = "royalblue3", aes(x=start,y=score,fill='c'),radius = 29, trackWidth = 10,size =0.1)

dev.off()
