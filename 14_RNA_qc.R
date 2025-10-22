setwd(save_wd)
#####loading RNA data####
ready_RH23_24h<-read.csv("DE_H23_2Gy_vs_0Gy.csv")
ready_RH23_24h<- ready_RH23_24h[which(ready_RH23_24h$gene_biotype=="protein_coding"),]
ready_RHCC44_24h<-read.csv("DE_HCC44_2Gy_vs_0Gy.csv")
ready_RHCC44_24h<- ready_RHCC44_24h[which(ready_RHCC44_24h$gene_biotype=="protein_coding"),]
ready_RRAD_24h<-read.csv("DE_HCC44_vs_H23_2Gy.csv")
ready_RRAD_24h<- ready_RRAD_24h[which(ready_RRAD_24h$gene_biotype=="protein_coding"),]
ready_RCon_24h<-read.csv("DE_HCC44_vs_H23_0Gy.csv")
ready_RCon_24h<- ready_RCon_24h[which(ready_RCon_24h$gene_biotype=="protein_coding"),]

RH23_24h<- RH23_24h[which(RH23_24h$gene_biotype=="protein_coding"),]
RHCC44_24h<- RHCC44_24h[which(RHCC44_24h$gene_biotype=="protein_coding"),]
#####Transcriptome - Boxplot Intensity####
box_HCC44.H23_24h <- cbind(RH23_24h[,3:9],RHCC44_24h[,3:10])
colnames(box_HCC44.H23_24h)<-  c("","H23 24h control","","H23 24h 2Gy","","","","HCC44 24h control","","","","","","","HCC44 24h 2Gy")
#####Transcriptome - Intensities and counts####
jpeg(filename = "Transcriptome_counts.jpeg",  width=20.5, height=14, unit="cm",quality=120,res=90)
layout(matrix(c(1,1,2, 2,2,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar = c(10, 10, 4, 4),xpd=TRUE)
vioplot(log10(2^box_HCC44.H23_24h[,c(8:15,1:7)]),pchMed=15,lwd=1.5,xaxt = "n",yaxt="n",cex.axis=4,ylog=FALSE,ann=FALSE,axes=FALSE,cex.lab=3.5,ylim=c(-1,6),boxwex=0.5, col=c(rep("#3288BD",4),rep("#D53E4F",4),rep("#66C2A5",3),rep("#FDAE61",4)))
box(which = "plot",lwd=2)
axis(1,lwd=2,at=c(1,5,9,12),labels = FALSE,lwd.ticks =2,cex=3, cex.axis=2.5,tcl = -1)
axis(2,lwd=2,at=c(0,5),lwd.ticks = 2,cex=3, cex.axis=2.5,tcl = -1.5,las=2,mgp=c(12,2,0))
text(x=1.1,y=-3.25,labels="NCI-H23 control",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=5.1,y=-3.25,labels="NCI-H23 2Gy",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=9.2,y=-3.25,labels="HCC44 control",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=12.2,y=-3.25,labels="HCC44 2Gy",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
lines(c(0,16),c(0,0),lty=3,col="black",lwd=2)
title(ylab=(log[10]~"(Count)"), cex.lab=2.5,cex.main=3,line=6,adj=c(0.7))
title(main="Transcriptome - HCC44 & NCI-H23", cex.lab=3,cex.main=2.5,line=1.5,adj=c(0.5))
minor.tick(nx =1,ny=1.00025,tick.ratio =1, x.args = list( lwd = 2,tcl = -0.75), y.args = list( lwd = 2,tcl = -0.75)) 
dev.off()
#####Transcriptome - Identified RNAs - counting####
n_Rcount<-NULL
n_Rcount["H23_24h_c"]<-data.frame(count(RH23_24h[complete.cases(RH23_24h$mean_c),]))
n_Rcount["H23_24h_2Gy"]<-data.frame(count(RH23_24h[complete.cases(RH23_24h$mean_2Gy),]))
n_Rcount["HCC44_24h_c"]<-data.frame(count(RHCC44_24h[complete.cases(RHCC44_24h$mean_c),]))
n_Rcount["HCC44_24h_2Gy"]<-data.frame(count(RHCC44_24h[complete.cases(RHCC44_24h$mean_2Gy),]))
#####Transcriptome - Identified RNAs - plotting####
jpeg(filename = "Transcriptome_identifiedRNAs.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(16, 14, 6, 4),xpd=TRUE)
barplot(as.numeric(n_Rcount[c(3,4,1,2)]), col=c("#3288BD","#D53E4F","#ABDDA4", "#FDAE61"),xlim=c(0,5),space=0.5,width=0.75,axes=FALSE,ylim=c(0,20000))
minor.tick(nx =1,ny=2,tick.ratio =1, x.args = list( lwd = 3), y.args = list( lwd = 3,tcl = -0.75)) 
axis(2,lwd=3,at = c(0,5000,10000,15000,20000),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,las=2,mgp=c(6,2,0))
axis(1,lwd=3,at = c(0,5.5),labels=c("",""),lwd.ticks = 3,cex=4, cex.axis=4,tcl = 0)
title( ylab="detected RNAs", cex.lab=3,cex.main=3,line=10,adj=1)
text(x = 0.8,y = -1500,adj=c(1,0),labels = "HCC44 24h c" ,srt = 45,cex = 2.5,col="black")
text(x = 2.0,y = -1500,adj=c(1,0),labels ="HCC44 24h 2Gy" ,srt = 45,cex = 2.5,col="black")
text(x = 3.15,y = -1500,adj=c(1,0),labels ="H23 24h c"  ,srt = 45,cex = 2.5,col="black")
text(x = 4.3,y =-1500,adj=c(1,0),labels ="H23 24h 2Gy"   ,srt = 45,cex = 2.5,col="black")
text(x = c(0.5),y = 15500,adj=c(0,0),labels =n_Rcount[3],srt = 45,cex = 2.5,col="black")
text(x = c(1.8),y = 15500,adj=c(0,0),labels =n_Rcount[4],srt = 45,cex = 2.5,col="black")
text(x = c(3.0),y = 15500,adj=c(0,0),labels =n_Rcount[1],srt = 45,cex = 2.5,col="black")
text(x = c(4.1),y = 15500,adj=c(0,0),labels =n_Rcount[2],srt = 45,cex = 2.5,col="black")
text(x = c(4.6),y = 14400,adj=c(0,0),labels ="13.751",srt = 0,cex = 2.25,col="black")
lines(c(0,6.5),c(13751,13751), lty=3,lwd=3,col="black")
dev.off()
#####Transcriptome - calculating PCA#####
intens_RNA_all_24h <- data.frame(log2(RH23_24h[,3:9]),log2(RHCC44_24h[,3:10]))
intens_RNA_all_24h<-data.frame(intens_RNA_all_24h[complete.cases(intens_RNA_all_24h),])
intens_RNA_all_24h_filt<- t(data.frame(intens_RNA_all_24h))
pc_RNA_all_24h_filt <- prcomp(intens_RNA_all_24h_filt,
                              center = TRUE,
                              scale = TRUE)
#####Transcriptome - plotting PCA####
jpeg(filename = "Transcriptome_PCA.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar = c(10, 10, 6, 4),xpd=TRUE)
leg.txt <- c("HCC44 c","HCC44 2Gy","NCI-H23 c", "NCI-H23 2Gy")
s_pc_RNA_all_24h <- summary(pc_RNA_all_24h_filt)
plot(as.data.frame(s_pc_RNA_all_24h["x"])[,1], as.data.frame(s_pc_RNA_all_24h["x"])[,2],pch=21, xlab=paste("PCA 1 (", round(s_pc_RNA_all_24h$importance[2]*100, 1), "%)", sep = ""), 
     ylab=paste("PCA 2 (", round(s_pc_RNA_all_24h$importance[5]*100, 1), "%)", sep = ""),bg=c(rep( "#66C2A5",3),rep("#FDAE61",4),rep("#3288BD",4), rep("#D53E4F",4)),  col="black",
     cex=5,cex.main=2,cex.lab=4,xlim=c(-125,125),ylim=c(-125,125), cex.axis=4,ann=FALSE,axes=FALSE,mgp=c(7,0,0),lwd=2,main="RNA")
axis(1,lwd=3,at = c(-100,0,100),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(6,3,0))
axis(2,lwd=3,at = c(-100,0,100),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(6,2,0),las=2)
lines(c(0,0),c(-130,130), lty=3,lwd=2)
lines(c(-105,105),c(0,0), lty=3,lwd=2)
minor.tick(nx =2/1,ny=2/1,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
title(main="Transcriptome (24h after radiation)",cex.lab=3,cex.main=3,line=2,adj=c(0.7))
title(xlab=paste("PCA 1 (", round(s_pc_RNA_all_24h$importance[2]*100, 1), "%)", sep = ""),cex.lab=3.5,cex.main=3,line=6,adj=c(0.5))
title(ylab=paste("PCA 2 (", round(s_pc_RNA_all_24h$importance[5]*100, 1), "%)", sep = ""),cex.lab=3.5,cex.main=3,line=7,adj=c(0.5))
box(lwd=3,lty="solid")
legend(x=c(-100,100),y=c(-75,-125),ncol=2, box.col="white" ,bty = "n", bg = NA ,leg.txt,pch=21,cex = 2.25,pt.cex = 3,pt.bg = c("#3288BD","#D53E4F","#66C2A5","#FDAE61"),  col="black",box.lwd =2,pt.lwd = 2,y.intersp=1,x.intersp=0.75)
dev.off()

#####Transcriptome - calculation COV####
COV_RHCC44_24h_2Gy<- RHCC44_24h[,c(1,2,15,16)]
COV_RHCC44_24h_2Gy[["%CV"]] <-(COV_RHCC44_24h_2Gy[,4]/COV_RHCC44_24h_2Gy[,3])*100
COV_RHCC44_24h_2Gy<-  COV_RHCC44_24h_2Gy[complete.cases(COV_RHCC44_24h_2Gy[,4]),]
COV_RHCC44_24h_2Gy<- COV_RHCC44_24h_2Gy[which(COV_RHCC44_24h_2Gy$mean_2Gy>COV_RHCC44_24h_2Gy$SD_2Gy),]
dens_RHCC44_24h_2Gy_cov<-density(as.numeric(COV_RHCC44_24h_2Gy[,5]),na.rm=TRUE)

COV_RHCC44_24h_c<- RHCC44_24h[,c(1,2,13,14)]
COV_RHCC44_24h_c[["%CV"]] <-(COV_RHCC44_24h_c[,4]/COV_RHCC44_24h_c[,3])*100
COV_RHCC44_24h_c<-  COV_RHCC44_24h_c[complete.cases(COV_RHCC44_24h_c[,4]),]
COV_RHCC44_24h_c<- COV_RHCC44_24h_c[which(COV_RHCC44_24h_c$mean_c>COV_RHCC44_24h_c$SD_c),]
dens_RHCC44_24h_c_cov<-density(COV_RHCC44_24h_c[,5],na.rm=TRUE)

COV_RH23_24h_2Gy<- RH23_24h[,c(1,2,14,15)]
COV_RH23_24h_2Gy[["%CV"]] <-(COV_RH23_24h_2Gy[,4]/COV_RH23_24h_2Gy[,3])*100
COV_RH23_24h_2Gy<-  COV_RH23_24h_2Gy[complete.cases(COV_RH23_24h_2Gy[,4]),]
COV_RH23_24h_2Gy<- COV_RH23_24h_2Gy[which(COV_RH23_24h_2Gy$mean_2Gy>COV_RH23_24h_2Gy$SD_2Gy),]
dens_RH23_24h_2Gy_cov<-density(COV_RH23_24h_2Gy[,5],na.rm=TRUE)

COV_RH23_24h_c<- RH23_24h[,c(1,2,12,13)]
COV_RH23_24h_c[["%CV"]] <-(COV_RH23_24h_c[,4]/COV_RH23_24h_c[,3])*100
COV_RH23_24h_c<-  COV_RH23_24h_c[complete.cases(COV_RH23_24h_c[,4]),]
COV_RH23_24h_c<- COV_RH23_24h_c[which(COV_RH23_24h_c$mean_c>COV_RH23_24h_c$SD_c),]
dens_RH23_24h_c_cov<-density(COV_RH23_24h_c[,5],na.rm=TRUE)
#####Transcriptome - plotting CoV####
jpeg(filename = "Transcriptome_CoV.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1),ncol=3,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt<-c("NCI-H23 c","NCI-H23 2Gy","HCC44 c","HCC44 2Gy")
plot(dens_RHCC44_24h_2Gy_cov$x,length(COV_RHCC44_24h_2Gy[,3])*dens_RHCC44_24h_2Gy_cov$y,ylim=c(0,4000),xlim=c(-10,120),type="l",lwd=2,
     col= "#D53E4F",axes=FALSE,cex.lab=3,mgp=c(5,0,1),ann=FALSE)
lines(dens_RHCC44_24h_c_cov$x,length(COV_RHCC44_24h_c[,3])*dens_RHCC44_24h_c_cov$y, col="#3288BD",lwd=2)
lines(dens_RH23_24h_2Gy_cov$x,length(COV_RH23_24h_2Gy[,3])*dens_RH23_24h_2Gy_cov$y, col="#66C2A5",lwd=2)
lines(dens_RH23_24h_c_cov$x,length(COV_RH23_24h_c[,3])*dens_RH23_24h_c_cov$y, col="#FDAE61",lwd=2)
axis(1,lwd=3,at=c(0,50,100),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,3,0))
axis(2,lwd=3,at = c(0,2000,4000),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,2,0),las=2)
lines(c(15,15),c(-10,4000), lty=3,lwd=3)
box(type="solid",lwd=3)
legend(x=c(50,120),y=c(3400,4000),leg.txt,bty = "n", bg = NA ,lty=1,box.col="white",lwd=2,adj=c(0,0.25),  col=c("#66C2A5","#FDAE61","#3288BD", "#D53E4F"),inset=0.02,cex=2.75,box.lwd = 2)
title(ylab="Density", cex.lab=3,cex.main=3,line=9,adj=c(0.5))
title(xlab="Coefficient of Variation[%]", cex.lab=3,cex.main=3,line=6,)
title(main="Transcriptome", cex.lab=3,cex.main=3,line=2,adj=c(0.5))
minor.tick(nx =2/1,ny=2/1,tick.ratio =0.75, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
dev.off()
setwd("C:/Users/annar/Desktop/Recherche/Ergebnisse/RNA-Seq/p1472/quality control")
RPKM<- data.frame(read.csv("RPKM_Bohnenberger.csv"))
TPM<- data.frame(read.csv("TPM_Bohnenberger.csv"))

#####Transcriptome - extracting RPKM#####
RPKM_H23_c <- cbind(RPKM[,c(1,2,8,13,21,22)])
RPKM_H23_c["Mean"] <- log2(rowMeans(RPKM_H23_c[,4:6],na.rm=TRUE))
RPKM_H23_c[sapply(RPKM_H23_c, is.infinite)] <- NA
RNA_type_RPKM_H23_c<-data.frame((count(RPKM_H23_c[complete.cases(RPKM_H23_c[,7]),], gene_biotype)))
RNA_type_RPKM_H23_c<-RNA_type_RPKM_H23_c[order(RNA_type_RPKM_H23_c[,2],decreasing=TRUE),]

RPKM_HCC44_2Gy <- cbind(RPKM[,c(1,2,8,10:12,20)])
RPKM_HCC44_2Gy["Mean"] <- log2(rowMeans(RPKM_HCC44_2Gy[,4:7],na.rm=TRUE))
RPKM_HCC44_2Gy[sapply(RPKM_HCC44_2Gy, is.infinite)] <- NA
RNA_type_RPKM_HCC44_2Gy<-data.frame((count(RPKM_HCC44_2Gy[complete.cases(RPKM_HCC44_2Gy[,8]),], gene_biotype)))
RNA_type_RPKM_HCC44_2Gy<-RNA_type_RPKM_HCC44_2Gy[order(RNA_type_RPKM_HCC44_2Gy[,2],decreasing=TRUE),]

RPKM_HCC44_c <- cbind(RPKM[,c(1,2,8,9,17:19)])
RPKM_HCC44_c["Mean"] <- log2(rowMeans(RPKM_HCC44_c[,4:7],na.rm=TRUE))
RPKM_HCC44_c[sapply(RPKM_HCC44_c, is.infinite)] <- NA
RNA_type_RPKM_HCC44_c<-data.frame((count(RPKM_HCC44_c[complete.cases(RPKM_HCC44_c[,8]),], gene_biotype)))
RNA_type_RPKM_HCC44_c<-RNA_type_RPKM_HCC44_c[order(RNA_type_RPKM_HCC44_c[,2],decreasing=TRUE),]

RPKM_H23_2Gy <- cbind(RPKM[,c(1,2,8,23,14:16)])
RPKM_H23_2Gy["Mean"] <- log2(rowMeans(RPKM_H23_2Gy[,4:7],na.rm=TRUE))
RPKM_H23_2Gy[sapply(RPKM_H23_2Gy, is.infinite)] <- NA
RNA_type_RPKM_H23_2Gy<-data.frame((count(RPKM_H23_2Gy[complete.cases(RPKM_H23_2Gy[,8]),], gene_biotype)))
RNA_type_RPKM_H23_2Gy<-RNA_type_RPKM_H23_2Gy[order(RNA_type_RPKM_H23_2Gy[,2],decreasing=TRUE),]
setwd(save_wd)
#####Transcriptome - plotting RPKM violine plot####
jpeg(filename = "Transcriptome_RPKM.jpeg",  width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(3,3,4,4,5, 3,3,4,4,5,1,1,2,2,5,  1,1,2,2,5),ncol=5,nrow = 4,byrow = TRUE))
par(mar =c(2, 10, 6, 2),xpd=TRUE)
vioplot(as.numeric(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[1,1],7]),as.numeric(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[2,1],7]),as.numeric(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[3,1],7]),
        as.numeric(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[4,1],7]),as.numeric(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[5,1],7]),rbind(as.matrix((RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[6,1],7])),
                                                                                                                                                                                  as.matrix((RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[7,1],7])),as.matrix(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[8,1],7]),as.matrix(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[9,1],7]),
                                                                                                                                                                                  as.matrix(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[10,1],7]),as.matrix(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[11,1],7]),as.matrix(RPKM_H23_c[RPKM_H23_c$gene_biotype %like% RNA_type_RPKM_H23_c[12,1],7])),
        col=brewer.pal(6, "Spectral"),ylim=c(-5,20),xlim=c(0.5,6.5),pchMed=15,lwd=3,yaxt = "n",xaxt = "n",cex.axis=2,ann=FALSE,cex.lab=3,boxwex=0.5,main="H23 control",cex.main=3)
axis(1,lwd=3, labels = FALSE,lwd.ticks = 4,cex=3, cex.axis=3,tcl = 0, mgp=c(0,0,0))
axis(2,lwd=3,at = c(0,10,20),lwd.ticks = 4,cex=3, cex.axis=2.5,tcl = -1.5,las=2, mgp=c(0,2,0))
minor.tick(nx =1,ny=2/1,tick.ratio =0.75, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
box(which = "plot",lwd=3)
title(main="H23 control",cex.main=3, cex.lab=3,line=2)
title( ylab=bquote(log[2]~"RPKM"),cex.main=3, cex.lab=3,line=5)

vioplot(as.numeric(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[1,1],8]),as.numeric(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[2,1],8]),as.numeric(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[3,1],8]),
        as.numeric(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[4,1],8]),as.numeric(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[5,1],8]),rbind(as.matrix((RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[6,1],8])),
                                                                                                                                                                                              as.matrix((RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[7,1],8])),as.matrix(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[8,1],8]),as.matrix(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[9,1],8]),
                                                                                                                                                                                              as.matrix(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[10,1],8]),as.matrix(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[11,1],8]),as.matrix(RPKM_H23_2Gy[RPKM_H23_2Gy$gene_biotype %like% RNA_type_RPKM_H23_2Gy[12,1],8])),
        col=brewer.pal(6, "Spectral"),ylim=c(-5,20),xlim=c(0.5,6.5),pchMed=15,lwd=3,yaxt = "n",xaxt = "n",cex.axis=2,ann=FALSE,cex.lab=3,boxwex=0.5,main="H23 control",cex.main=3)
axis(1,lwd=3, labels = FALSE,lwd.ticks = 3,cex=3, cex.axis=3,tcl = 0, mgp=c(0,0,0))
axis(2,lwd=3,at = c(0,10,20),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,las=2, mgp=c(0,2,0))
minor.tick(nx =1,ny=2/1,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
box(which = "plot",lwd=3)
title(main="H23 2Gy",cex.main=3, cex.lab=3,line=2)
title( ylab=bquote(log[2]~"RPKM"),cex.main=3, cex.lab=3,line=5)

vioplot(as.numeric(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[1,1],8]),as.numeric(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[2,1],8]),as.numeric(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[3,1],8]),
        as.numeric(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[4,1],8]),as.numeric(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[5,1],8]),rbind(as.matrix((RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[6,1],8])),
                                                                                                                                                                                              as.matrix((RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[7,1],8])),as.matrix(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[8,1],8]),as.matrix(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[9,1],8]),
                                                                                                                                                                                              as.matrix(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[10,1],8]),as.matrix(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[11,1],8]),as.matrix(RPKM_HCC44_c[RPKM_HCC44_c$gene_biotype %like% RNA_type_RPKM_HCC44_c[12,1],8])),
        col=brewer.pal(6, "Spectral"),ylim=c(-5,20),xlim=c(0.5,6.5),pchMed=15,lwd=3,yaxt = "n",xaxt = "n",cex.axis=2,ann=FALSE,cex.lab=3,boxwex=0.5,main="H23 control",cex.main=3)
axis(1,lwd=3, labels = FALSE,lwd.ticks = 3,cex=3, cex.axis=3,tcl = 0, mgp=c(0,0,0))
axis(2,lwd=3,at = c(0,10,20),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,las=2, mgp=c(0,2,0))
minor.tick(nx =1,ny=2/1,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
box(which = "plot",lwd=3)
title(main="HCC44 control",cex.main=3, cex.lab=2.5,line=2)
title( ylab=bquote(log[2]~"RPKM"),cex.main=3, cex.lab=3,line=5)

vioplot(as.numeric(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[1,1],8]),as.numeric(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[2,1],8]),as.numeric(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[3,1],8]),
        as.numeric(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[4,1],8]),as.numeric(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[5,1],8]),rbind(as.matrix((RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[6,1],8])),
                                                                                                                                                                                                          as.matrix((RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[7,1],8])),as.matrix(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[8,1],8]),as.matrix(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[9,1],8]),
                                                                                                                                                                                                          as.matrix(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[10,1],8]),as.matrix(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[11,1],8]),as.matrix(RPKM_HCC44_2Gy[RPKM_HCC44_2Gy$gene_biotype %like% RNA_type_RPKM_HCC44_2Gy[12,1],8])),
        col=brewer.pal(6, "Spectral"),ylim=c(-5,20),xlim=c(0.5,6.5),pchMed=15,lwd=3,yaxt = "n",xaxt = "n",cex.axis=2,ann=FALSE,cex.lab=3,boxwex=0.5,main="H23 control",cex.main=3)
axis(1,lwd=3, labels = FALSE,lwd.ticks = 3,cex=3, cex.axis=3,tcl = 0, mgp=c(0,0,0))
axis(2,lwd=3,at = c(0,10,20),lwd.ticks =3 ,cex=3, cex.axis=2.5,tcl = -1.5,las=2, mgp=c(0,2,0))
minor.tick(nx =1,ny=2/1,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
box(which = "plot",lwd=3)
title(main="HCC44 2Gy",cex.main=3, cex.lab=2.5,line=2)
title( ylab=bquote(log[2]~"RPKM"),cex.main=3, cex.lab=3,line=5)
par(mar =c(2, 4, 2, 4),xpd=TRUE)

par(mar =c(2, 4, 6, 4),xpd=TRUE)
lgd=as.vector(c("protein \n coding","lncRNA","processed \n pseudogene", "TEC" , "unprocessed \n pseudogene","rest"))
plot(1:10,1:10,ann=FALSE,axes=FALSE,col="white")
legend(x=c(-0.5,8),y=c(2,10.25),ncol = 1,adj=c(0,0.25),legend=lgd, pch=16, col=brewer.pal(6, "Spectral"),pt.cex=3, cex = 2.25,box.lwd=2,y.intersp = 1.3)
dev.off()
#####Transcriptome - Venn diagramm####
venn.diagram(list(RH23_24h[complete.cases(RH23_24h[,c(13)]),][,3],RH23_24h[complete.cases(RH23_24h[,c(15)]),][,3],RHCC44_24h[complete.cases(RHCC44_24h[,c(14)]),][,3],RHCC44_24h[complete.cases(RHCC44_24h[,c(16)]),][,3]),
             category.names = c("H23 24h c","H23 24h 2Gy","HCC44 24h c","HCC44 24h 2Gy"),
             filename = 'Transcriptome_VennDiagramm.jpeg',
             print.mode = "raw",na="none",
             fill = c("#66C2A5","#FDAE61","#3288BD","#D53E4F"),
             ext.percent=c(2,2,0), main="RNA",main.cex=6,main.lwd=4,fontfamily= 1,fontface=5,cat.fontfamily=1,
             cex=3,  cat.cex=4,sub.cex = 0.75,cat.pos=c(352.5,7.5,340,20),main.pos=c(0.5,0.8),
             disable.logging = TRUE,width=6000, height=4000,cat.dist=c(0.35,0.35,0.3,0.3),margin = 0.6)
#####Transcriptome - heatmap####
ann_colors <- list(
  Radiation = c("control" = "grey", "2 Gy" = "#EF3B2C"),
  Cell_line = c(HCC44 =  "#6BAED6", "NCI-H23" = "#41AB5D")
) 

my_colors <-rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))
annotation_col = data.frame(
  Radiation = factor(c((rep("control",3)),rep("2 Gy",4),rep("control",4),rep("2 Gy",4))),
  CellLine = factor(c(rep("NCI-H23",7),rep("HCC44",8))))

rownames(annotation_col) <-c("H23 c1","H23 c2","H23 c3","H23 2Gy1","H23 2Gy2","H23 2Gy3","H23 2Gy4","HCC44 c1","HCC44 c2","HCC44 c3","HCC44 c4","HCC44 2Gy1","HCC44 2Gy2","HCC44 2Gy3","HCC44 2Gy4")

intens_RNA_all_24h <- data.frame(RH23_24h[,1:9],RHCC44_24h[,3:10])
intens_RNA_all_24h<-data.frame(intens_RNA_all_24h[complete.cases(intens_RNA_all_24h),])
colnames(intens_RNA_all_24h)<- c("Genes","ID","H23 c1","H23 c2","H23 c3","H23 2Gy1","H23 2Gy2","H23 2Gy3","H23 2Gy4","HCC44 c1","HCC44 c2","HCC44 c3","HCC44 c4","HCC44 2Gy1","HCC44 2Gy2","HCC44 2Gy3","HCC44 2Gy4")
intens_RNA_all_24h_hm1<- t(intens_RNA_all_24h[,3:17])
intens_RNA_all_24h_hm1<- (apply(intens_RNA_all_24h_hm1, 2, cal_z_score))
Graph_QC_RNA_2<-pheatmap(intens_RNA_all_24h_hm1,fontsize=25, names.arg=intens_RNA_all_24h_hm1[,1],annotation_row= annotation_col, labels_col = colnames(intens_RNA_all_24h_hm1), color=my_colors,annotation_colors = ann_colors, fontsize_row = 3,  fontsize_col = 15, show_rownames=FALSE, show_colnames=FALSE)
save_pheatmap_png(Graph_QC_RNA_2,"Transcriptome_Heatmap.png")
setwd(source_wd)
