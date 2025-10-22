setwd(save_wd)
#####Proteome - intensity#####
box_HCC44_all <- cbind(raw_HCC44_0h_c[,4:7],raw_HCC44_0h_2Gy[,4:7],raw_HCC44_2h_c[,4:7],raw_HCC44_2h_2Gy[,4:7],raw_HCC44_6h_c[,4:7],raw_HCC44_6h_2Gy[,4:7],raw_HCC44_24h_c[,4:7],raw_HCC44_24h_2Gy[,4:7],raw_HCC44_48h_c[,4:7],raw_HCC44_48h_2Gy[,4:7])
box_H23_all <- cbind(raw_H23_0h_c[,4:7],raw_H23_0h_2Gy[,4:7],raw_H23_2h_c[,4:7],raw_H23_2h_2Gy[,4:7],raw_H23_6h_c[,4:7],raw_H23_6h_2Gy[,4:7],raw_H23_24h_c[,4:7],raw_H23_24h_2Gy[,4:7],raw_H23_48h_c[,4:7],raw_H23_48h_2Gy[,4:7])
dat_HCC44_a<-  c("HCC44 0h c","HCC44 0h 2Gy","HCC44 2h c","HCC44 2h 2Gy","HCC44 6h c","HCC44 6h 2Gy","HCC44 24h c","HCC44 24h 2Gy","HCC44 48h c","HCC44 48h 2Gy")
dat_H23_a<-  c("H23 0h c","H23 0h 2Gy","H23 2h c","H23 2h 2Gy","H23 6h c","H23 6h 2Gy","H23 24h c","H23 24h 2Gy","H23 48h c","H23 48h 2Gy")
#####Proteome - plotting intensity####
jpeg(filename = "Proteome_qc_Intensity.jpeg",  width=41, height=21, unit="cm",quality=100,res=100)
layout(matrix(c(2,2,2,2,3,3, 1,1,1,1,3,3),ncol=6,nrow = 2,byrow = TRUE))
par(mar =c(10, 10, 4, 4),xpd=TRUE)
vioplot(log10(box_H23_all),lwd=1.5,axis=FALSE,ann=FALSE,xaxt = "n",yaxt = "n",pchMed=15,cex.axis=1,cex.lab=1,ylim=c(-5,10),boxwex=0.5, col=c(rep("#66C2A5",4),rep("#FDAE61",4)))
box(which = "plot",lwd=2)
axis(1,lwd=2,at=c(1,5,9,13, 17,21,25,29, 33,37), labels=FALSE,lwd.ticks = 2,cex=3, cex.axis=3,tcl = -1)
axis(2,lwd=2,at=c(-5,0,5),lwd.ticks = 2,cex=4, cex.axis=2.5,tcl = -1,las=2,mgp=c(12,2,0))
text(x=1.0,y=-7,labels="H23 0h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=9,y=-7,labels="H23 2h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=17,y=-7,labels="H23 6h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=25,y=-7,labels="H23 24h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=33,y=-7,labels="H23 48h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
title(ylab=(log[10]~"(Intensity)"), cex.lab=2.5,cex.main=3,line=6,adj=c(0.7))
title(main="Proteome - NCI-H23", cex.lab=3,cex.main=3.5,line=1.5,adj=c(0.5))
lines(c(-1,42),c(0,0),lty=3,col="black",lwd=2)
minor.tick(nx =1,ny=2,tick.ratio =1, x.args = list( lwd = 2,tcl = -0.75), y.args = list( lwd = 2,tcl = -0.75)) 
legend(x=c(0,15),y=c(-2.5,-4),adj=c(0,0.5),bty = "n", bg = NA ,ncol=2,c("control","2 Gy"), pch=16, col=c("#66C2A5","#FDAE61"),inset=0.05,pt.cex=2.5, cex =2,box.lwd=2,x.intersp = 0.45,text.width = 0.6)

vioplot(log10(box_HCC44_all),pchMed=15,axis=FALSE,ann=FALSE,lwd=1.5,xaxt = "n",yaxt = "n",cex.axis=4,cex.lab=4,ylim=c(-5,10),boxwex=0.5, col=c(rep("#3288BD",4),rep("#D53E4F",4)))
box(which = "plot",lwd=2)
axis(1,lwd=2,at=c(1,5,9,13, 17,21,25,29, 33,37),labels = FALSE,lwd.ticks = 2,cex=3, cex.axis=3,tcl = -1)
axis(2,lwd=2,at=c(-5,0,5),lwd.ticks = 2,cex=3, cex.axis=2.5,tcl = -1,las=2,mgp=c(12,2,0))
text(x=1,y=-7,labels="HCC44 0h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=9,y=-7,labels="HCC44 2h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=17,y=-7,labels="HCC44 6h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=25,y=-7,labels="HCC44 24h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
text(x=33,y=-7,labels="HCC44 48h",srt = 35,cex = 2,xpd=TRUE,adj=c(1,0.5))
title(ylab=(log[10]~"(Intensity)"), cex.lab=2.5,cex.main=3,line=6,adj=c(0.7))
title(main="Proteome - HCC44", cex.lab=2.5,cex.main=3.5,line=1.5,adj=c(0.5))
lines(c(-1,42),c(0,0),lty=3,col="black",lwd=2)
minor.tick(nx =1,ny=2,tick.ratio =1, x.args = list( lwd =2,tcl = -0.75), y.args = list( lwd = 2,tcl = -0.75)) 
legend(x=c(0,15),y=c(-2.5,-4),adj=c(0,0.5),bty = "n", bg = NA ,ncol=2,c("control","2 Gy"), pch=16, col=c("#3288BD","#D53E4F"),inset=0.05,pt.cex=2.5, cex =2,box.lwd=2,x.intersp = 0.45,text.width = 0.6)
dev.off()
#####Proteome - PCA####
pbase_pca_filt<- data.frame(log2(raw_H23_0h_c[,c(4:7)]),log2(raw_H23_0h_2Gy[,c(4:7)]),log2(raw_HCC44_0h_c[,c(4:7)]),log2(raw_HCC44_0h_2Gy[,c(4:7)]))
pbase_pca_filt <- pbase_pca_filt %>% mutate_all(~ifelse(is.nan(.), NA, .))
pbase_pca_filt<-t(pbase_pca_filt[complete.cases(pbase_pca_filt),])
pc_pbase_pca_filt <- prcomp(pbase_pca_filt,
                            center = TRUE,
                            scale = TRUE)
#####Proteome - PCA plot#####
jpeg(filename = "Proteome_qc_PCA.jpeg",  width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2,3,3,3),ncol=6,nrow = 1,byrow = TRUE))
par(mar = c(10, 10, 6, 4),xpd=TRUE)
leg.txt <- c("HCC44 c", "HCC44 2 Gy","NCI-H23 c","NCI-H23 2 Gy")
s_pc_pbase_pca_filt <- summary(pc_pbase_pca_filt)
plot(as.data.frame(s_pc_pbase_pca_filt["x"])[,1],ann=FALSE, as.data.frame(s_pc_pbase_pca_filt["x"])[,2],pch=21, xlab=paste("PCA 1 (", round(s_pc_pbase_pca_filt$importance[2]*100, 1), "%)", sep = ""), 
     ylab=paste("PCA 2 (", round(s_pc_pbase_pca_filt$importance[5]*100, 1), "%)", sep = ""),bg=c(rep("#66C2A5",4),rep("#FDAE61",4),rep("#3288BD",4),rep("#D53E4F",4)),  col="black",
     ,cex=5,cex.main=3,cex.lab=4,xlim=c(-125,125),ylim=c(-125,125), cex.axis=4,axes=FALSE,mgp=c(8,0,0),lwd=2,main="Proteome - Control")
minor.tick(nx =2,ny=2,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
axis(1,lwd=3,at = c(-100,0,100),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(6,3,0))
axis(2,lwd=3,at = c(-100,0,100),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(6,2,0),las=2)
lines(c(0,0),c(-125,125), lty=3,lwd=2)
lines(c(-125,125),c(0,0), lty=3,lwd=2)
title(main="Proteome (0h after radiation)",cex.lab=3,cex.main=3,line=2,adj=c(0.5))
title(xlab=paste("PCA 1 (", round(s_pc_pbase_pca_filt$importance[2]*100, 1), "%)", sep = ""),cex.lab=3.5,cex.main=3,line=6,adj=c(0.5))
title(ylab=paste("PCA 2 (", round(s_pc_pbase_pca_filt$importance[5]*100, 1), "%)", sep = ""),cex.lab=3.5,cex.main=3,line=7,adj=c(0.5))
box(lwd=3,lty="solid")
#legend(x=c(-135,125),y=c(-150,-70),ncol=2, box.col="white" ,bty = "n", bg = NA ,leg.txt,pch=21,cex =2.25,pt.cex = 3,pt.bg = c("#3288BD","#D53E4F","#66C2A5","#FDAE61"),  col="black",box.lwd =2,pt.lwd = 2,y.intersp=1,x.intersp=0.75)

par(mar = c(10, 0, 6, 4),xpd=TRUE)
plot(1:10,1:10,ann=FALSE,axes=FALSE,col="white")
legend(x=c(0,8),y=c(2,10.25),ncol = 1,adj=c(0,0.25),legend=leg.txt ,bty = "n", bg = NA , pch=21, pt.bg= c("#3288BD","#D53E4F","#66C2A5","#FDAE61"),  col="black",pt.cex=3, cex = 2.75,box.lwd=2,y.intersp = 1.3)
dev.off()
#####Proteome - COV #####
pHCC44_0h_2Gy_cov<- raw_HCC44_0h_2Gy[,1:2]
pHCC44_0h_2Gy_cov<- cbind(pHCC44_0h_2Gy_cov,(raw_HCC44_0h_2Gy[,4:7]))
pHCC44_0h_2Gy_cov<-pHCC44_0h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_0h_2Gy_cov[["mean"]]<-rowMeans(pHCC44_0h_2Gy_cov[,3:6], na.rm=TRUE)
pHCC44_0h_2Gy_cov[["SD"]] <-apply(pHCC44_0h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_0h_2Gy_cov[["%CV"]] <-(pHCC44_0h_2Gy_cov[["SD"]]/pHCC44_0h_2Gy_cov[["mean"]])*100
dens_pHCC44_0h_2Gy_cov<-density(pHCC44_0h_2Gy_cov[,9],na.rm=TRUE)

pHCC44_2h_2Gy_cov<- raw_HCC44_2h_2Gy[,1:2]
pHCC44_2h_2Gy_cov<- cbind(pHCC44_2h_2Gy_cov,(raw_HCC44_2h_2Gy[,4:7]))
pHCC44_2h_2Gy_cov<-pHCC44_2h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_2h_2Gy_cov[["mean"]]<-rowMeans(pHCC44_2h_2Gy_cov[,3:6], na.rm=TRUE)
pHCC44_2h_2Gy_cov[["SD"]] <-apply(pHCC44_2h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_2h_2Gy_cov[["%CV"]] <-(pHCC44_2h_2Gy_cov[["SD"]]/pHCC44_2h_2Gy_cov[["mean"]])*100
dens_pHCC44_2h_2Gy_cov<-density(pHCC44_2h_2Gy_cov[,9],na.rm=TRUE)

pHCC44_6h_2Gy_cov<- raw_HCC44_6h_2Gy[,1:2]
pHCC44_6h_2Gy_cov<- cbind(pHCC44_6h_2Gy_cov,(raw_HCC44_6h_2Gy[,4:7]))
pHCC44_6h_2Gy_cov<-pHCC44_6h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_6h_2Gy_cov[["mean"]]<-rowMeans(pHCC44_6h_2Gy_cov[,3:6], na.rm=TRUE)
pHCC44_6h_2Gy_cov[["SD"]] <-apply(pHCC44_6h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_6h_2Gy_cov[["%CV"]] <-(pHCC44_6h_2Gy_cov[["SD"]]/pHCC44_6h_2Gy_cov[["mean"]])*100
dens_pHCC44_6h_2Gy_cov<-density(pHCC44_6h_2Gy_cov[,9],na.rm=TRUE)

pHCC44_24h_2Gy_cov<- raw_HCC44_24h_2Gy[,1:2]
pHCC44_24h_2Gy_cov<- cbind(pHCC44_24h_2Gy_cov,(raw_HCC44_24h_2Gy[,4:7]))
pHCC44_24h_2Gy_cov<-pHCC44_24h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_24h_2Gy_cov[["mean"]]<-rowMeans(pHCC44_24h_2Gy_cov[,3:6], na.rm=TRUE)
pHCC44_24h_2Gy_cov[["SD"]] <-apply(pHCC44_24h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_24h_2Gy_cov[["%CV"]] <-(pHCC44_24h_2Gy_cov[["SD"]]/pHCC44_24h_2Gy_cov[["mean"]])*100
dens_pHCC44_24h_2Gy_cov<-density(pHCC44_24h_2Gy_cov[,9],na.rm=TRUE)

pHCC44_48h_2Gy_cov<- raw_HCC44_48h_2Gy[,1:2]
pHCC44_48h_2Gy_cov<- cbind(pHCC44_48h_2Gy_cov,(raw_HCC44_48h_2Gy[,4:7]))
pHCC44_48h_2Gy_cov<-pHCC44_48h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_48h_2Gy_cov[["mean"]]<-rowMeans(pHCC44_48h_2Gy_cov[,3:6], na.rm=TRUE)
pHCC44_48h_2Gy_cov[["SD"]] <-apply(pHCC44_48h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_48h_2Gy_cov[["%CV"]] <-(pHCC44_48h_2Gy_cov[["SD"]]/pHCC44_48h_2Gy_cov[["mean"]])*100
dens_pHCC44_48h_2Gy_cov<-density(pHCC44_48h_2Gy_cov[,9],na.rm=TRUE)

pHCC44_0h_c_cov<- raw_HCC44_0h_c[,1:2]
pHCC44_0h_c_cov<- cbind(pHCC44_0h_c_cov,(raw_HCC44_0h_c[,4:7]))
pHCC44_0h_c_cov<-pHCC44_0h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_0h_c_cov[["mean"]]<-rowMeans(pHCC44_0h_c_cov[,3:6], na.rm=TRUE)
pHCC44_0h_c_cov[["SD"]] <-apply(pHCC44_0h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_0h_c_cov[["%CV"]] <-(pHCC44_0h_c_cov[["SD"]]/pHCC44_0h_c_cov[["mean"]])*100
dens_pHCC44_0h_c_cov<-density(pHCC44_0h_c_cov[,9],na.rm=TRUE)

pHCC44_2h_c_cov<- raw_HCC44_2h_c[,1:2]
pHCC44_2h_c_cov<- cbind(pHCC44_2h_c_cov,(raw_HCC44_2h_c[,4:7]))
pHCC44_2h_c_cov<-pHCC44_2h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_2h_c_cov[["mean"]]<-rowMeans(pHCC44_2h_c_cov[,3:6], na.rm=TRUE)
pHCC44_2h_c_cov[["SD"]] <-apply(pHCC44_2h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_2h_c_cov[["%CV"]] <-(pHCC44_2h_c_cov[["SD"]]/pHCC44_2h_c_cov[["mean"]])*100
dens_pHCC44_2h_c_cov<-density(pHCC44_2h_c_cov[,9],na.rm=TRUE)

pHCC44_6h_c_cov<- raw_HCC44_6h_c[,1:2]
pHCC44_6h_c_cov<- cbind(pHCC44_6h_c_cov,(raw_HCC44_6h_c[,4:7]))
pHCC44_6h_c_cov<-pHCC44_6h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_6h_c_cov[["mean"]]<-rowMeans(pHCC44_6h_c_cov[,3:6], na.rm=TRUE)
pHCC44_6h_c_cov[["SD"]] <-apply(pHCC44_6h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_6h_c_cov[["%CV"]] <-(pHCC44_6h_c_cov[["SD"]]/pHCC44_6h_c_cov[["mean"]])*100
dens_pHCC44_6h_c_cov<-density(pHCC44_6h_c_cov[,9],na.rm=TRUE)

pHCC44_24h_c_cov<- raw_HCC44_24h_c[,1:2]
pHCC44_24h_c_cov<- cbind(pHCC44_24h_c_cov,(raw_HCC44_24h_c[,4:7]))
pHCC44_24h_c_cov<-pHCC44_24h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_24h_c_cov[["mean"]]<-rowMeans(pHCC44_24h_c_cov[,3:6], na.rm=TRUE)
pHCC44_24h_c_cov[["SD"]] <-apply(pHCC44_24h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_24h_c_cov[["%CV"]] <-(pHCC44_24h_c_cov[["SD"]]/pHCC44_24h_c_cov[["mean"]])*100
dens_pHCC44_24h_c_cov<-density(pHCC44_24h_c_cov[,9],na.rm=TRUE)

pHCC44_48h_c_cov<- raw_HCC44_48h_c[,1:2]
pHCC44_48h_c_cov<- cbind(pHCC44_48h_c_cov,(raw_HCC44_48h_c[,4:7]))
pHCC44_48h_c_cov<-pHCC44_48h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pHCC44_48h_c_cov[["mean"]]<-rowMeans(pHCC44_48h_c_cov[,3:6], na.rm=TRUE)
pHCC44_48h_c_cov[["SD"]] <-apply(pHCC44_48h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pHCC44_48h_c_cov[["%CV"]] <-(pHCC44_48h_c_cov[["SD"]]/pHCC44_48h_c_cov[["mean"]])*100
dens_pHCC44_48h_c_cov<-density(pHCC44_48h_c_cov[,9],na.rm=TRUE)

pH23_0h_2Gy_cov<- raw_H23_0h_2Gy[,1:2]
pH23_0h_2Gy_cov<- cbind(pH23_0h_2Gy_cov,(raw_H23_0h_2Gy[,4:7]))
pH23_0h_2Gy_cov<-pH23_0h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_0h_2Gy_cov[["mean"]]<-rowMeans(pH23_0h_2Gy_cov[,3:6], na.rm=TRUE)
pH23_0h_2Gy_cov[["SD"]] <-apply(pH23_0h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_0h_2Gy_cov[["%CV"]] <-(pH23_0h_2Gy_cov[["SD"]]/pH23_0h_2Gy_cov[["mean"]])*100
dens_pH23_0h_2Gy_cov<-density(pH23_0h_2Gy_cov[,9],na.rm=TRUE)

pH23_2h_2Gy_cov<- raw_H23_2h_2Gy[,1:2]
pH23_2h_2Gy_cov<- cbind(pH23_2h_2Gy_cov,(raw_H23_2h_2Gy[,4:7]))
pH23_2h_2Gy_cov<-pH23_2h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_2h_2Gy_cov[["mean"]]<-rowMeans(pH23_2h_2Gy_cov[,3:6], na.rm=TRUE)
pH23_2h_2Gy_cov[["SD"]] <-apply(pH23_2h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_2h_2Gy_cov[["%CV"]] <-(pH23_2h_2Gy_cov[["SD"]]/pH23_2h_2Gy_cov[["mean"]])*100
dens_pH23_2h_2Gy_cov<-density(pH23_2h_2Gy_cov[,9],na.rm=TRUE)

pH23_6h_2Gy_cov<- raw_H23_6h_2Gy[,1:2]
pH23_6h_2Gy_cov<- cbind(pH23_6h_2Gy_cov,(raw_H23_6h_2Gy[,4:7]))
pH23_6h_2Gy_cov<-pH23_6h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_6h_2Gy_cov[["mean"]]<-rowMeans(pH23_6h_2Gy_cov[,3:6], na.rm=TRUE)
pH23_6h_2Gy_cov[["SD"]] <-apply(pH23_6h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_6h_2Gy_cov[["%CV"]] <-(pH23_6h_2Gy_cov[["SD"]]/pH23_6h_2Gy_cov[["mean"]])*100
dens_pH23_6h_2Gy_cov<-density(pH23_6h_2Gy_cov[,9],na.rm=TRUE)

pH23_24h_2Gy_cov<- raw_H23_24h_2Gy[,1:2]
pH23_24h_2Gy_cov<- cbind(pH23_24h_2Gy_cov,(raw_H23_24h_2Gy[,4:7]))
pH23_24h_2Gy_cov<-pH23_24h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_24h_2Gy_cov[["mean"]]<-rowMeans(pH23_24h_2Gy_cov[,3:6], na.rm=TRUE)
pH23_24h_2Gy_cov[["SD"]] <-apply(pH23_24h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_24h_2Gy_cov[["%CV"]] <-(pH23_24h_2Gy_cov[["SD"]]/pH23_24h_2Gy_cov[["mean"]])*100
dens_pH23_24h_2Gy_cov<-density(pH23_24h_2Gy_cov[,9],na.rm=TRUE)

pH23_48h_2Gy_cov<- raw_H23_48h_2Gy[,1:2]
pH23_48h_2Gy_cov<- cbind(pH23_48h_2Gy_cov,(raw_H23_48h_2Gy[,4:7]))
pH23_48h_2Gy_cov<-pH23_48h_2Gy_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_48h_2Gy_cov[["mean"]]<-rowMeans(pH23_48h_2Gy_cov[,3:6], na.rm=TRUE)
pH23_48h_2Gy_cov[["SD"]] <-apply(pH23_48h_2Gy_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_48h_2Gy_cov[["%CV"]] <-(pH23_48h_2Gy_cov[["SD"]]/pH23_48h_2Gy_cov[["mean"]])*100
dens_pH23_48h_2Gy_cov<-density(pH23_48h_2Gy_cov[,9],na.rm=TRUE)

pH23_0h_c_cov<- raw_H23_0h_c[,1:2]
pH23_0h_c_cov<- cbind(pH23_0h_c_cov,(raw_H23_0h_c[,4:7]))
pH23_0h_c_cov<-pH23_0h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_0h_c_cov[["mean"]]<-rowMeans(pH23_0h_c_cov[,3:6], na.rm=TRUE)
pH23_0h_c_cov[["SD"]] <-apply(pH23_0h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_0h_c_cov[["%CV"]] <-(pH23_0h_c_cov[["SD"]]/pH23_0h_c_cov[["mean"]])*100
dens_pH23_0h_c_cov<-density(pH23_0h_c_cov[,9],na.rm=TRUE)

pH23_2h_c_cov<- raw_H23_2h_c[,1:2]
pH23_2h_c_cov<- cbind(pH23_2h_c_cov,(raw_H23_2h_c[,4:7]))
pH23_2h_c_cov<-pH23_2h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_2h_c_cov[["mean"]]<-rowMeans(pH23_2h_c_cov[,3:6], na.rm=TRUE)
pH23_2h_c_cov[["SD"]] <-apply(pH23_2h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_2h_c_cov[["%CV"]] <-(pH23_2h_c_cov[["SD"]]/pH23_2h_c_cov[["mean"]])*100
dens_pH23_2h_c_cov<-density(pH23_2h_c_cov[,9],na.rm=TRUE)

pH23_6h_c_cov<- raw_H23_6h_c[,1:2]
pH23_6h_c_cov<- cbind(pH23_6h_c_cov,(raw_H23_6h_c[,4:7]))
pH23_6h_c_cov<-pH23_6h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_6h_c_cov[["mean"]]<-rowMeans(pH23_6h_c_cov[,3:6], na.rm=TRUE)
pH23_6h_c_cov[["SD"]] <-apply(pH23_6h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_6h_c_cov[["%CV"]] <-(pH23_6h_c_cov[["SD"]]/pH23_6h_c_cov[["mean"]])*100
dens_pH23_6h_c_cov<-density(pH23_6h_c_cov[,9],na.rm=TRUE)

pH23_24h_c_cov<- raw_H23_24h_c[,1:2]
pH23_24h_c_cov<- cbind(pH23_24h_c_cov,(raw_H23_24h_c[,4:7]))
pH23_24h_c_cov<-pH23_24h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_24h_c_cov[["mean"]]<-rowMeans(pH23_24h_c_cov[,3:6], na.rm=TRUE)
pH23_24h_c_cov[["SD"]] <-apply(pH23_24h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_24h_c_cov[["%CV"]] <-(pH23_24h_c_cov[["SD"]]/pH23_24h_c_cov[["mean"]])*100
dens_pH23_24h_c_cov<-density(pH23_24h_c_cov[,9],na.rm=TRUE)

pH23_48h_c_cov<- raw_H23_48h_c[,1:2]
pH23_48h_c_cov<- cbind(pH23_48h_c_cov,(raw_H23_48h_c[,4:7]))
pH23_48h_c_cov<-pH23_48h_c_cov%>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
pH23_48h_c_cov[["mean"]]<-rowMeans(pH23_48h_c_cov[,3:6], na.rm=TRUE)
pH23_48h_c_cov[["SD"]] <-apply(pH23_48h_c_cov[,3:6], 1, sd, na.rm=TRUE)
pH23_48h_c_cov[["%CV"]] <-(pH23_48h_c_cov[["SD"]]/pH23_48h_c_cov[["mean"]])*100
dens_pH23_48h_c_cov<-density(pH23_48h_c_cov[,9],na.rm=TRUE)

#####Proteome - plotting CoV####
jpeg(filename = "Proteome_qc_CoV_HCC44.jpeg",  width=20.5, height=28, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1,2,2,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt<-c("0h","2h","6h","24h","48h")
plot(dens_pHCC44_0h_c_cov$x,length(raw_HCC44_0h_c[,4])*dens_pHCC44_0h_c_cov$y,ylim=c(0,400),xlim=c(-10,215),type="l",lwd=3,
     xlab="Coefficient of Variation[%]",ylab="Density", col=color[1,],ann=FALSE,axes=FALSE,cex.lab=4,mgp=c(6,0,1))
lines(dens_pHCC44_2h_c_cov$x,length(raw_HCC44_2h_c[,4])*dens_pHCC44_2h_c_cov$y, col=color[3,],lwd=3)
lines(dens_pHCC44_6h_c_cov$x,length(raw_HCC44_6h_c[,4])*dens_pHCC44_6h_c_cov$y, col=color[5,],lwd=3)
lines(dens_pHCC44_24h_c_cov$x,length(raw_HCC44_24h_c[,4])*dens_pHCC44_24h_c_cov$y, col=color[7,],lwd=3)
lines(dens_pHCC44_48h_c_cov$x,length(raw_HCC44_48h_c[,4])*dens_pHCC44_48h_c_cov$y, col=color[9,],lwd=3)
axis(1,lwd=3,at = c(0,100,200),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,3,0))
axis(2,lwd=3,at = c(0,200,400),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,2,0),las=2)
minor.tick(nx =2/1,ny=2/1,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
lines(c(15,15),c(-10,410), lty=3,lwd=2)
box(which = "plot",type="solid",lwd=3)
title(ylab="Density", cex.lab=3,cex.main=3,line=8,adj=c(0.5))
title(xlab="Coefficient of Variation[%]", cex.lab=3,cex.main=3,line=6,)
title(main="HCC44 - Control", cex.lab=3,cex.main=3,line=2,adj=c(0.5))
legend(x=c(140,215),y=c(160,400),leg.txt,ncol=1,lty=1,cex=2.75,lwd=4,  col=color[c(1,3,5,7,9),],inset=0.02,box.lwd =3,x.intersp = 0.7)

plot(dens_pHCC44_0h_2Gy_cov$x,length(raw_HCC44_0h_2Gy[,4])*dens_pHCC44_0h_2Gy_cov$y,ylim=c(0,400),xlim=c(-10,215),type="l",lwd=3,
     xlab="Coefficient of Variation[%]",ylab="Density", col=color[1,],ann=FALSE,axes=FALSE,cex.lab=4,mgp=c(6,0,1))
lines(dens_pHCC44_2h_2Gy_cov$x,length(raw_HCC44_2h_2Gy[,4])*dens_pHCC44_2h_2Gy_cov$y, col=color[3,],lwd=3)
lines(dens_pHCC44_6h_2Gy_cov$x,length(raw_HCC44_6h_2Gy[,4])*dens_pHCC44_6h_2Gy_cov$y, col=color[5,],lwd=3)
lines(dens_pHCC44_24h_2Gy_cov$x,length(raw_HCC44_24h_2Gy[,4])*dens_pHCC44_24h_2Gy_cov$y, col=color[7,],lwd=3)
lines(dens_pHCC44_48h_2Gy_cov$x,length(raw_HCC44_48h_2Gy[,4])*dens_pHCC44_48h_2Gy_cov$y, col=color[9,],lwd=3)
axis(1,lwd=3,at = c(0,100,200),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,3,0))
axis(2,lwd=3,at = c(0,200,400),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,2,0),las=2)
minor.tick(nx =2/1,ny=2/1,tick.ratio =1, x.args = list( lwd = 3,tcl=-0.75), y.args = list( lwd = 3,tcl=-0.75)) 
lines(c(15,15),c(-10,410), lty=3,lwd=2)
box(which = "plot",type="solid",lwd=3)
title(ylab="Density", cex.lab=3,cex.main=3,line=8,adj=c(0.5))
title(xlab="Coefficient of Variation[%]", cex.lab=3,cex.main=3,line=6,)
title(main="HCC44 - Radiation", cex.lab=3,cex.main=3,line=2,adj=c(0.5))
#legend(x=c(45,250),y=c(235,450),leg.txt,ncol=1,lty=1,cex=2.75,lwd=4,  col=c("dodgerblue","red","darkblue","darkred"),inset=0.02,box.lwd =3,x.intersp = 0.7)
dev.off()

jpeg(filename = "Proteome_qc_CoV_H23.jpeg",  width=20.5, height=28, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1,2,2,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt<-c("0h","2h","6h","24h","48h")
plot(dens_pH23_0h_c_cov$x,length(raw_H23_0h_c[,4])*dens_pH23_0h_c_cov$y,ylim=c(0,450),xlim=c(-10,215),type="l",lwd=3,
     xlab="Coefficient of Variation[%]",ylab="Density",ann=FALSE, col=color[1,],axes=FALSE,cex.lab=4,mgp=c(6,0,1))
lines(dens_pH23_2h_c_cov$x,length(raw_H23_2h_c[,4])*dens_pH23_2h_c_cov$y, col=color[3,],lwd=3)
lines(dens_pH23_6h_c_cov$x,length(raw_H23_6h_c[,4])*dens_pH23_6h_c_cov$y, col=color[5,],lwd=3)
lines(dens_pH23_24h_c_cov$x,length(raw_H23_24h_c[,4])*dens_pH23_24h_c_cov$y, col=color[7,],lwd=3)
lines(dens_pH23_48h_c_cov$x,length(raw_H23_48h_c[,4])*dens_pH23_48h_c_cov$y, col=color[9,],lwd=3)
axis(1,lwd=3,at = c(0,100,200),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,mgp=c(4,3,0))
axis(2,lwd=3,at = c(0,150,300,450),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,mgp=c(4,2,0),las=2)
minor.tick(nx =2/1,ny=2/1,tick.ratio =0.75, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
lines(c(15,15),c(-10,470), lty=3,lwd=2)
box(which = "plot",type="solid",lwd=3)
title(ylab="Density", cex.lab=3,cex.main=3,line=8,adj=c(0.5))
title(xlab="Coefficient of Variation[%]", cex.lab=3,cex.main=3,line=6,)
title(main="NCI-H23 - Control", cex.lab=3,cex.main=3,line=2,adj=c(0.5))
legend(x=c(140,215),y=c(180,450),leg.txt,ncol=1,lty=1,cex=2.75,lwd=4,  col=color[c(1,3,5,7,9),],inset=0.02,box.lwd =3,x.intersp = 0.7)

plot(dens_pH23_0h_2Gy_cov$x,length(raw_H23_0h_2Gy[,4])*dens_pH23_0h_2Gy_cov$y,ylim=c(0,450),xlim=c(-10,215),type="l",lwd=3,
     xlab="Coefficient of Variation[%]",ylab="Density",ann=FALSE, col=color[1,],axes=FALSE,cex.lab=4,mgp=c(6,0,1))
lines(dens_pH23_2h_2Gy_cov$x,length(raw_H23_2h_2Gy[,4])*dens_pH23_2h_2Gy_cov$y, col=color[3,],lwd=3)
lines(dens_pH23_6h_2Gy_cov$x,length(raw_H23_6h_2Gy[,4])*dens_pH23_6h_2Gy_cov$y, col=color[5,],lwd=3)
lines(dens_pH23_24h_2Gy_cov$x,length(raw_H23_24h_2Gy[,4])*dens_pH23_24h_2Gy_cov$y, col=color[7,],lwd=3)
lines(dens_pH23_48h_2Gy_cov$x,length(raw_H23_48h_2Gy[,4])*dens_pH23_48h_2Gy_cov$y, col=color[91,],lwd=3)
axis(1,lwd=3,at = c(0,100,200),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,3,0))
axis(2,lwd=3,at = c(0,150,300,450),lwd.ticks = 3,cex=2.5, cex.axis=2.5,tcl = -1.5,mgp=c(4,2,0),las=2)
minor.tick(nx =2/1,ny=2/1,tick.ratio =0.75, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
lines(c(15,15),c(-10,470), lty=3,lwd=2)
box(which = "plot",type="solid",lwd=3)
title(ylab="Density", cex.lab=3,cex.main=3,line=8,adj=c(0.7))
title(xlab="Coefficient of Variation[%]", cex.lab=3,cex.main=3,line=6,)
title(main="NCI-H23 - Radiation", cex.lab=3,cex.main=3,line=2,adj=c(0.5))
#legend(x=c(45,250),y=c(235,450),leg.txt,ncol=1,lty=1,cex=2.75,lwd=4,  col=c("limegreen","gold","darkgreen","darkgoldenrod3"),inset=0.02,box.lwd =3,x.intersp = 0.7)
#leg.txt<-c("0h","2h","6h","24h","48h")
#par(mar =c(2, 1, 2, 1),xpd=TRUE)
#plot(1:10,1:10,ann=FALSE,axes=FALSE, col="white")
#legend(x=c(0,10),y=c(1.5,10),leg.txt,ncol=1,title="HCC44",lty=1,cex=2.75,lwd=4,  col=c("#C6DBEF",cols_HCC44[2:5]),inset=0.02,box.lwd =3,x.intersp = 0.7)
#plot(1:10,1:10,ann=FALSE,axes=FALSE, col="white")
#legend(x=c(0,10),y=c(1.5,10),leg.txt,ncol=1,title="NCI-H23",lty=1,cex=2.75,lwd=4,  col=c("#C7E9C0",cols_Con[2:5]),inset=0.02,box.lwd =3,x.intersp = 0.7)
dev.off()
#####Proteome - identified Proteins#####
n_pcount<-NULL
n_pcount["H23_0h_c"]<-as.numeric(data.frame(count(raw_H23_0h_c[complete.cases(raw_H23_0h_c$mean),])))
n_pcount["H23_0h_2Gy"]<-as.numeric(data.frame(count(data.frame((raw_H23_0h_2Gy[complete.cases(raw_H23_0h_2Gy$mean),])))))
n_pcount["H23_24h_c"]<-as.numeric(data.frame(count(data.frame((raw_H23_24h_c[complete.cases(raw_H23_24h_c$mean),])))))
n_pcount["H23_24h_2Gy"]<-as.numeric(data.frame(count(data.frame((raw_H23_24h_2Gy[complete.cases(raw_H23_24h_2Gy$mean),])))))
n_pcount["HCC44_0h_c"]<-as.numeric(data.frame(count(data.frame((raw_HCC44_0h_c[complete.cases(raw_HCC44_0h_c$mean),])))))
n_pcount["HCC44_0h_2Gy"]<-as.numeric(data.frame(count(data.frame((raw_HCC44_0h_2Gy[complete.cases(raw_HCC44_0h_2Gy$mean),])))))
n_pcount["HCC44_24h_c"]<-as.numeric(data.frame(count(data.frame((raw_HCC44_24h_c[complete.cases(raw_HCC44_24h_c$mean),])))))
n_pcount["HCC44_24h_2Gy"]<-as.numeric(data.frame(count(data.frame((raw_HCC44_24h_2Gy[complete.cases(raw_HCC44_24h_2Gy$mean),])))))



#####Proteome - plotting identified proteins#####
jpeg(filename = "Proteome_qc_IndentifiedProteins.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar = c(14, 14, 4, 6),xpd=TRUE)
barplot(as.numeric(n_pcount[c(5:8,1:4)]), col=c("#3288BD","#D53E4F","#5E4FA2","#9E0142","#ABDDA4", "#FEE08B","#66C2A5","#FDAE61"),axes=FALSE,ylim=c(0,11000))
minor.tick(nx =1,ny=2,tick.ratio =1, x.args = list( lwd = 3), y.args = list( lwd = 3,tcl = -1)) 
axis(2,lwd=3,at = c(0,5000,10000),labels=c(0,"5.000","10.000"),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,las=2,mgp=c(6,2,0))
axis(1,lwd=3,at = c(0,10),labels=c("",""),lwd.ticks = 3,cex=4, cex.axis=4,tcl = 0)
title( ylab="detected proteins", cex.lab=3,cex.main=3,line=10,adj=1)
text(x = 7,y = -600,adj=c(1,0),labels =  "NCI-H23 0h c",srt = 45,cex =2.25,col="black")
text(x = 8.2,y = -600,adj=c(1,0),labels ="NCI-H23 0h 2Gy"  ,srt = 45,cex =2.25,col="black")
text(x = 9.45,y = -600,adj=c(1,0),labels ="NCI-H23 24h c"  ,srt = 45,cex =2.25,col="black")
text(x = 5.75,y = -600,adj=c(1,0),labels ="NCI-H23 24h 2Gy"  ,srt = 45,cex =2.25,col="black")
text(x = 0.7,y = -600,adj=c(1,0),labels ="HCC44 0h c"  ,srt = 45,cex =2.25,col="black")
text(x = 1.9, y = -600,adj=c(1,0),labels ="HCC44 0h 2Gy"  ,srt = 45,cex =2.25,col="black")
text(x = 3.2,y = -600,adj=c(1,0),labels ="HCC44 24h c"  ,srt = 45,cex =2.25,col="black")
text(x = 4.4,y =-600,adj=c(1,0),labels ="HCC44 24h 2Gy"  ,srt = 45,cex =2.25,col="black")
text(x = c(0.7),y = 9900,adj=c(0,0),labels =n_pcount[5],srt = 45,cex =2.25,col="black")
text(x = c(1.9),y = 9900,adj=c(0,0),labels =n_pcount[6],srt = 45,cex =2.25,col="black")
text(x = c(3.1),y = 9900,adj=c(0,0),labels =n_pcount[7],srt = 45,cex =2.25,col="black")
text(x = c(4.3),y = 9900,adj=c(0,0),labels =n_pcount[8],srt = 45,cex =2.25,col="black")
text(x = c(5.5),y = 9900,adj=c(0,0),labels =n_pcount[1],srt = 45,cex =2.25,col="black")
text(x = c(6.75),y = 9900,adj=c(0,0),labels =n_pcount[2],srt = 45,cex =2.25,col="black")
text(x = c(7.9),y = 9900,adj=c(0,0),labels =n_pcount[3],srt = 45,cex =2.25,col="black")
text(x = c(9.1),y = 9900,adj=c(0,0),labels =n_pcount[4],srt = 45,cex =2.25,col="black")
text(x = c(9.9),y = 7800,adj=c(0,0),labels ="7557",srt = 0,cex = 2.5,col="black")
lines(c(0,11),c(7557,7557), lty=3,lwd=3,col="black")
dev.off()
#####Proteome - heatmap####
annotation_col = data.frame(
  Radiation = factor(rep(c((rep("control",4)),rep("2 Gy",4)),2)),
  Cell_line = factor(c(rep("NCI-H23",8),rep("HCC44",8)))
)
rownames(annotation_col) <-c("NCI-H23 c 0h.1","NCI-H23 c 0h.2","NCI-H23 c 0h.3","NCI-H23 c 0h.4","NCI-H23 2Gy 0h.1","NCI-H23 2Gy 0h.2","NCI-H23 2Gy 0h.3","NCI-H23 2Gy 0h.4",
                             "HCC44 c 0h.1","HCC44 c 0h.2","HCC44 c 0h.3","HCC44 c 0h.4","HCC44 2Gy 0h.1","HCC44 2Gy 0h.2","HCC44 2Gy 0h.3","HCC44 2Gy 0h.4")
pH23_pHCC44_0h_all_filt<- data.frame(log2(raw_H23_0h_c[,c(5:8)]),log2(raw_H23_0h_2Gy[,c(5:8)]),log2(raw_HCC44_0h_c[,c(5:8)]),log2(raw_HCC44_0h_2Gy[,c(5:8)]))  
pH23_pHCC44_0h_all_filt <- pH23_pHCC44_0h_all_filt %>% mutate_all(~ifelse(is.nan(.), NA, .))
pH23_pHCC44_0h_all_filt <- cbind( H23_0h[,2:3],pH23_pHCC44_0h_all_filt)
pH23_pHCC44_0h_all_filt<-(pH23_pHCC44_0h_all_filt[complete.cases(pH23_pHCC44_0h_all_filt[,3:18]),])
colnames(pH23_pHCC44_0h_all_filt)<- c("ID","Suggested.Symbol","NCI-H23 c 0h.1","NCI-H23 c 0h.2","NCI-H23 c 0h.3","NCI-H23 c 0h.4","NCI-H23 2Gy 0h.1","NCI-H23 2Gy 0h.2","NCI-H23 2Gy 0h.3","NCI-H23 2Gy 0h.4",
                                      "HCC44 c 0h.1","HCC44 c 0h.2","HCC44 c 0h.3","HCC44 c 0h.4","HCC44 2Gy 0h.1","HCC44 2Gy 0h.2","HCC44 2Gy 0h.3","HCC44 2Gy 0h.4")
ready_pH23_pHCC44_0h_intens_mean_hm1<- (apply(pH23_pHCC44_0h_all_filt[,3:18], 1, cal_z_score))
Graph_QC_Proteom_2<-pheatmap(ready_pH23_pHCC44_0h_intens_mean_hm1,fontsize=25, names.arg=pH23_pHCC44_0h_all_filt[,2],annotation_row = annotation_col, labels_row = colnames(pH23_pHCC44_0h_all_filt[,3:18]),color=my_colors,annotation_colors = ann_colors, fontsize_row = 3,  fontsize_col = 15, show_colnames = FALSE, show_rownames=FALSE)
save_pheatmap_png(Graph_QC_Proteom_2,"Proteome_qc_Heatmap.png")
#####Proteome - Venn diagramm####
venn.diagram(list(raw_H23_0h_c[complete.cases(raw_H23_0h_c[,c(8)]),][,3],raw_H23_0h_2Gy[complete.cases(raw_H23_0h_2Gy[,c(8)]),][,3],raw_HCC44_0h_c[complete.cases(raw_HCC44_0h_c[,c(8)]),][,3],raw_HCC44_0h_2Gy[complete.cases(raw_HCC44_0h_2Gy[,c(8)]),][,3]),
             category.names = c("H23 0h c","H23 0h 2Gy","HCC44 0h c","HCC44 0h 2Gy"),
             filename = "Proteome_qc_vennDiagramm.jpeg",
             print.mode = "raw",na="none",cat.fontfamily=1,main.fontfamily= 1,
             fill = c("#ABDDA4","#FDAE61","#3288BD","#D53E4F"),ext.line.lwd=4,ext.percent=c(2,2,0), main="Proteom",main.cex=6,lwd=4,
             cex=3, fontface=5, cat.cex=4,sub.cex = 0.75,cat.pos=c(350,10,335,25),main.pos=c(0.5,0.8),
             disable.logging = TRUE,width=6000, height=4000,cat.dist=c(0.35,0.35,0.3,0.3),margin = 0.6)
setwd(source_wd)