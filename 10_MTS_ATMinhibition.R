#####set directory####
setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/MTS assays/Exp3_ATMinhibition/processed")
time <- as.vector(c(0, 24, 48, 72, 96, 120))
#####loading data####
mean_H23_DMSO_2Gy<- read.csv("ATMinhibition_mean_H23_DMSO_2Gy.csv")
mean_H23_DMSO_c<- read.csv("ATMinhibition_mean_H23_DMSO_c.csv")
mean_H23_AZD1390_2Gy<- read.csv("ATMinhibition_mean_H23_AZD1390_2Gy.csv")
mean_H23_AZD1390_c<- read.csv("ATMinhibition_mean_H23_AZD1390_c.csv")

raw_H23_1 <- read.csv("ATMinhibition_H23_rep1.csv")
raw_H23_1 <- data.frame(lapply(raw_H23_1, as.numeric))
raw_H23_2 <- read.csv("ATMinhibition_H23_rep2.csv")
raw_H23_2 <- data.frame(lapply(raw_H23_2, as.numeric))
raw_H23_3 <- read.csv("ATMinhibition_H23_rep3.csv")
raw_H23_3 <- data.frame(lapply(raw_H23_3, as.numeric))
raw_H23_4 <- read.csv("ATMinhibition_H23_rep4.csv")
raw_H23_4 <- data.frame(lapply(raw_H23_4, as.numeric))

mean_HCC44_DMSO_2Gy<- read.csv("ATMinhibition_mean_HCC44_DMSO_2Gy.csv")
mean_HCC44_DMSO_c<- read.csv("ATMinhibition_mean_HCC44_DMSO_c.csv")
mean_HCC44_AZD1390_2Gy<- read.csv("ATMinhibition_mean_HCC44_AZD1390_2Gy.csv")
mean_HCC44_AZD1390_c<- read.csv("ATMinhibition_mean_HCC44_AZD1390_c.csv")

raw_HCC44_1 <- read.csv("ATMinhibition_HCC44_rep1.csv")
raw_HCC44_1 <- data.frame(lapply(raw_HCC44_1, as.numeric))
raw_HCC44_2 <- read.csv("ATMinhibition_HCC44_rep2.csv")
raw_HCC44_2 <- data.frame(lapply(raw_HCC44_2, as.numeric))
raw_HCC44_3 <- read.csv("ATMinhibition_HCC44_rep3.csv")
raw_HCC44_3 <- data.frame(lapply(raw_HCC44_3, as.numeric))
raw_HCC44_4 <- read.csv("ATMinhibition_HCC44_rep4.csv")
raw_HCC44_4 <- data.frame(lapply(raw_HCC44_4, as.numeric))
####HCC44 - binding for ANOVA####
AOV_HCC44_DMSO_c <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_HCC44_1[11:15,3:8]))),range01(as.numeric(unlist(raw_HCC44_2[11:15,3:8]))),range01(as.numeric(unlist(raw_HCC44_3[11:15,3:8]))),range01(as.numeric(unlist(raw_HCC44_4[11:15,3:8])))),
                                      Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                                      Treatment =rep("DMSO control",120),
                                      BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                                      Cell_line=rep("HCC44",120))
AOV_HCC44_DMSO_2Gy <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_HCC44_1[3:7,3:8]))),range01(as.numeric(unlist(raw_HCC44_2[3:7,3:8]))),range01(as.numeric(unlist(raw_HCC44_3[3:7,3:8]))),range01(as.numeric(unlist(raw_HCC44_4[3:7,3:8])))),
                               Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                               Treatment =rep("DMSO 2Gy",120),
                               BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                               Cell_line=rep("HCC44",120))
AOV_HCC44_AZD1390_c <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_HCC44_1[11:15,11:16]))),range01(as.numeric(unlist(raw_HCC44_2[11:15,11:16]))),range01(as.numeric(unlist(raw_HCC44_3[11:15,11:16]))),range01(as.numeric(unlist(raw_HCC44_4[11:15,11:16])))),
                               Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                               Treatment =rep("DMSO control",120),
                               BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                               Cell_line=rep("HCC44",120))
AOV_HCC44_AZD1390_2Gy <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_HCC44_1[3:7,11:16]))),range01(as.numeric(unlist(raw_HCC44_2[3:7,11:16]))),range01(as.numeric(unlist(raw_HCC44_3[3:7,11:16]))),range01(as.numeric(unlist(raw_HCC44_4[3:7,11:16])))),
                                 Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                                 Treatment =rep("DMSO 2Gy",120),
                                 BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                                 Cell_line=rep("HCC44",120))
####H23 - binding for ANOVA####
AOV_H23_DMSO_c <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_H23_1[11:15,3:8]))),range01(as.numeric(unlist(raw_H23_2[11:15,3:8]))),range01(as.numeric(unlist(raw_H23_3[11:15,3:8]))),range01(as.numeric(unlist(raw_H23_4[11:15,3:8])))),
                               Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                               Treatment =rep("DMSO control",120),
                               BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                               Cell_line=rep("H23",120))
AOV_H23_DMSO_2Gy <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_H23_1[3:7,3:8]))),range01(as.numeric(unlist(raw_H23_2[3:7,3:8]))),range01(as.numeric(unlist(raw_H23_3[3:7,3:8]))),range01(as.numeric(unlist(raw_H23_4[3:7,3:8])))),
                                 Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                                 Treatment =rep("DMSO 2Gy",120),
                                 BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                                 Cell_line=rep("H23",120))
AOV_H23_AZD1390_c <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_H23_1[11:15,11:16]))),range01(as.numeric(unlist(raw_H23_2[11:15,11:16]))),range01(as.numeric(unlist(raw_H23_3[11:15,11:16]))),range01(as.numeric(unlist(raw_H23_4[11:15,11:16])))),
                                  Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                                  Treatment =rep("DMSO control",120),
                                  BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                                  Cell_line=rep("H23",120))
AOV_H23_AZD1390_2Gy <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_H23_1[3:7,11:16]))),range01(as.numeric(unlist(raw_H23_2[3:7,11:16]))),range01(as.numeric(unlist(raw_H23_3[3:7,11:16]))),range01(as.numeric(unlist(raw_H23_4[3:7,11:16])))),
                                    Time=rep(c(rep("0h",5),rep("24h",5),rep("48h",5),rep("72h",5),rep("96h",5),rep("120h",5)),4),
                                    Treatment =rep("DMSO 2Gy",120),
                                    BioReplicate=c(rep("Rep1",30),rep("Rep2",30),rep("Rep3",30),rep("Rep4",30)),
                                    Cell_line=rep("H23",120))
#####ANOVA####
HCC44_Rad <- rbind(AOV_HCC44_DMSO_2Gy,AOV_HCC44_AZD1390_2Gy)
#model_HCC44_2Gy <- lmer(Intensity ~ Treatment + (1|BioReplicate), data = HCC44_Rad)
HCC44_con <- rbind(AOV_HCC44_DMSO_c,AOV_HCC44_AZD1390_c)
#model_HCC44_c <- lmer(Intensity ~ Treatment + (1|BioReplicate), data = HCC44_Rad)

H23_Rad <- rbind(AOV_H23_DMSO_2Gy,AOV_H23_AZD1390_2Gy)
#model_H23_2Gy <- lmer(Intensity ~ Treatment + (1|BioReplicate), data = H23_Rad)
H23_con <- rbind(AOV_H23_DMSO_c,AOV_H23_AZD1390_c)
#model_H23_c <- lmer(Intensity ~ Treatment + (1|BioReplicate), data = H23_Rad)
#####HCC44####
plot_HCC44 <- data.frame(mean_HCC44_DMSO_c[,c(3,5,7,9)],mean_HCC44_DMSO_2Gy[,c(3,5,7,9)],mean_HCC44_AZD1390_c[,c(3,5,7,9)],mean_HCC44_AZD1390_2Gy[,c(3,5,7,9)])
plot_HCC44 <- range01(plot_HCC44)
#####H23####
plot_H23 <- data.frame(mean_H23_DMSO_c[,c(3,5,7,9)],mean_H23_DMSO_2Gy[,c(3,5,7,9)],mean_H23_AZD1390_c[,c(3,5,7,9)],mean_H23_AZD1390_2Gy[,c(3,5,7,9)])
plot_H23 <- range01(plot_H23)
####Graphes####
setwd(save_wd)
jpeg(filename = "MTS_AZD1390_DoseResponse.jpeg",  width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1,2,2,2),ncol=6,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt <- c("control DMSO", "control AZD1390","2Gy DMSO", "2Gy AZD1390")
plot(time,rowMeans(plot_HCC44[,1:4]),pch=16,col="#66C2A5", type="l", mgp=c(7,1.5, 0),cex.lab=4,axes=FALSE,lwd=3, ylab="", xlab="", ylim=c(0,1.25))
lines(time,rowMeans(plot_HCC44[,5:8]),pch=16,col="#FDAE61",lwd=3)
lines(time,rowMeans(plot_HCC44[,9:12]),pch=16,col="#3288BD",lwd=3)
lines(time,rowMeans(plot_HCC44[,13:16]),pch=16,col="#D53E4F",lwd=3)
arrows(x0=time, y0=rowMeans(plot_HCC44[,1:4])-apply(plot_HCC44[,1:4],1 , sd), x1=time, y1=rowMeans(plot_HCC44[,1:4])+apply(plot_HCC44[,1:4],1 , sd), code=3, angle=90, length=0.1,lwd=3,col="grey50")
arrows(x0=time, y0=rowMeans(plot_HCC44[,5:8])-apply(plot_HCC44[,5:8],1 , sd), x1=time, y1=rowMeans(plot_HCC44[,5:8])+apply(plot_HCC44[,5:8],1 , sd), code=3, angle=90, length=0.1,lwd=3)
arrows(x0=time, y0=rowMeans(plot_HCC44[,9:12])-apply(plot_HCC44[,9:12],1 , sd), x1=time, y1=rowMeans(plot_HCC44[,9:12])+apply(plot_HCC44[,9:12],1 , sd), code=3, angle=90, length=0.1,lwd=3,col="grey50")
arrows(x0=time, y0=rowMeans(plot_HCC44[,13:16])-apply(plot_HCC44[,13:16],1 , sd), x1=time, y1=rowMeans(plot_HCC44[,13:16])+apply(plot_HCC44[,13:16],1 , sd), code=3, angle=90, length=0.1,lwd=3)
axis(1,lwd=2,at = c(0,60,120),labels=c(0,60,120),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5, mgp=c(5,3.5, 0))
axis(2,lwd=2,at = c(0,0.5,1),labels=c(0,50,100),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5, mgp=c(5,2.5, 0),las=2)
legend(x=c(-2,122),y=c(1.2,1.1),leg.txt,bty = "n", bg = NA,ncol=1,adj=c(0.05,0.1),lwd=3,pt.cex=3,pt.lwd=3, cex = 2.75,col=c("#ABDDA4","#3288BD","#FDAE61","#D53E4F"),box.lwd=2,y.intersp = 0.9,x.intersp =1)
box(lty="solid",lwd=3)
title( ylab="Viability [%]", cex.lab=3.5,cex.main=6,line=9)
title( xlab="Time[h]", cex.lab=3.5,cex.main=6,line=7)
minor.tick(nx =2/1,ny=4/2,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 

plot(time,rowMeans(plot_H23[,1:4]),pch=16,col="#66C2A5", type="l", mgp=c(7,1.5, 0),cex.lab=4,axes=FALSE,lwd=3, ylab="", xlab="", ylim=c(0,1.25))
lines(time,rowMeans(plot_H23[,5:8]),pch=16,col="#FDAE61",lwd=3)
lines(time,rowMeans(plot_H23[,9:12]),pch=16,col="#3288BD",lwd=3)
lines(time,rowMeans(plot_H23[,13:16]),pch=16,col="#D53E4F",lwd=3)
arrows(x0=time, y0=rowMeans(plot_H23[,1:4])-apply(plot_H23[,1:4],1 , sd), x1=time, y1=rowMeans(plot_H23[,1:4])+apply(plot_H23[,1:4],1 , sd), code=3, angle=90, length=0.1,lwd=3,col="grey50")
arrows(x0=time, y0=rowMeans(plot_H23[,5:8])-apply(plot_H23[,5:8],1 , sd), x1=time, y1=rowMeans(plot_H23[,5:8])+apply(plot_H23[,5:8],1 , sd), code=3, angle=90, length=0.1,lwd=3)
arrows(x0=time, y0=rowMeans(plot_H23[,9:12])-apply(plot_H23[,9:12],1 , sd), x1=time, y1=rowMeans(plot_H23[,9:12])+apply(plot_H23[,9:12],1 , sd), code=3, angle=90, length=0.1,lwd=3,col="grey50")
arrows(x0=time, y0=rowMeans(plot_H23[,13:16])-apply(plot_H23[,13:16],1 , sd), x1=time, y1=rowMeans(plot_H23[,13:16])+apply(plot_H23[,13:16],1 , sd), code=3, angle=90, length=0.1,lwd=3)
axis(1,lwd=2,at = c(0,60,120),labels=c(0,60,120),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5, mgp=c(5,3.5, 0))
axis(2,lwd=2,at = c(0,0.5,1),labels=c(0,50,100),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5, mgp=c(5,2.5, 0),las=2)
legend(x=c(-2,122),y=c(1.2,1.1),leg.txt,bty = "n", bg = NA,ncol=1,adj=c(0.05,0.1),lwd=3,pt.cex=3,pt.lwd=3, cex = 2.75,col=c("#ABDDA4","#3288BD","#FDAE61","#D53E4F"),box.lwd=2,y.intersp = 0.9,x.intersp =1)
box(lty="solid",lwd=3)
title( ylab="Viability [%]", cex.lab=3.5,cex.main=6,line=9)
title( xlab="Time[h]", cex.lab=3.5,cex.main=6,line=7)
minor.tick(nx =2/1,ny=4/2,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
dev.off()
setwd(source_wd)