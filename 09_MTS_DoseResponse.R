####importing Data####
setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/MTS assays/Exp1_DoseResponse/processed")
HCC44_1<- read.csv("DoseResponse_mean_HCC44_Rep1.csv")
HCC44_2<- read.csv("DoseResponse_mean_HCC44_Rep2.csv")
HCC44_3<- read.csv("DoseResponse_mean_HCC44_Rep3.csv")
H23_1<- read.csv("DoseResponse_mean_H23_Rep1.csv")
H23_2<- read.csv("DoseResponse_mean_H23_Rep2.csv")
H23_3<- read.csv("DoseResponse_mean_H23_Rep3.csv")

raw_HCC44_1<- read.csv("DoseResponse_HCC44_Rep1.csv")
raw_HCC44_2<- read.csv("DoseResponse_HCC44_Rep2.csv")
raw_HCC44_3<- read.csv("DoseResponse_HCC44_Rep3.csv")
raw_H23_1<- read.csv("DoseResponse_H23_Rep1.csv")
raw_H23_2<- read.csv("DoseResponse_H23_Rep2.csv")
raw_H23_3<- read.csv("DoseResponse_H23_Rep3.csv")
time <- as.vector(c(0, 24, 48, 72, 96, 120))
####HCC44 control - extraction####
HCC44_c <- data.frame(HCC44_1[,c(2,3:5)],HCC44_2[,3:5],HCC44_3[,3:5])
HCC44_c["mean_HCC44_c_all"]<- rowMeans(HCC44_c[,2:10])
HCC44_c["SD_HCC44_c_all"]<- apply(HCC44_c[,2:10], 1, sd)
HCC44_c["SEM_HCC44_c_all"]<- HCC44_c["SD_HCC44_c_all"]/sqrt(3)
####HCC44 2Gy - extraction####
HCC44_2Gy <- data.frame(HCC44_1[,c(2,6:8)],HCC44_2[,6:8],HCC44_3[,6:8])
HCC44_2Gy["mean_HCC44_2Gy_all"]<- rowMeans(HCC44_2Gy[,2:10])
HCC44_2Gy["SD_HCC44_2Gy_all"]<- apply(HCC44_2Gy[,2:10], 1, sd)
HCC44_2Gy["SEM_HCC44_2Gy_all"]<- HCC44_2Gy["SD_HCC44_2Gy_all"]/sqrt(3)
####H23 control - extraction ####
H23_c <- data.frame(H23_1[,c(2,3:5)],H23_2[,3:5],H23_3[,3:5])
H23_c["mean_H23_c_all"]<- rowMeans(H23_c[,2:4])
H23_c["SD_H23_c_all"]<- apply(H23_c[,2:4], 1, sd)
H23_c["SEM_H23_c_all"]<- H23_c["SD_H23_c_all"]/sqrt(3)
####H23 2Gy - extraction####
H23_2Gy <- data.frame(H23_1[,c(2,6:8)],H23_2[,6:8],H23_3[,6:8])
H23_2Gy["mean_H23_2Gy_all"]<- rowMeans(H23_2Gy[,2:10])
H23_2Gy["SD_H23_2Gy_all"]<- apply(H23_2Gy[,2:10], 1, sd)
H23_2Gy["SEM_H23_2Gy_all"]<- H23_2Gy["SD_H23_2Gy_all"]/sqrt(3)
####HCC44 - binding for ANOVA####
AOV_HCC44_c <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_HCC44_1[1:10,2:7]))),range01(as.numeric(unlist(raw_HCC44_1[13:22,2:7]))),range01(as.numeric(unlist(raw_HCC44_1[25:34,2:7]))),
                                      range01(as.numeric(unlist(raw_HCC44_2[1:10,2:7]))),range01(as.numeric(unlist(raw_HCC44_2[13:22,2:7]))),range01(as.numeric(unlist(raw_HCC44_2[25:34,2:7]))),
                                      range01(as.numeric(unlist(raw_HCC44_3[1:10,2:7]))),range01(as.numeric(unlist(raw_HCC44_3[13:22,2:7]))),range01(as.numeric(unlist(raw_HCC44_3[25:34,2:7])))),
                          Time=rep(c(rep("0h",10),rep("24h",10),rep("48h",10),rep("72h",10),rep("96h",10),rep("120h",10)),9),
                          Treatment =rep("control",540),
                          BioReplicate=c(rep("Rep1",180),rep("Rep2",180),rep("Rep3",180)),
                          Plate=rep(c(rep("Plate1",60),rep("Plate2",60),rep("Plate3",60)),3),
                          Cell_line=c(rep("HCC44",540)))
AOV_HCC44_2Gy <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_HCC44_1[1:10,10:15]))),range01(as.numeric(unlist(raw_HCC44_1[13:22,10:15]))),range01(as.numeric(unlist(raw_HCC44_1[25:34,10:15]))),
                                      range01(as.numeric(unlist(raw_HCC44_2[1:10,10:15]))),range01(as.numeric(unlist(raw_HCC44_2[13:22,10:15]))),range01(as.numeric(unlist(raw_HCC44_2[25:34,10:15]))),
                                      range01(as.numeric(unlist(raw_HCC44_3[1:10,10:15]))),range01(as.numeric(unlist(raw_HCC44_3[13:22,10:15]))),range01(as.numeric(unlist(raw_HCC44_3[25:34,10:15])))),
                          Time=rep(c(rep("0h",10),rep("24h",10),rep("48h",10),rep("72h",10),rep("96h",10),rep("120h",10)),9),
                          Treatment =rep("2Gy",540),
                          BioReplicate=c(rep("Rep1",180),rep("Rep2",180),rep("Rep3",180)),
                          Plate=rep(c(rep("Plate1",60),rep("Plate2",60),rep("Plate3",60)),3),
                          Cell_line=rep("HCC44",540))
AOV_HCC44 <- rbind(AOV_HCC44_c,AOV_HCC44_2Gy)
#model_HCC44 <- lmer(Intensity ~ Treatment + (1|BioReplicate/Plate), data = AOV_HCC44)
#summary(model_HCC44)
####H23 - binding for ANOVA####
AOV_H23_c <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_H23_1[1:10,2:7]))),range01(as.numeric(unlist(raw_H23_1[13:22,2:7]))),range01(as.numeric(unlist(raw_H23_1[25:34,2:7]))),
                                    range01(as.numeric(unlist(raw_H23_2[1:10,2:7]))),range01(as.numeric(unlist(raw_H23_2[13:22,2:7]))),range01(as.numeric(unlist(raw_H23_2[25:34,2:7]))),
                                    range01(as.numeric(unlist(raw_H23_3[1:10,2:7]))),range01(as.numeric(unlist(raw_H23_3[13:22,2:7]))),range01(as.numeric(unlist(raw_H23_3[25:34,2:7])))),
                          Time=rep(c(rep("0h",10),rep("24h",10),rep("48h",10),rep("72h",10),rep("96h",10),rep("120h",10)),9),
                          Treatment =rep("control",540),
                          BioReplicate=c(rep("Rep1",180),rep("Rep2",180),rep("Rep3",180)),
                          Plate=rep(c(rep("Plate1",60),rep("Plate2",60),rep("Plate3",60)),3),
                          Cell_line=rep("H23",540))
AOV_H23_2Gy <- data.frame(Intensity=c(range01(as.numeric(unlist(raw_H23_1[1:10,10:15]))),range01(as.numeric(unlist(raw_H23_1[13:22,10:15]))),range01(as.numeric(unlist(raw_H23_1[25:34,10:15]))),
                                      range01( as.numeric(unlist(raw_H23_2[1:10,10:15]))),range01(as.numeric(unlist(raw_H23_2[13:22,10:15]))),range01(as.numeric(unlist(raw_H23_2[25:34,10:15]))),
                                        range01(as.numeric(unlist(raw_H23_3[1:10,10:15]))),range01(as.numeric(unlist(raw_H23_3[13:22,10:15]))),range01(as.numeric(unlist(raw_H23_3[25:34,10:15])))),
                            Time=rep(c(rep("0h",10),rep("24h",10),rep("48h",10),rep("72h",10),rep("96h",10),rep("120h",10)),9),
                            Treatment =rep("2Gy",540),
                            BioReplicate=c(rep("Rep1",180),rep("Rep2",180),rep("Rep3",180)),
                            Plate=rep(c(rep("Plate1",60),rep("Plate2",60),rep("Plate3",60)),3),
                            Cell_line=rep("H23",540))
AOV_H23 <- rbind(AOV_H23_c,AOV_H23_2Gy)
#model_H23 <- lmer(Intensity ~ Treatment + (1|BioReplicate/Plate), data = AOV_H23)
#summary(model_H23)
####Plotting H23 and HCC44####
plot_all <- cbind(HCC44_1[,3:5],HCC44_2[,3:5],HCC44_3[,3:5],HCC44_1[,6:8],HCC44_2[,6:8],HCC44_3[,6:8],H23_1[,3:5],H23_2[,3:5],H23_3[,3:5],H23_1[,6:8],H23_2[,6:8],H23_3[,6:8])
plot_all<-range01(plot_all[,1:36])
setwd(save_wd)
jpeg(filename = "MTS_doseResponse.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1),ncol=3,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt<- c("HCC44 c","HCC44 2 Gy","NCI-H23 c","NCI-H23 2 Gy")
plot(time,rowMeans(plot_all[,19:27]),pch=16,col="#66C2A5", type="l", mgp=c(7,1.5, 0),cex.lab=4,axes=FALSE,lwd=3, ylab="", xlab="", ylim=c(0,1.25))
lines(time,rowMeans(plot_all[,28:36]),pch=16,col="#FDAE61",lwd=3)
lines(time,rowMeans(plot_all[,1:9]),pch=16,col="#3288BD",lwd=3)
lines(time,rowMeans(plot_all[,10:18]),pch=16,col="#D53E4F",lwd=3)
arrows(time, y0=rowMeans(plot_all[,1:9])-apply(plot_all[,1:9],1 , sd), x1=time, y1=rowMeans(plot_all[,1:9])+apply(plot_all[,1:9],1 , sd), code=3, angle=90, length=0.1,lwd=3,col="grey50")
arrows(time, y0=rowMeans(plot_all[,10:18])-apply(plot_all[,10:18],1 , sd), x1=time, y1=rowMeans(plot_all[,10:18])+apply(plot_all[,10:18],1 , sd), code=3, angle=90, length=0.1,lwd=3)
arrows(time, y0=rowMeans(plot_all[,19:27])-apply(plot_all[,19:27],1 , sd), x1=time, y1=rowMeans(plot_all[,19:27])+apply(plot_all[,19:27],1 , sd), code=3, angle=90, length=0.1,lwd=3,col="grey50")
arrows(time, y0=rowMeans(plot_all[,28:36])-apply(plot_all[,28:36],1 , sd), x1=time, y1=rowMeans(plot_all[,28:36])+apply(plot_all[,28:36],1 , sd), code=3, angle=90, length=0.1,lwd=3)
axis(1,lwd=2,at = c(0,60,120),labels=c(0,60,120),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5, mgp=c(5,3.5, 0))
axis(2,lwd=2,at = c(0,0.5,1),labels=c(0,50,100),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5, mgp=c(5,2.5, 0),las=2)
legend(x=c(0,120),y=c(1.25,0.9), box.col="white",leg.txt,ncol=2,bty = "n", bg = NA,adj=c(0,0.25),lwd=3,pt.cex=3,pt.lwd=3, cex = 2.5,col=c("#3288BD","#D53E4F","#66C2A5","#FDAE61"),box.lwd=2,y.intersp = 0.9)
box(lty="solid",lwd=3)
title( ylab="Viability [%]", cex.lab=3.5,cex.main=6,line=9)
title( xlab="Time[h]", cex.lab=3.5,cex.main=6,line=7)
minor.tick(nx =2/1,ny=4/2,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
dev.off()
setwd(source_wd)
