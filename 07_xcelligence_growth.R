setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/xCelligence/Exp1_Dose response/processed")
#####Extracting H23#####
xCell_Data_H23<- read.csv("extracted_H23_rep1.csv")
xCell_Data_H23 <- data.frame(xCell_Data_H23[-c(1:37),])
xCell_Data_H23 <- data.frame(xCell_Data_H23[,2:8])
xCell_name <- c("time","2Gy", "1.5Gy", "1Gy", "05Gy", "0Gy", "Medium")
colnames(xCell_Data_H23) <- xCell_name
zeit<-stringr::str_pad(xCell_Data_H23$time, 9, side = "left", pad = 0)

h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_H23 <-(h+min+sec)/3600
#####Extracting HCC44####
xCell_Data_HCC44<- read.csv("extracted_HCC44_rep1.csv")
xCell_Data_HCC44 <- data.frame(xCell_Data_HCC44[-c(1:37),])
xCell_Data_HCC44 <- data.frame(xCell_Data_HCC44[,2:8])
xCell_name <- c("time","2Gy", "1.5Gy", "1Gy", "05Gy", "0Gy", "Medium")
colnames(xCell_Data_HCC44) <- xCell_name
zeit<-stringr::str_pad(xCell_Data_HCC44$time, 9, side = "left", pad = 0)

h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_HCC44 <-(h+min+sec)/3600
#####Graph####
setwd(save_wd)
jpeg(filename = "xCelligence_doseResponse.jpeg", width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(2,2,1,1,3,3),ncol=6,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt <- c("Control","0.5Gy", "1Gy","1.5Gy","2Gy", "Medium")
plot(xCell_inmin_H23,xCell_Data_H23$"0Gy", pch=16, cex=0.4,cex.main=4, axes=FALSE,mgp=c(7,1.5,0),xlim=c(0,175),  ylim = c(-0.5,4),cex.axis=3,cex.lab=4, col= "#5E4FA2", ylab="", xlab="")
points(xCell_inmin_H23,xCell_Data_H23$"05Gy", pch=16, col=  "#3288BD", cex=0.4)
points(xCell_inmin_H23,xCell_Data_H23$"1Gy", pch=16, col= "#ABDDA4", cex=0.4)
points(xCell_inmin_H23,xCell_Data_H23$"1.5Gy", pch=16, col="#FDAE61", cex=0.4)
points(xCell_inmin_H23,xCell_Data_H23$"2Gy", pch=16, col="#D53E4F", cex=0.4)
points(xCell_inmin_H23,xCell_Data_H23$Medium, pch=16, col="grey60", cex=0.4)
box(which="plot",lty="solid",lwd=3)
title(main="NCI-H23", cex.lab=3,cex.main=3,line=2)
title( ylab="Cell Index [a.u.]", cex.lab=3.5,cex.main=6,line=7,adj=0.75)
title( xlab="Time [h]", cex.lab=3.5,cex.main=6,line=7)
axis(1,lwd=3,at = c(0,100,175),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 3.5,0))
axis(2,lwd=3,at = c(0,2,4),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 3,0),las=2)
#legend(c(0,180),c(4,2.6),leg.txt,inset=c(-0.6,-0.85),bty = "n", bg = NA,box.col="white",box.lwd = 2,cex=2.5, ncol=2,pch=16,adj=c(0,0), col=c("#5E4FA2","#3288BD","#ABDDA4","#FDAE61","#D53E4F","grey60"))
minor.tick(nx =2/1,ny=2/1,tick.ratio =1.5,  y.args = list( lwd = 3,tcl = -1),  x.args = list( lwd = 3,tcl = -1)) 

plot(xCell_inmin_HCC44,xCell_Data_HCC44$"0Gy", pch=16, cex=0.4,cex.main=4, axes=FALSE,mgp=c(7,1.5,0),xlim=c(0,190),  ylim = c(-0.5,4),cex.axis=3,cex.lab=4, col="#5E4FA2", ylab="", xlab="")
points(xCell_inmin_HCC44,xCell_Data_HCC44$"05Gy", pch=16, col="#3288BD", cex=0.4)
points(xCell_inmin_HCC44,xCell_Data_HCC44$"1Gy", pch=16, col="#ABDDA4", cex=0.4)
points(xCell_inmin_HCC44,xCell_Data_HCC44$"1.5Gy", pch=16, col="#FDAE61", cex=0.4)
points(xCell_inmin_HCC44,xCell_Data_HCC44$"2Gy", pch=16, col="#D53E4F", cex=0.4)
points(xCell_inmin_HCC44,xCell_Data_HCC44$Medium, pch=16, col="grey60", cex=0.4)
box(which="plot",lty="solid",lwd=3)
title(main="HCC44", cex.lab=3.5,cex.main=3,line=2)
title( ylab="Cell Index [a.u.]", cex.lab=3.5,cex.main=6,line=7,adj=0.75)
title( xlab="Time [h]", cex.lab=3.5,cex.main=6,line=7)
axis(1,lwd=3,at = c(0,100,175),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 3.5,0))
axis(2,lwd=3,at = c(0,2,4),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 3,0),las=2)
minor.tick(nx =2/1,ny=2/1,tick.ratio =1.5,  y.args = list( lwd = 3,tcl = -1),  x.args = list( lwd = 3,tcl = -1)) 

par(mar = c(10, 0, 6, 4),xpd=TRUE)
plot(1:10,1:10, col="white",ann=FALSE, axes=FALSE)
legend(x=c(1,12),y=c(1,11),leg.txt,box.lwd = 3,cex=3.5,bty = "n", bg = NA,pt.cex = 4, ncol=1,pch=16, col=c("#5E4FA2", "#3288BD","#ABDDA4","#FDAE61","#D53E4F","grey60"),y.intersp = 0.95,x.intersp=0.5)
dev.off()
setwd(source_wd)