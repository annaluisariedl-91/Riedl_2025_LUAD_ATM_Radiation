setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/xCelligence/Exp2_ATMinhibition/processed")
#####Extracting H23#####
xCell_name <- c("time","2Gy_AZD1390", "2Gy_DMSO", "05Gy_AZD1390", "05Gy_DMSO","0Gy_AZD1390", "0Gy_DMSO", "Medium")
xCell_Data_H23_1<- read.csv("extracted_mean_H23_rep1.csv",sep="\t")
xCell_Data_H23_1 <- data.frame(xCell_Data_H23_1[-(1:15),c(2,9:3)])
colnames(xCell_Data_H23_1) <- xCell_name
xCell_Data_H23_2<- read.csv("extracted_mean_H23_rep2.csv",sep="\t")
xCell_Data_H23_2 <- data.frame(xCell_Data_H23_2[-(1:15),2:9])
colnames(xCell_Data_H23_2) <- xCell_name
xCell_Data_H23_3<- read.csv("extracted_mean_H23_rep3.csv",sep="\t")
xCell_Data_H23_3 <- data.frame(xCell_Data_H23_3[-(1:15),c(2,4:9,3)])
colnames(xCell_Data_H23_3) <- xCell_name
xCell_Data_H23_4<- read.csv("extracted_mean_H23_rep4.csv",sep="\t")
xCell_Data_H23_4 <- data.frame(xCell_Data_H23_4[-(1:15),2:9])
colnames(xCell_Data_H23_4) <- xCell_name

zeit<-stringr::str_pad(xCell_Data_H23_1$time, 9, side = "left", pad = 0)
h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_H23_1 <-(h+min+sec)/3600

zeit<-stringr::str_pad(xCell_Data_H23_2$time, 9, side = "left", pad = 0)
h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_H23_2 <-(h+min+sec)/3600

zeit<-stringr::str_pad(xCell_Data_H23_3$time, 9, side = "left", pad = 0)
h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_H23_3 <-(h+min+sec)/3600

zeit<-stringr::str_pad(xCell_Data_H23_4$time, 9, side = "left", pad = 0)
h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_H23_4 <-(h+min+sec)/3600
#####
mean_data_H23<-data.frame("time[h]"=xCell_Data_H23_1[,2])
mean_data_H23$mean_H23_DMSO_2Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,3])),data.frame(as.numeric(unlist(xCell_Data_H23_2[,3]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,3]))))),na.rm=FALSE)
mean_data_H23$sem_H23_DMSO_2Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,3]))),data.frame(as.numeric(unlist(xCell_Data_H23_2[,3]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,3])))), 1, function(x) sd(x))/sqrt(3)
mean_data_H23$mean_H23_10nM_AZD1390_2Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,2])),data.frame(as.numeric(unlist(xCell_Data_H23_2[,2]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,2]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,2]))))),na.rm=FALSE)
mean_data_H23$sem_H23_10nM_AZD1390_2Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,2]))),data.frame(as.numeric(unlist(xCell_Data_H23_2[,2]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,2]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,2])))), 1, function(x) sd(x))/sqrt(4)
mean_data_H23$mean_H23_DMSO_05Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,5])),data.frame(as.numeric(unlist(xCell_Data_H23_2[,5]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,5]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,5]))))),na.rm=FALSE)
mean_data_H23$sem_H23_DMSO_05Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,5]))),data.frame(as.numeric(unlist(xCell_Data_H23_2[,5]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,5]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,5])))), 1, function(x) sd(x))/sqrt(4)
mean_data_H23$mean_H23_10nM_AZD1390_05Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,4])),data.frame(as.numeric(unlist(xCell_Data_H23_2[,4]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,4]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,4]))))),na.rm=FALSE)
mean_data_H23$sem_H23_10nM_AZD1390_05Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,4]))),data.frame(as.numeric(unlist(xCell_Data_H23_2[,4]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,4]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,4])))), 1, function(x) sd(x))/sqrt(4)
mean_data_H23$mean_H23_DMSO_0Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,7])),data.frame(as.numeric(unlist(xCell_Data_H23_2[,7]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,7]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,7]))))),na.rm=FALSE)
mean_data_H23$sem_H23_DMSO_0Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,7]))),data.frame(as.numeric(unlist(xCell_Data_H23_2[,7]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,7]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,7])))), 1, function(x) sd(x))/sqrt(4)
mean_data_H23$mean_H23_10nM_AZD1390_0Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,6])),data.frame(as.numeric(unlist(xCell_Data_H23_2[,6]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,6]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,6]))))),na.rm=FALSE)
mean_data_H23$sem_H23_10nM_AZD1390_0Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,6]))),data.frame(as.numeric(unlist(xCell_Data_H23_2[,6]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,6]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,6])))), 1, function(x) sd(x))/sqrt(4)
mean_data_H23$mean_H23_med<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,8])),data.frame(as.numeric(unlist(xCell_Data_H23_2[,8]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,8]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,8]))))),na.rm=FALSE)
mean_data_H23$sem_H23_med<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_H23_1[,8]))),data.frame(as.numeric(unlist(xCell_Data_H23_2[,8]))),data.frame(as.numeric(unlist(xCell_Data_H23_3[,8]))),data.frame(as.numeric(unlist(xCell_Data_H23_4[,8])))), 1, function(x) sd(x))/sqrt(4)
#####Extracting HCC44####
xCell_name <- c("time","4Gy_AZD1390", "4Gy_DMSO", "2Gy_AZD1390", "2Gy_DMSO","0Gy_AZD1390", "0Gy_DMSO", "Medium")
xCell_Data_HCC44_1<- read.csv("extracted_mean_HCC44_rep1.csv", sep="\t")
xCell_Data_HCC44_1 <- data.frame(xCell_Data_HCC44_1[,c(2,9:3)])
xCell_Data_HCC44_1 <- data.frame(xCell_Data_HCC44_1[-c(1:15),])
colnames(xCell_Data_HCC44_1) <- xCell_name
xCell_Data_HCC44_2<- read.csv("extracted_mean_HCC44_rep2.csv",sep="\t")
xCell_Data_HCC44_2 <- data.frame(xCell_Data_HCC44_2[-c(1:15),])
xCell_Data_HCC44_2 <- data.frame(xCell_Data_HCC44_2[,2:9])
colnames(xCell_Data_HCC44_2) <- xCell_name
xCell_Data_HCC44_3<- read.csv("extracted_mean_HCC44_rep3.csv", sep="\t")
xCell_Data_HCC44_3 <- data.frame(xCell_Data_HCC44_3[-c(1:15),])
xCell_Data_HCC44_3 <- data.frame(xCell_Data_HCC44_3[,c(2:9)])
colnames(xCell_Data_HCC44_3) <- xCell_name

zeit<-stringr::str_pad(xCell_Data_HCC44_1$time, 9, side = "left", pad = 0)
h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_HCC44_1 <-(h+min+sec)/3600

zeit<-stringr::str_pad(xCell_Data_HCC44_2$time, 9, side = "left", pad = 0)
h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_HCC44_2 <-(h+min+sec)/3600

zeit<-stringr::str_pad(xCell_Data_HCC44_3$time, 9, side = "left", pad = 0)
h<-substr(zeit, 1,3)
h <- as.double(h)
h <- h*3600
min <-substr(zeit, 5,6)
min <-as.double(min)
min <- min*60
sec<-substr(zeit, 8,9)
sec <- as.double(sec)
xCell_inmin_HCC44_3 <-(h+min+sec)/3600
#####
mean_data_HCC44<-data.frame("time[h]"=xCell_Data_HCC44_1[,2])
mean_data_HCC44$mean_HCC44_DMSO_2Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,3])),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,3]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3]))))),na.rm=FALSE)
mean_data_HCC44$sem_HCC44_DMSO_2Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,3]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,3]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3])))), 1, function(x) sd(x))/sqrt(3)
mean_data_HCC44$mean_HCC44_10nM_AZD1390_2Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,2])),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,2]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3]))))),na.rm=FALSE)
mean_data_HCC44$sem_HCC44_10nM_AZD1390_2Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,2]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,2]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,2])))), 1, function(x) sd(x))/sqrt(3)
mean_data_HCC44$mean_HCC44_DMSO_05Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,5])),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,5]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3]))))),na.rm=FALSE)
mean_data_HCC44$sem_HCC44_DMSO_05Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,5]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,5]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,5])))), 1, function(x) sd(x))/sqrt(3)
mean_data_HCC44$mean_HCC44_10nM_AZD1390_05Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,4])),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,4]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3]))))),na.rm=FALSE)
mean_data_HCC44$sem_HCC44_10nM_AZD1390_05Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,4]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,4]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,4])))), 1, function(x) sd(x))/sqrt(3)
mean_data_HCC44$mean_HCC44_DMSO_0Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,7])),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,7]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3]))))),na.rm=FALSE)
mean_data_HCC44$sem_HCC44_DMSO_0Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,7]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,7]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,7])))), 1, function(x) sd(x))/sqrt(3)
mean_data_HCC44$mean_HCC44_10nM_AZD1390_0Gy<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,6])),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,6]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3]))))),na.rm=FALSE)
mean_data_HCC44$sem_HCC44_10nM_AZD1390_0Gy<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,6]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,6]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,6])))), 1, function(x) sd(x))/sqrt(3)
mean_data_HCC44$mean_HCC44_med<- rowMeans(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,8])),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,8]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,3]))))),na.rm=FALSE)
mean_data_HCC44$sem_HCC44_med<- apply(cbind(data.frame(as.numeric(unlist(xCell_Data_HCC44_1[,8]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_2[,8]))),data.frame(as.numeric(unlist(xCell_Data_HCC44_3[,8])))), 1, function(x) sd(x))/sqrt(3)

#####Graphes####
setwd(save_wd)
jpeg(filename = "xCelligence_DoseResponse_AZD1390.jpeg", width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1,2,2,2),ncol=6,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 10),xpd=TRUE)
leg.txt <- c("HCC44 DMSO","HCC44 AZD1390", "NCI-H23 DMSO","NCI-H23 AZD1390","Medium")
plot(xCell_inmin_H23_1,mean_data_H23[,10], pch=16, cex=0.5,lwd.tick=3,cex.main=4,  ylim = c(-0.5,6), axes=FALSE, col="#FDAE61", ann=FALSE,cex.lab=4,cex.axis=3,mgp=c(7,1,0))
points(xCell_inmin_H23_1,mean_data_H23[,12], pch=16, col="#66C2A5", cex=0.5)
points(xCell_inmin_H23_1,mean_data_H23[,13], pch=16, col="grey", cex=0.2)
points(xCell_inmin_H23_1,mean_data_HCC44[,12], pch=16, cex=0.5, col="#D53E4F")
points(xCell_inmin_H23_1,mean_data_HCC44[,10], pch=16, col="#3288BD", cex=0.5)
points(xCell_inmin_H23_1,mean_data_HCC44[,13], pch=16, col="grey", cex=0.2)
box(which="plot",lty="solid",lwd=3)
title( ylab="Cell Index [a.u.]", cex.lab=3.5,cex.main=6,line=8)
title( xlab="Time [h]", cex.lab=3.5,cex.main=6,line=7)
title(main="AZD1390 unradiated", cex.lab=3,cex.main=3,line=2)
axis(1,lwd=3,at = c(0,50,100),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 4,0))
axis(2,lwd=3,at = c(0,3,6),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 3,0),las=2)
minor.tick(nx =2/1,ny=2/1,tick.ratio =1,  y.args = list( lwd = 3,tcl = -1),  x.args = list( lwd = 3,tcl = -1)) 

plot(xCell_inmin_H23_1,mean_data_H23[,6], pch=16, cex=0.5,  ylim = c(-0.5,6),axes=FALSE, col="#FDAE61", cex.lab=4,cex.axis=3,mgp=c(7,1,0),ann=FALSE,cex.main=4)
points(xCell_inmin_H23_1,mean_data_H23[,8], pch=16, col="#66C2A5", cex=0.5)
points(xCell_inmin_H23_1,mean_data_H23[,13], pch=16, col="grey", cex=0.2)
points(xCell_inmin_H23_1,mean_data_HCC44[,6], pch=16, col="#D53E4F", cex=0.5)
points(xCell_inmin_H23_1,mean_data_HCC44[,8], pch=16, col="#3288BD", cex=0.5)
points(xCell_inmin_H23_1,mean_data_HCC44[,13], pch=16, col="grey50", cex=0.2)
box(which="plot",lty="solid",lwd=3)
axis(1,lwd=3,at = c(0,50,100),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 4,0))
axis(2,lwd=3,at = c(0,3,6),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 3,0),las=2)
title(main="AZD1390 2Gy", cex.lab=3,cex.main=3,line=2)
title( ylab="Cell Index [a.u.]", cex.lab=3.5,cex.main=6,line=8)
title( xlab="Time [h]", cex.lab=3.5,cex.main=6,line=7)
legend(c(-5,65),c(2.6,6.1),leg.txt,box.lwd=2, ncol=1,pch=16,cex=2.5,bty = "n", bg = NA,adj = c(0, 0.3), col=c("#66C2A5","#FDAE61","#3288BD","#D53E4F","grey"),y.intersp = 0.85,x.intersp = 0.75)
minor.tick(nx =2/1,ny=2/1,tick.ratio =1,  y.args = list( lwd = 3,tcl = -1),  x.args = list( lwd = 3,tcl = -1)) 
dev.off()
setwd(source_wd)