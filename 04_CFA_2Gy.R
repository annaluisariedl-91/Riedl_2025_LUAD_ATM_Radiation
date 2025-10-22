####Set directory#####
setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/CFA/Exp2_DoseResponse2Gy/raw")
dose <-c(0,0.5,1,1.5,2)
tdose<-sapply(dose,c)
#####loading raw data####
HCC44_1<- read.csv("CFA_1_2Gy_HCC44_rep1.csv",sep=",")
HCC44_2<- read.csv("CFA_1_2Gy_HCC44_rep2.csv",sep=",")
HCC44_4<- read.csv("CFA_1_2Gy_HCC44_rep3.csv",sep=",")
H23_1<- read.csv("CFA_1_2Gy_H23_rep1.csv",sep=",")
H23_2<- read.csv("CFA_1_2Gy_H23_rep2.csv",sep=",")
H23_3<- read.csv("CFA_1_2Gy_H23_rep3.csv",sep=",")
H23_4<- read.csv("CFA_1_2Gy_H23_rep4.csv",sep=",")
#####Caclulation HCC44 I####
HCC44_1 <- HCC44_1[2:31,1:3]
colnames(HCC44_1)<- c("Dose","N0","N")
HCC44_1["N0/N1"]<-as.numeric(HCC44_1[,3])/as.numeric(HCC44_1[,2])
HCC44_1["Plating effeciency"]<- (HCC44_1[,4]/mean(as.numeric(HCC44_1[1:6,4])))
cells_HCC44_1 <- data.frame("N0"=as.numeric(HCC44_1[c(1,7,13,19,25),2]))

HCC44_1_0Gy <- (as.numeric(HCC44_1[1:6,3:3]))
HCC44_1_05Gy <- (as.numeric(HCC44_1[7:12,3:3]))
HCC44_1_1Gy <- (as.numeric(HCC44_1[13:18,3:3]))
HCC44_1_1.5Gy <- (as.numeric(HCC44_1[19:24,3:3]))
HCC44_1_2Gy <- as.numeric(HCC44_1[25:30,3:3])
PE_HCC44_1<- mean(mean(HCC44_1_0Gy)/as.numeric(HCC44_1[1,2]))
mean_HCC44_1 <- data.frame("dose_Gy"=c(0,0.5,1,1.5,2),
                           "N0"=c(cells_HCC44_1),
                           "N_mean"= c(mean(HCC44_1_0Gy,na.rm=TRUE),mean(HCC44_1_05Gy,na.rm=TRUE),mean(HCC44_1_1Gy,na.rm=TRUE),mean(HCC44_1_1.5Gy,na.rm=TRUE),mean(HCC44_1_2Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC44_1_0Gy,na.rm=TRUE),sd(HCC44_1_05Gy,na.rm=TRUE),sd(HCC44_1_1Gy,na.rm=TRUE),sd(HCC44_1_1.5Gy,na.rm=TRUE),sd(HCC44_1_2Gy,na.rm=TRUE)))
mean_HCC44_1["SEM"] <- mean_HCC44_1[,4]/sqrt(6)
mean_HCC44_1["N/N0"] <- mean_HCC44_1[,3]/mean_HCC44_1[,2]
mean_HCC44_1["PE"] <- mean_HCC44_1[1,3]/mean_HCC44_1[1,2]
mean_HCC44_1["SF"] <- mean_HCC44_1["N_mean"]/(mean_HCC44_1["PE"]*mean_HCC44_1["N0"])

HCC44_1["SF"]<- as.numeric(unlist(HCC44_1["N"]))/(mean_HCC44_1[1,7]*as.numeric(unlist(HCC44_1["N0"])))
#####Caclulation HCC44 II####
HCC44_2 <- HCC44_2[2:31,1:3]
colnames(HCC44_2)<- c("Dose","N0","N")
HCC44_2["N0/N1"]<-as.numeric(HCC44_2[,3])/as.numeric(HCC44_2[,2])
HCC44_2["Plating effeciency"]<- (HCC44_2[,4]/mean(as.numeric(HCC44_2[1:6,4])))
cells_HCC44_2 <- data.frame("N0"=as.numeric(HCC44_2[c(1,7,13,19,25),2]))

HCC44_2_0Gy <- (as.numeric(HCC44_2[1:6,3:3]))
HCC44_2_05Gy <- (as.numeric(HCC44_2[7:12,3:3]))
HCC44_2_1Gy <- (as.numeric(HCC44_2[13:18,3:3]))
HCC44_2_1.5Gy <- (as.numeric(HCC44_2[19:24,3:3]))
HCC44_2_2Gy <- as.numeric(HCC44_2[25:30,3:3])
PE_HCC44_2<- mean(mean(HCC44_2_0Gy)/as.numeric(HCC44_2[1,2]))
mean_HCC44_2 <- data.frame("dose_Gy"=c(0,0.5,1,1.5,2),
                           "N0"=c(cells_HCC44_2),
                           "N_mean"= c(mean(HCC44_2_0Gy,na.rm=TRUE),mean(HCC44_2_05Gy,na.rm=TRUE),mean(HCC44_2_1Gy,na.rm=TRUE),mean(HCC44_2_1.5Gy,na.rm=TRUE),mean(HCC44_2_2Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC44_2_0Gy,na.rm=TRUE),sd(HCC44_2_05Gy,na.rm=TRUE),sd(HCC44_2_1Gy,na.rm=TRUE),sd(HCC44_2_1.5Gy,na.rm=TRUE),sd(HCC44_2_2Gy,na.rm=TRUE)))
mean_HCC44_2["SEM"] <- mean_HCC44_2[,4]/sqrt(6)
mean_HCC44_2["N/N0"] <- mean_HCC44_2[,3]/mean_HCC44_2[,2]
mean_HCC44_2["PE"] <- mean_HCC44_2[1,3]/mean_HCC44_2[1,2]
mean_HCC44_2["SF"] <- mean_HCC44_2["N_mean"]/(mean_HCC44_2["PE"]*mean_HCC44_2["N0"])

HCC44_2["SF"]<- as.numeric(unlist(HCC44_2["N"]))/(mean_HCC44_2[1,7]*as.numeric(unlist(HCC44_2["N0"])))
#####Caclulation HCC44 III####
HCC44_3 <- HCC44_3[2:31,1:3]
colnames(HCC44_3)<- c("Dose","N0","N")
HCC44_3["N0/N1"]<-as.numeric(HCC44_3[,3])/as.numeric(HCC44_3[,2])
HCC44_3["Plating effeciency"]<- (HCC44_3[,4]/mean(as.numeric(HCC44_3[1:6,4])))
cells_HCC44_3 <- data.frame("N0"=as.numeric(HCC44_3[c(1,7,13,19,25),2]))

HCC44_3_0Gy <- (as.numeric(HCC44_3[1:6,3:3]))
HCC44_3_05Gy <- (as.numeric(HCC44_3[7:12,3:3]))
HCC44_3_1Gy <- (as.numeric(HCC44_3[13:18,3:3]))
HCC44_3_1.5Gy <- (as.numeric(HCC44_3[19:24,3:3]))
HCC44_3_2Gy <- as.numeric(HCC44_3[25:30,3:3])
PE_HCC44_3<- mean(mean(HCC44_3_0Gy)/as.numeric(HCC44_3[1,2]))
mean_HCC44_3 <- data.frame("dose_Gy"=c(0,0.5,1,1.5,2),
                           "N0"=c(cells_HCC44_3),
                           "N_mean"= c(mean(HCC44_3_0Gy,na.rm=TRUE),mean(HCC44_3_05Gy,na.rm=TRUE),mean(HCC44_3_1Gy,na.rm=TRUE),mean(HCC44_3_1.5Gy,na.rm=TRUE),mean(HCC44_3_2Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC44_3_0Gy,na.rm=TRUE),sd(HCC44_3_05Gy,na.rm=TRUE),sd(HCC44_3_1Gy,na.rm=TRUE),sd(HCC44_3_1.5Gy,na.rm=TRUE),sd(HCC44_3_2Gy,na.rm=TRUE)))
mean_HCC44_3["SEM"] <- mean_HCC44_3[,4]/sqrt(6)
mean_HCC44_3["N/N0"] <- mean_HCC44_3[,3]/mean_HCC44_3[,2]
mean_HCC44_3["PE"] <- mean_HCC44_3[1,3]/mean_HCC44_3[1,2]
mean_HCC44_3["SF"] <- mean_HCC44_3["N_mean"]/(mean_HCC44_3["PE"]*mean_HCC44_3["N0"])

HCC44_3["SF"]<- as.numeric(unlist(HCC44_3["N"]))/(mean_HCC44_3[1,7]*as.numeric(unlist(HCC44_3["N0"])))
#####binding HCC44 results####
all_HCC44 <- cbind(mean_HCC44_1[,c(1:3,8)],mean_HCC44_2[,c(3,8)],mean_HCC44_3[,c(3,8)])
colnames(all_HCC44)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3")
all_HCC44["SF_mean"]<- rowMeans(all_HCC44[,c(4,6,8)],na.rm=TRUE)
all_HCC44["SF_sd"]<- apply(all_HCC44[,c(4,6,8)], 1, sd,na.rm=TRUE)
all_HCC44["SF_SEM"]<- all_HCC44["SF_sd"]/sqrt(3)

#####Caclulation NCI-H23 I####
H23_1 <- H23_1[2:15,1:3]
colnames(H23_1)<- c("Dose","N0","N")
H23_1["N0/N1"]<-as.numeric(H23_1[,3])/as.numeric(H23_1[,2])
H23_1["Plating effeciency"]<- (H23_1[,4]/mean(as.numeric(H23_1[1:6,4])))
cells_H23_1 <- data.frame("N0"=as.numeric(H23_1[c(1,4,7,10,13),2]))

H23_1_0Gy <- (as.numeric(H23_1[1:3,3:3]))
H23_1_05Gy <- (as.numeric(H23_1[4:6,3:3]))
H23_1_1Gy <- (as.numeric(H23_1[7:9,3:3]))
H23_1_1.5Gy <- (as.numeric(H23_1[10:12,3:3]))
H23_1_2Gy <- as.numeric(H23_1[13:15,3:3])
PE_H23_1<- mean(mean(H23_1_0Gy)/as.numeric(H23_1[1,2]))
mean_H23_1 <- data.frame("dose_Gy"=c(0,0.5,1,1.5,2),
                         "N0"=c(cells_H23_1),
                         "N_mean"= c(mean(H23_1_0Gy,na.rm=TRUE),mean(H23_1_05Gy,na.rm=TRUE),mean(H23_1_1Gy,na.rm=TRUE),mean(H23_1_1.5Gy,na.rm=TRUE),mean(H23_1_2Gy,na.rm=TRUE)),
                         "N_SD"= c(sd(H23_1_0Gy,na.rm=TRUE),sd(H23_1_05Gy,na.rm=TRUE),sd(H23_1_1Gy,na.rm=TRUE),sd(H23_1_1.5Gy,na.rm=TRUE),sd(H23_1_2Gy,na.rm=TRUE)))
mean_H23_1["SEM"] <- mean_H23_1[,4]/sqrt(3)
mean_H23_1["N/N0"] <- mean_H23_1[,3]/mean_H23_1[,2]
mean_H23_1["PE"] <- mean_H23_1[1,3]/mean_H23_1[1,2]
mean_H23_1["SF"] <- mean_H23_1["N_mean"]/(mean_H23_1["PE"]*mean_H23_1["N0"])

H23_1["SF"]<- as.numeric(unlist(H23_1["N"]))/(mean_H23_1[1,7]*as.numeric(unlist(H23_1["N0"])))
#####Caclulation NCI-H23 II####
H23_2 <- H23_2[2:31,1:3]
colnames(H23_2)<- c("Dose","N0","N")
H23_2["N0/N1"]<-as.numeric(H23_2[,3])/as.numeric(H23_2[,2])
H23_2["Plating effeciency"]<- (H23_2[,4]/mean(as.numeric(H23_2[1:6,4])))
cells_H23_2 <- data.frame("N0"=as.numeric(H23_2[c(1,7,10,17,25),2]))

H23_2_0Gy <- (as.numeric(H23_2[1:6,3:3]))
H23_2_05Gy <- (as.numeric(H23_2[7:12,3:3]))
H23_2_1Gy <- (as.numeric(H23_2[13:18,3:3]))
H23_2_1.5Gy <- (as.numeric(H23_2[19:24,3:3]))
H23_2_2Gy <- as.numeric(H23_2[25:30,3:3])
PE_H23_2<- mean(mean(H23_2_0Gy)/as.numeric(H23_2[1,2]))
mean_H23_2 <- data.frame("dose_Gy"=c(0,0.5,1,1.5,2),
                         "N0"=c(cells_H23_2),
                         "N_mean"= c(mean(H23_2_0Gy,na.rm=TRUE),mean(H23_2_05Gy,na.rm=TRUE),mean(H23_2_1Gy,na.rm=TRUE),mean(H23_2_1.5Gy,na.rm=TRUE),mean(H23_2_2Gy,na.rm=TRUE)),
                         "N_SD"= c(sd(H23_2_0Gy,na.rm=TRUE),sd(H23_2_05Gy,na.rm=TRUE),sd(H23_2_1Gy,na.rm=TRUE),sd(H23_2_1.5Gy,na.rm=TRUE),sd(H23_2_2Gy,na.rm=TRUE)))
mean_H23_2["SEM"] <- mean_H23_2[,4]/sqrt(6)
mean_H23_2["N/N0"] <- mean_H23_2[,3]/mean_H23_2[,2]
mean_H23_2["PE"] <- mean_H23_2[1,3]/mean_H23_2[1,2]
mean_H23_2["SF"] <- mean_H23_2["N_mean"]/(mean_H23_2["PE"]*mean_H23_2["N0"])

H23_2["SF"]<- as.numeric(unlist(H23_2["N"]))/(mean_H23_2[1,7]*as.numeric(unlist(H23_2["N0"])))
#####Caclulation NCI-H23 III####
H23_3 <- H23_3[2:31,1:3]
colnames(H23_3)<- c("Dose","N0","N")
H23_3["N0/N1"]<-as.numeric(H23_3[,3])/as.numeric(H23_3[,2])
H23_3["Plating effeciency"]<- (H23_3[,4]/mean(as.numeric(H23_3[1:6,4])))
cells_H23_3 <- data.frame("N0"=as.numeric(H23_3[c(1,7,10,17,25),2]))

H23_3_0Gy <- (as.numeric(H23_3[1:6,3:3]))
H23_3_05Gy <- (as.numeric(H23_3[7:12,3:3]))
H23_3_1Gy <- (as.numeric(H23_3[13:18,3:3]))
H23_3_1.5Gy <- (as.numeric(H23_3[19:24,3:3]))
H23_3_2Gy <- as.numeric(H23_3[25:30,3:3])
PE_H23_3<- mean(mean(H23_3_0Gy)/as.numeric(H23_3[1,2]))
mean_H23_3 <- data.frame("dose_Gy"=c(0,0.5,1,1.5,2),
                         "N0"=c(cells_H23_3),
                         "N_mean"= c(mean(H23_3_0Gy,na.rm=TRUE),mean(H23_3_05Gy,na.rm=TRUE),mean(H23_3_1Gy,na.rm=TRUE),mean(H23_3_1.5Gy,na.rm=TRUE),mean(H23_3_2Gy,na.rm=TRUE)),
                         "N_SD"= c(sd(H23_3_0Gy,na.rm=TRUE),sd(H23_3_05Gy,na.rm=TRUE),sd(H23_3_1Gy,na.rm=TRUE),sd(H23_3_1.5Gy,na.rm=TRUE),sd(H23_3_2Gy,na.rm=TRUE)))
mean_H23_3["SEM"] <- mean_H23_3[,4]/sqrt(6)
mean_H23_3["N/N0"] <- mean_H23_3[,3]/mean_H23_3[,2]
mean_H23_3["PE"] <- mean_H23_3[1,3]/mean_H23_3[1,2]
mean_H23_3["SF"] <- mean_H23_3["N_mean"]/(mean_H23_3["PE"]*mean_H23_3["N0"])

H23_3["SF"]<- as.numeric(unlist(H23_3["N"]))/(mean_H23_3[1,7]*as.numeric(unlist(H23_3["N0"])))
#####Caclulation NCI-H23 IV####
H23_4 <- H23_4[2:31,1:3]
colnames(H23_4)<- c("Dose","N0","N")
H23_4["N0/N1"]<-as.numeric(H23_4[,3])/as.numeric(H23_4[,2])
H23_4["Plating effeciency"]<- (H23_4[,4]/mean(as.numeric(H23_4[1:6,4])))
cells_H23_4 <- data.frame("N0"=as.numeric(H23_4[c(1,7,13,19,25),2]))

H23_4_0Gy <- (as.numeric(H23_4[1:6,3:3]))
H23_4_05Gy <- (as.numeric(H23_4[7:12,3:3]))
H23_4_1Gy <- (as.numeric(H23_4[13:18,3:3]))
H23_4_1.5Gy <- (as.numeric(H23_4[19:24,3:3]))
H23_4_2Gy <- as.numeric(H23_4[25:30,3:3])
PE_H23_4<- mean(mean(H23_4_0Gy)/as.numeric(H23_4[1,2]))
mean_H23_4 <- data.frame("dose_Gy"=c(0,0.5,1,1.5,2),
                         "N0"=c(cells_H23_4),
                         "N_mean"= c(mean(H23_4_0Gy,na.rm=TRUE),mean(H23_4_05Gy,na.rm=TRUE),mean(H23_4_1Gy,na.rm=TRUE),mean(H23_4_1.5Gy,na.rm=TRUE),mean(H23_4_2Gy,na.rm=TRUE)),
                         "N_SD"= c(sd(H23_4_0Gy,na.rm=TRUE),sd(H23_4_05Gy,na.rm=TRUE),sd(H23_4_1Gy,na.rm=TRUE),sd(H23_4_1.5Gy,na.rm=TRUE),sd(H23_4_2Gy,na.rm=TRUE)))
mean_H23_4["SEM"] <- mean_H23_4[,4]/sqrt(6)
mean_H23_4["N/N0"] <- mean_H23_4[,3]/mean_H23_4[,2]
mean_H23_4["PE"] <- mean_H23_4[1,3]/mean_H23_4[1,2]
mean_H23_4["SF"] <- mean_H23_4["N_mean"]/(mean_H23_4["PE"]*mean_H23_4["N0"])

H23_4["SF"]<- as.numeric(unlist(H23_4["N"]))/(mean_H23_4[1,7]*as.numeric(unlist(H23_4["N0"])))
#####binding H23 results####
all_H23 <- cbind(mean_H23_1[,c(1:3,8)],mean_H23_2[,c(3,8)],mean_H23_3[,c(3,8)],mean_H23_4[,c(3,8)])
colnames(all_H23)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_H23["SF_mean"]<- rowMeans(all_H23[,c(4,6,8,10)],na.rm=TRUE)
all_H23["SF_sd"]<- apply(all_H23[,c(4,6,8,10)], 1, sd,na.rm=TRUE)
all_H23["SF_SEM"]<- all_H23["SF_sd"]/sqrt(4)
#####Modelling HCC44#####
fit_HCC44 <- nlsLM(all_HCC44[,9] ~ exp(- (alpha * dose + beta * dose^2)), 
                                start = list(alpha = 0.3, beta = 0.03), 
                                lower = c(0,0),
                                data = data.frame(dose, all_HCC44))
alpha_HCC44 <- coef(fit_HCC44)["alpha"]
beta_HCC44<- coef(fit_HCC44)["beta"]
#####Modelling H23#####
fit_H23 <- nlsLM(all_H23[,11] ~ exp(- (alpha * dose + beta * dose^2)), 
                   start = list(alpha = 0.3, beta = 0.03), 
                   lower = c(0,0),
                   data = data.frame(dose, all_H23))
alpha_H23 <- coef(fit_H23)["alpha"]
beta_H23<- coef(fit_H23)["beta"]
####H23 and HCC44 - Plotting####
setwd(save_wd)
jpeg(filename = "CFA_2Gy_base.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1),ncol=3,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt<-c("HCC44","NCI-H23")
plot(x=dose, y=all_H23[,11],cex=5,ann=TRUE,axes=FALSE,lwd=2,cex.axis=3, mgp=c(7,1.5, 0), ylim=c(0.001,5),pch=21,log="y", col="black", bg="#66C2A5", lty=2,cex.main=4,xlab="",  ylab="",cex.lab=4)
points(x=dose, y=all_HCC44[,9],xaxt="n", lwd=2,ylim=c(0.01,50), cex=5, pch=22, col="black", bg="#3288BD", lty=2)
lines(dose,exp(-(dose^2*(alpha_H23))+dose*(beta_H23)),lwd=3,lty=3)
lines(dose,exp(-(dose^2*(alpha_HCC44))+dose*(beta_HCC44)),lwd=3,lty=3)
axis(1,lwd=3,at = c(0,0.5,1,1.5, 2),lwd.ticks = 4,cex=2, cex.axis=3,tcl = -1.5,mgp=c(3, 4,0))
axis(2,lwd=3,at = c(1, 0.1, 0.01,0.001),labels = c("1","","",parse(text="10^-3")),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,mgp=c(3, 2,0),las=2)
box(lty="solid",lwd=3)
arrows(x0=dose, y0=all_H23[,11]-all_H23[,13], x1=dose, y1=all_H23[,11]+all_H23[,13], code=3, angle=90, length=0.1,lwd=3)
arrows(x0=dose, y0=all_HCC44[,9]-all_HCC44[,11], x1=dose, y1=all_HCC44[,9]+all_HCC44[,11], code=3, angle=90, length=0.1,lwd=3)
title( ylab="Surviving Fraction", cex.lab=3,cex.main=6,line=9,adj=0.75)
title( xlab="Dose [Gy]", cex.lab=3,cex.main=6,line=7)
legend(c(0,1.6),c(0.004,0.009),leg.txt,box.col="white",bty = "n", bg = NA, ncol=1, pch=c(21,22),adj=c(0,0.5),pt.cex=3,pt.lwd=3, cex = 2.5,pt.bg=c("#3288BD","#66C2A5"),box.lwd=2,y.intersp = 1.1)
minor.tick(nx =2/1,ny=0,tick.ratio =1.5, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3)) 
dev.off()
setwd(source_wd)
