#####set directory#####
setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/CFA/Exp3_ATMinhibition/raw/HCC44")
dose_HCC44 <-c(0,1,2,4,8)
dose_H23<-c(0,0.5,1,2,4)
######loading raw data####
HCC44_DMSO_1<- read.csv("CFA_HCC44_DMSO_Rep1.csv",sep=",")
HCC44_DMSO_2<- read.csv("CFA_HCC44_DMSO_Rep2.csv",sep=",")
HCC44_1nM_AZD1390_1<- read.csv("CFA_HCC44_1nM_AZD1390_Rep1.csv",sep=",")
HCC44_1nM_AZD1390_2<- read.csv("CFA_HCC44_1nM_AZD1390_Rep2.csv",sep=",")
HCC44_10nM_AZD1390_1<- read.csv("CFA_HCC44_10nM_AZD1390_Rep1.csv",sep=",")
HCC44_10nM_AZD1390_2<- read.csv("CFA_HCC44_10nM_AZD1390_Rep2.csv",sep=",")

setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/CFA/Exp3_ATMinhibition/raw/NCI-H23")
H23_DMSO_1<- read.csv("CFA_H23_DMSO_Rep1.csv",sep=",")
H23_DMSO_2<- read.csv("CFA_H23_DMSO_Rep2.csv",sep=",")
H23_DMSO_3<- read.csv("CFA_H23_DMSO_Rep3.csv",sep=",")
H23_1nM_AZD1390_1<- read.csv("CFA_H23_1nM_AZD1390_Rep1.csv",sep=",")
H23_1nM_AZD1390_2<- read.csv("CFA_H23_1nM_AZD1390_Rep2.csv",sep=",")
H23_1nM_AZD1390_3<- read.csv("CFA_H23_1nM_AZD1390_Rep3.csv",sep=",")
H23_10nM_AZD1390_1<- read.csv("CFA_H23_10nM_AZD1390_Rep1.csv",sep=",")
H23_10nM_AZD1390_2<- read.csv("CFA_H23_10nM_AZD1390_Rep1.csv",sep=",")
H23_10nM_AZD1390_3<- read.csv("CFA_H23_10nM_AZD1390_Rep1.csv",sep=",")
#####Caclulation HCC44 DMSO I####
HCC44_DMSO_1 <- HCC44_DMSO_1[2:21,1:3]
colnames(HCC44_DMSO_1)<- c("Dose","N0","N")
HCC44_DMSO_1["N0/N1"]<-as.numeric(HCC44_DMSO_1[,3])/as.numeric(HCC44_DMSO_1[,2])
HCC44_DMSO_1["Plating effeciency"]<- (HCC44_DMSO_1[,4]/mean(as.numeric(HCC44_DMSO_1[1:4,4])))
cells_HCC44_DMSO_1 <- data.frame("N0"=as.numeric(HCC44_DMSO_1[c(1,5,9,13,17),2]))

HCC44_DMSO_1_0Gy <- (as.numeric(HCC44_DMSO_1[1:4,3:3]))
HCC44_DMSO_1_1Gy <- (as.numeric(HCC44_DMSO_1[5:8,3:3]))
HCC44_DMSO_1_2Gy <- (as.numeric(HCC44_DMSO_1[9:12,3:3]))
HCC44_DMSO_1_4Gy <- (as.numeric(HCC44_DMSO_1[13:16,3:3]))
HCC44_DMSO_1_8Gy <- as.numeric(HCC44_DMSO_1[17:20,3:3])
PE_HCC44_DMSO_1<- mean(mean(HCC44_DMSO_1_0Gy)/as.numeric(HCC44_DMSO_1[1,2]))
mean_HCC44_DMSO_1 <- data.frame("dose_Gy"=c(0,1,2,4,8),
                                "N0"=c(cells_HCC44_DMSO_1),
                                "N_mean"= c(mean(HCC44_DMSO_1_0Gy,na.rm=TRUE),mean(HCC44_DMSO_1_1Gy,na.rm=TRUE),mean(HCC44_DMSO_1_2Gy,na.rm=TRUE),mean(HCC44_DMSO_1_4Gy,na.rm=TRUE),mean(HCC44_DMSO_1_8Gy,na.rm=TRUE)),
                                "N_SD"= c(sd(HCC44_DMSO_1_0Gy,na.rm=TRUE),sd(HCC44_DMSO_1_1Gy,na.rm=TRUE),sd(HCC44_DMSO_1_2Gy,na.rm=TRUE),sd(HCC44_DMSO_1_4Gy,na.rm=TRUE),sd(HCC44_DMSO_1_8Gy,na.rm=TRUE)))
mean_HCC44_DMSO_1["SEM"] <- mean_HCC44_DMSO_1[,4]/sqrt(4)
mean_HCC44_DMSO_1["N/N0"] <- mean_HCC44_DMSO_1[,3]/mean_HCC44_DMSO_1[,2]
mean_HCC44_DMSO_1["PE"] <- mean_HCC44_DMSO_1[1,3]/mean_HCC44_DMSO_1[1,2]
mean_HCC44_DMSO_1["SF"] <- mean_HCC44_DMSO_1["N_mean"]/(mean_HCC44_DMSO_1["PE"]*mean_HCC44_DMSO_1["N0"])

HCC44_DMSO_1["SF"]<- as.numeric(unlist(HCC44_DMSO_1["N"]))/(mean_HCC44_DMSO_1[1,7]*as.numeric(unlist(HCC44_DMSO_1["N0"])))
#####Caclulation HCC44 1nM AZD1390 I####
HCC44_1nM_AZD1390_1 <- HCC44_1nM_AZD1390_1[2:21,1:3]
colnames(HCC44_1nM_AZD1390_1)<- c("Dose","N0","N")
HCC44_1nM_AZD1390_1["N0/N1"]<-as.numeric(HCC44_1nM_AZD1390_1[,3])/as.numeric(HCC44_1nM_AZD1390_1[,2])
HCC44_1nM_AZD1390_1["Plating effeciency"]<- (HCC44_1nM_AZD1390_1[,4]/mean(as.numeric(HCC44_1nM_AZD1390_1[1:4,4])))
cells_HCC44_1nM_AZD1390_1 <- data.frame("N0"=as.numeric(HCC44_1nM_AZD1390_1[c(1,5,9,13,17),2]))

HCC44_1nM_AZD1390_1_0Gy <- (as.numeric(HCC44_1nM_AZD1390_1[1:4,3:3]))
HCC44_1nM_AZD1390_1_1Gy <- (as.numeric(HCC44_1nM_AZD1390_1[5:8,3:3]))
HCC44_1nM_AZD1390_1_2Gy <- (as.numeric(HCC44_1nM_AZD1390_1[9:12,3:3]))
HCC44_1nM_AZD1390_1_4Gy <- (as.numeric(HCC44_1nM_AZD1390_1[13:16,3:3]))
HCC44_1nM_AZD1390_1_8Gy <- as.numeric(HCC44_1nM_AZD1390_1[17:20,3:3])
PE_HCC44_1nM_AZD1390_1<- mean(mean(HCC44_1nM_AZD1390_1_0Gy)/as.numeric(HCC44_1nM_AZD1390_1[1,2]))
mean_HCC44_1nM_AZD1390_1 <- data.frame("dose_Gy"=c(0,1,2,4,8),
                                       "N0"=c(cells_HCC44_1nM_AZD1390_1),
                                       "N_mean"= c(mean(HCC44_1nM_AZD1390_1_0Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_1_1Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_1_2Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_1_4Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_1_8Gy,na.rm=TRUE)),
                                       "N_SD"= c(sd(HCC44_1nM_AZD1390_1_0Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_1_1Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_1_2Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_1_4Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_1_8Gy,na.rm=TRUE)))
mean_HCC44_1nM_AZD1390_1["SEM"] <- mean_HCC44_1nM_AZD1390_1[,4]/sqrt(4)
mean_HCC44_1nM_AZD1390_1["N/N0"] <- mean_HCC44_1nM_AZD1390_1[,3]/mean_HCC44_1nM_AZD1390_1[,2]
mean_HCC44_1nM_AZD1390_1["PE"] <- mean_HCC44_1nM_AZD1390_1[1,3]/mean_HCC44_1nM_AZD1390_1[1,2]
mean_HCC44_1nM_AZD1390_1["SF"] <- mean_HCC44_1nM_AZD1390_1["N_mean"]/(mean_HCC44_1nM_AZD1390_1["PE"]*mean_HCC44_1nM_AZD1390_1["N0"])

HCC44_1nM_AZD1390_1["SF"]<- as.numeric(unlist(HCC44_1nM_AZD1390_1["N"]))/(mean_HCC44_1nM_AZD1390_1[1,7]*as.numeric(unlist(HCC44_1nM_AZD1390_1["N0"])))
#####Caclulation HCC44 10nM AZD1390 I####
HCC44_10nM_AZD1390_1 <- HCC44_10nM_AZD1390_1[2:21,1:3]
colnames(HCC44_10nM_AZD1390_1)<- c("Dose","N0","N")
HCC44_10nM_AZD1390_1["N0/N1"]<-as.numeric(HCC44_10nM_AZD1390_1[,3])/as.numeric(HCC44_10nM_AZD1390_1[,2])
HCC44_10nM_AZD1390_1["Plating effeciency"]<- (HCC44_10nM_AZD1390_1[,4]/mean(as.numeric(HCC44_10nM_AZD1390_1[1:4,4])))
cells_HCC44_10nM_AZD1390_1 <- data.frame("N0"=as.numeric(HCC44_10nM_AZD1390_1[c(1,5,9,13,17),2]))

HCC44_10nM_AZD1390_1_0Gy <- (as.numeric(HCC44_10nM_AZD1390_1[1:4,3:3]))
HCC44_10nM_AZD1390_1_1Gy <- (as.numeric(HCC44_10nM_AZD1390_1[5:8,3:3]))
HCC44_10nM_AZD1390_1_2Gy <- (as.numeric(HCC44_10nM_AZD1390_1[9:12,3:3]))
HCC44_10nM_AZD1390_1_4Gy <- (as.numeric(HCC44_10nM_AZD1390_1[13:16,3:3]))
HCC44_10nM_AZD1390_1_8Gy <- as.numeric(HCC44_10nM_AZD1390_1[17:20,3:3])
PE_HCC44_10nM_AZD1390_1<- mean(mean(HCC44_10nM_AZD1390_1_0Gy)/as.numeric(HCC44_10nM_AZD1390_1[1,2]))
mean_HCC44_10nM_AZD1390_1 <- data.frame("dose_Gy"=c(0,1,2,4,8),
                                        "N0"=c(cells_HCC44_10nM_AZD1390_1),
                                        "N_mean"= c(mean(HCC44_10nM_AZD1390_1_0Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_1_1Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_1_2Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_1_4Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_1_8Gy,na.rm=TRUE)),
                                        "N_SD"= c(sd(HCC44_10nM_AZD1390_1_0Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_1_1Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_1_2Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_1_4Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_1_8Gy,na.rm=TRUE)))
mean_HCC44_10nM_AZD1390_1["SEM"] <- mean_HCC44_10nM_AZD1390_1[,4]/sqrt(4)
mean_HCC44_10nM_AZD1390_1["N/N0"] <- mean_HCC44_10nM_AZD1390_1[,3]/mean_HCC44_10nM_AZD1390_1[,2]
mean_HCC44_10nM_AZD1390_1["PE"] <- mean_HCC44_10nM_AZD1390_1[1,3]/mean_HCC44_10nM_AZD1390_1[1,2]
mean_HCC44_10nM_AZD1390_1["SF"] <- mean_HCC44_10nM_AZD1390_1["N_mean"]/(mean_HCC44_10nM_AZD1390_1["PE"]*mean_HCC44_10nM_AZD1390_1["N0"])

HCC44_10nM_AZD1390_1["SF"]<- as.numeric(unlist(HCC44_10nM_AZD1390_1["N"]))/(mean_HCC44_10nM_AZD1390_1[1,7]*as.numeric(unlist(HCC44_10nM_AZD1390_1["N0"])))
#####Caclulation HCC44 DMSO II####
HCC44_DMSO_2 <- HCC44_DMSO_2[2:21,1:3]
colnames(HCC44_DMSO_2)<- c("Dose","N0","N")
HCC44_DMSO_2["N0/N1"]<-as.numeric(HCC44_DMSO_2[,3])/as.numeric(HCC44_DMSO_2[,2])
HCC44_DMSO_2["Plating effeciency"]<- (HCC44_DMSO_2[,4]/mean(as.numeric(HCC44_DMSO_2[1:4,4])))
cells_HCC44_DMSO_2 <- data.frame("N0"=as.numeric(HCC44_DMSO_2[c(1,5,9,13,17),2]))

HCC44_DMSO_2_0Gy <- (as.numeric(HCC44_DMSO_2[1:4,3:3]))
HCC44_DMSO_2_1Gy <- (as.numeric(HCC44_DMSO_2[5:8,3:3]))
HCC44_DMSO_2_2Gy <- (as.numeric(HCC44_DMSO_2[9:12,3:3]))
HCC44_DMSO_2_4Gy <- (as.numeric(HCC44_DMSO_2[13:16,3:3]))
HCC44_DMSO_2_8Gy <- as.numeric(HCC44_DMSO_2[17:20,3:3])
PE_HCC44_DMSO_2<- mean(mean(HCC44_DMSO_2_0Gy)/as.numeric(HCC44_DMSO_2[1,2]))
mean_HCC44_DMSO_2 <- data.frame("dose_Gy"=c(0,1,2,4,8),
                                "N0"=c(cells_HCC44_DMSO_2),
                                "N_mean"= c(mean(HCC44_DMSO_2_0Gy,na.rm=TRUE),mean(HCC44_DMSO_2_1Gy,na.rm=TRUE),mean(HCC44_DMSO_2_2Gy,na.rm=TRUE),mean(HCC44_DMSO_2_4Gy,na.rm=TRUE),mean(HCC44_DMSO_2_8Gy,na.rm=TRUE)),
                                "N_SD"= c(sd(HCC44_DMSO_2_0Gy,na.rm=TRUE),sd(HCC44_DMSO_2_1Gy,na.rm=TRUE),sd(HCC44_DMSO_2_2Gy,na.rm=TRUE),sd(HCC44_DMSO_2_4Gy,na.rm=TRUE),sd(HCC44_DMSO_2_8Gy,na.rm=TRUE)))
mean_HCC44_DMSO_2["SEM"] <- mean_HCC44_DMSO_2[,4]/sqrt(4)
mean_HCC44_DMSO_2["N/N0"] <- mean_HCC44_DMSO_2[,3]/mean_HCC44_DMSO_2[,2]
mean_HCC44_DMSO_2["PE"] <- mean_HCC44_DMSO_2[1,3]/mean_HCC44_DMSO_2[1,2]
mean_HCC44_DMSO_2["SF"] <- mean_HCC44_DMSO_2["N_mean"]/(mean_HCC44_DMSO_2["PE"]*mean_HCC44_DMSO_2["N0"])

HCC44_DMSO_2["SF"]<- as.numeric(unlist(HCC44_DMSO_2["N"]))/(mean_HCC44_DMSO_2[1,7]*as.numeric(unlist(HCC44_DMSO_2["N0"])))
#####Caclulation HCC44 1nM AZD1390 II####
HCC44_1nM_AZD1390_2 <- HCC44_1nM_AZD1390_2[2:21,1:3]
colnames(HCC44_1nM_AZD1390_2)<- c("Dose","N0","N")
HCC44_1nM_AZD1390_2["N0/N1"]<-as.numeric(HCC44_1nM_AZD1390_2[,3])/as.numeric(HCC44_1nM_AZD1390_2[,2])
HCC44_1nM_AZD1390_2["Plating effeciency"]<- (HCC44_1nM_AZD1390_2[,4]/mean(as.numeric(HCC44_1nM_AZD1390_2[1:4,4])))
cells_HCC44_1nM_AZD1390_2 <- data.frame("N0"=as.numeric(HCC44_1nM_AZD1390_2[c(1,5,9,13,17),2]))

HCC44_1nM_AZD1390_2_0Gy <- (as.numeric(HCC44_1nM_AZD1390_2[1:4,3:3]))
HCC44_1nM_AZD1390_2_1Gy <- (as.numeric(HCC44_1nM_AZD1390_2[5:8,3:3]))
HCC44_1nM_AZD1390_2_2Gy <- (as.numeric(HCC44_1nM_AZD1390_2[9:12,3:3]))
HCC44_1nM_AZD1390_2_4Gy <- (as.numeric(HCC44_1nM_AZD1390_2[13:16,3:3]))
HCC44_1nM_AZD1390_2_8Gy <- as.numeric(HCC44_1nM_AZD1390_2[17:20,3:3])
PE_HCC44_1nM_AZD1390_2<- mean(mean(HCC44_1nM_AZD1390_2_0Gy)/as.numeric(HCC44_1nM_AZD1390_2[1,2]))
mean_HCC44_1nM_AZD1390_2 <- data.frame("dose_Gy"=c(0,1,2,4,8),
                                       "N0"=c(cells_HCC44_1nM_AZD1390_2),
                                       "N_mean"= c(mean(HCC44_1nM_AZD1390_2_0Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_2_1Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_2_2Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_2_4Gy,na.rm=TRUE),mean(HCC44_1nM_AZD1390_2_8Gy,na.rm=TRUE)),
                                       "N_SD"= c(sd(HCC44_1nM_AZD1390_2_0Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_2_1Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_2_2Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_2_4Gy,na.rm=TRUE),sd(HCC44_1nM_AZD1390_2_8Gy,na.rm=TRUE)))
mean_HCC44_1nM_AZD1390_2["SEM"] <- mean_HCC44_1nM_AZD1390_2[,4]/sqrt(4)
mean_HCC44_1nM_AZD1390_2["N/N0"] <- mean_HCC44_1nM_AZD1390_2[,3]/mean_HCC44_1nM_AZD1390_2[,2]
mean_HCC44_1nM_AZD1390_2["PE"] <- mean_HCC44_1nM_AZD1390_2[1,3]/mean_HCC44_1nM_AZD1390_2[1,2]
mean_HCC44_1nM_AZD1390_2["SF"] <- mean_HCC44_1nM_AZD1390_2["N_mean"]/(mean_HCC44_1nM_AZD1390_2["PE"]*mean_HCC44_1nM_AZD1390_2["N0"])

HCC44_1nM_AZD1390_2["SF"]<- as.numeric(unlist(HCC44_1nM_AZD1390_2["N"]))/(mean_HCC44_1nM_AZD1390_2[1,7]*as.numeric(unlist(HCC44_1nM_AZD1390_2["N0"])))
#####Caclulation HCC44 10nM AZD1390 II####
HCC44_10nM_AZD1390_2 <- HCC44_10nM_AZD1390_2[2:21,1:3]
colnames(HCC44_10nM_AZD1390_2)<- c("Dose","N0","N")
HCC44_10nM_AZD1390_2["N0/N1"]<-as.numeric(HCC44_10nM_AZD1390_2[,3])/as.numeric(HCC44_10nM_AZD1390_2[,2])
HCC44_10nM_AZD1390_2["Plating effeciency"]<- (HCC44_10nM_AZD1390_2[,4]/mean(as.numeric(HCC44_10nM_AZD1390_2[1:4,4])))
cells_HCC44_10nM_AZD1390_2 <- data.frame("N0"=as.numeric(HCC44_10nM_AZD1390_2[c(1,5,9,13,17),2]))

HCC44_10nM_AZD1390_2_0Gy <- (as.numeric(HCC44_10nM_AZD1390_2[1:4,3:3]))
HCC44_10nM_AZD1390_2_1Gy <- (as.numeric(HCC44_10nM_AZD1390_2[5:8,3:3]))
HCC44_10nM_AZD1390_2_2Gy <- (as.numeric(HCC44_10nM_AZD1390_2[9:12,3:3]))
HCC44_10nM_AZD1390_2_4Gy <- (as.numeric(HCC44_10nM_AZD1390_2[13:16,3:3]))
HCC44_10nM_AZD1390_2_8Gy <- as.numeric(HCC44_10nM_AZD1390_2[17:20,3:3])
PE_HCC44_10nM_AZD1390_2<- mean(mean(HCC44_10nM_AZD1390_2_0Gy)/as.numeric(HCC44_10nM_AZD1390_2[1,2]))
mean_HCC44_10nM_AZD1390_2 <- data.frame("dose_Gy"=c(0,1,2,4,8),
                                        "N0"=c(cells_HCC44_10nM_AZD1390_2),
                                        "N_mean"= c(mean(HCC44_10nM_AZD1390_2_0Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_2_1Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_2_2Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_2_4Gy,na.rm=TRUE),mean(HCC44_10nM_AZD1390_2_8Gy,na.rm=TRUE)),
                                        "N_SD"= c(sd(HCC44_10nM_AZD1390_2_0Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_2_1Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_2_2Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_2_4Gy,na.rm=TRUE),sd(HCC44_10nM_AZD1390_1_8Gy,na.rm=TRUE)))
mean_HCC44_10nM_AZD1390_2["SEM"] <- mean_HCC44_10nM_AZD1390_2[,4]/sqrt(4)
mean_HCC44_10nM_AZD1390_2["N/N0"] <- mean_HCC44_10nM_AZD1390_2[,3]/mean_HCC44_10nM_AZD1390_2[,2]
mean_HCC44_10nM_AZD1390_2["PE"] <- mean_HCC44_10nM_AZD1390_2[1,3]/mean_HCC44_10nM_AZD1390_2[1,2]
mean_HCC44_10nM_AZD1390_2["SF"] <- mean_HCC44_10nM_AZD1390_2["N_mean"]/(mean_HCC44_10nM_AZD1390_2["PE"]*mean_HCC44_10nM_AZD1390_2["N0"])

HCC44_10nM_AZD1390_2["SF"]<- as.numeric(unlist(HCC44_10nM_AZD1390_2["N"]))/(mean_HCC44_10nM_AZD1390_2[1,7]*as.numeric(unlist(HCC44_10nM_AZD1390_2["N0"])))
#####binding HCC44 results####
all_HCC44_DMSO<- cbind(mean_HCC44_DMSO_1[,c(1:3,8)],mean_HCC44_DMSO_2[,c(3,8)])
colnames(all_HCC44_DMSO)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2")
all_HCC44_DMSO["SF_mean"]<- rowMeans(all_HCC44_DMSO[,c(4,6)])
all_HCC44_DMSO["SF_sd"]<- apply(all_HCC44_DMSO[,c(4,6)], 1, sd)
all_HCC44_DMSO["SF_SEM"]<- all_HCC44_DMSO["SF_sd"]/sqrt(2)

all_HCC44_1nM_AZD1390<- cbind(mean_HCC44_1nM_AZD1390_1[,c(1:3,8)],mean_HCC44_1nM_AZD1390_2[,c(3,8)])
colnames(all_HCC44_1nM_AZD1390)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2")
all_HCC44_1nM_AZD1390["SF_mean"]<- rowMeans(all_HCC44_1nM_AZD1390[,c(4,6)])
all_HCC44_1nM_AZD1390["SF_sd"]<- apply(all_HCC44_1nM_AZD1390[,c(4,6)], 1, sd)
all_HCC44_1nM_AZD1390["SF_SEM"]<- all_HCC44_1nM_AZD1390["SF_sd"]/sqrt(2)

all_HCC44_10nM_AZD1390<- cbind(mean_HCC44_10nM_AZD1390_1[,c(1:3,8)],mean_HCC44_10nM_AZD1390_2[,c(3,8)])
colnames(all_HCC44_10nM_AZD1390)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2")
all_HCC44_10nM_AZD1390["SF_mean"]<- rowMeans(all_HCC44_10nM_AZD1390[,c(4,6)])
all_HCC44_10nM_AZD1390["SF_sd"]<- apply(all_HCC44_10nM_AZD1390[,c(4,6)], 1, sd)
all_HCC44_10nM_AZD1390["SF_SEM"]<- all_HCC44_10nM_AZD1390["SF_sd"]/sqrt(2)
#####Caclulation H23 DMSO I####
H23_DMSO_1 <- H23_DMSO_1[2:21,1:3]
colnames(H23_DMSO_1)<- c("Dose","N0","N")
H23_DMSO_1["N0/N1"]<-as.numeric(H23_DMSO_1[,3])/as.numeric(H23_DMSO_1[,2])
H23_DMSO_1["Plating effeciency"]<- (H23_DMSO_1[,4]/mean(as.numeric(H23_DMSO_1[1:4,4])))
cells_H23_DMSO_1 <- data.frame("N0"=as.numeric(H23_DMSO_1[c(1,5,9,13,17),2]))

H23_DMSO_1_0Gy <- (as.numeric(H23_DMSO_1[1:4,3:3]))
H23_DMSO_1_0.5Gy <- (as.numeric(H23_DMSO_1[5:8,3:3]))
H23_DMSO_1_1Gy <- (as.numeric(H23_DMSO_1[9:12,3:3]))
H23_DMSO_1_2Gy <- (as.numeric(H23_DMSO_1[13:16,3:3]))
H23_DMSO_1_4Gy <- as.numeric(H23_DMSO_1[17:20,3:3])
PE_H23_DMSO_1<- mean(mean(H23_DMSO_1_0Gy)/as.numeric(H23_DMSO_1[1,2]))
mean_H23_DMSO_1 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                              "N0"=c(cells_H23_DMSO_1),
                              "N_mean"= c(mean(H23_DMSO_1_0Gy,na.rm=TRUE),mean(H23_DMSO_1_0.5Gy,na.rm=TRUE),mean(H23_DMSO_1_1Gy,na.rm=TRUE),mean(H23_DMSO_1_2Gy,na.rm=TRUE),mean(H23_DMSO_1_4Gy,na.rm=TRUE)),
                              "N_SD"= c(sd(H23_DMSO_1_0Gy,na.rm=TRUE),sd(H23_DMSO_1_0.5Gy,na.rm=TRUE),sd(H23_DMSO_1_1Gy,na.rm=TRUE),sd(H23_DMSO_1_2Gy,na.rm=TRUE),sd(H23_DMSO_1_4Gy,na.rm=TRUE)))
mean_H23_DMSO_1["SEM"] <- mean_H23_DMSO_1[,4]/sqrt(4)
mean_H23_DMSO_1["N/N0"] <- mean_H23_DMSO_1[,3]/mean_H23_DMSO_1[,2]
mean_H23_DMSO_1["PE"] <- mean_H23_DMSO_1[1,3]/mean_H23_DMSO_1[1,2]
mean_H23_DMSO_1["SF"] <- mean_H23_DMSO_1["N_mean"]/(mean_H23_DMSO_1["PE"]*mean_H23_DMSO_1["N0"])

H23_DMSO_1["SF"]<- as.numeric(unlist(H23_DMSO_1["N"]))/(mean_H23_DMSO_1[1,7]*as.numeric(unlist(H23_DMSO_1["N0"])))
#####Caclulation H23 1nM AZD1390 I####
H23_1nM_AZD1390_1 <- H23_1nM_AZD1390_1[2:21,1:3]
colnames(H23_1nM_AZD1390_1)<- c("Dose","N0","N")
H23_1nM_AZD1390_1["N0/N1"]<-as.numeric(H23_1nM_AZD1390_1[,3])/as.numeric(H23_1nM_AZD1390_1[,2])
H23_1nM_AZD1390_1["Plating effeciency"]<- (H23_1nM_AZD1390_1[,4]/mean(as.numeric(H23_1nM_AZD1390_1[1:4,4])))
cells_H23_1nM_AZD1390_1 <- data.frame("N0"=as.numeric(H23_1nM_AZD1390_1[c(1,5,9,13,17),2]))

H23_1nM_AZD1390_1_0Gy <- (as.numeric(H23_1nM_AZD1390_1[1:4,3:3]))
H23_1nM_AZD1390_1_0.5Gy <- (as.numeric(H23_1nM_AZD1390_1[5:8,3:3]))
H23_1nM_AZD1390_1_1Gy <- (as.numeric(H23_1nM_AZD1390_1[9:12,3:3]))
H23_1nM_AZD1390_1_2Gy <- (as.numeric(H23_1nM_AZD1390_1[13:16,3:3]))
H23_1nM_AZD1390_1_4Gy <- as.numeric(H23_1nM_AZD1390_1[17:20,3:3])
PE_H23_1nM_AZD1390_1<- mean(mean(H23_1nM_AZD1390_1_0Gy)/as.numeric(H23_1nM_AZD1390_1[1,2]))
mean_H23_1nM_AZD1390_1 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                                     "N0"=c(cells_H23_1nM_AZD1390_1),
                                     "N_mean"= c(mean(H23_1nM_AZD1390_1_0Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_1_0.5Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_1_1Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_1_2Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_1_4Gy,na.rm=TRUE)),
                                     "N_SD"= c(sd(H23_1nM_AZD1390_1_0Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_1_0.5Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_1_1Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_1_2Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_1_4Gy,na.rm=TRUE)))
mean_H23_1nM_AZD1390_1["SEM"] <- mean_H23_1nM_AZD1390_1[,4]/sqrt(4)
mean_H23_1nM_AZD1390_1["N/N0"] <- mean_H23_1nM_AZD1390_1[,3]/mean_H23_1nM_AZD1390_1[,2]
mean_H23_1nM_AZD1390_1["PE"] <- mean_H23_1nM_AZD1390_1[1,3]/mean_H23_1nM_AZD1390_1[1,2]
mean_H23_1nM_AZD1390_1["SF"] <- mean_H23_1nM_AZD1390_1["N_mean"]/(mean_H23_1nM_AZD1390_1["PE"]*mean_H23_1nM_AZD1390_1["N0"])

H23_1nM_AZD1390_1["SF"]<- as.numeric(unlist(H23_1nM_AZD1390_1["N"]))/(mean_H23_1nM_AZD1390_1[1,7]*as.numeric(unlist(H23_1nM_AZD1390_1["N0"])))
#####Caclulation H23 10nM AZD1390 I####
H23_10nM_AZD1390_1 <- H23_10nM_AZD1390_1[2:21,1:3]
colnames(H23_10nM_AZD1390_1)<- c("Dose","N0","N")
H23_10nM_AZD1390_1["N0/N1"]<-as.numeric(H23_10nM_AZD1390_1[,3])/as.numeric(H23_10nM_AZD1390_1[,2])
H23_10nM_AZD1390_1["Plating effeciency"]<- (H23_10nM_AZD1390_1[,4]/mean(as.numeric(H23_10nM_AZD1390_1[1:4,4])))
cells_H23_10nM_AZD1390_1 <- data.frame("N0"=as.numeric(H23_10nM_AZD1390_1[c(1,5,9,13,17),2]))

H23_10nM_AZD1390_1_0Gy <- (as.numeric(H23_10nM_AZD1390_1[1:4,3:3]))
H23_10nM_AZD1390_1_0.5Gy <- (as.numeric(H23_10nM_AZD1390_1[5:8,3:3]))
H23_10nM_AZD1390_1_1Gy <- (as.numeric(H23_10nM_AZD1390_1[9:12,3:3]))
H23_10nM_AZD1390_1_2Gy <- (as.numeric(H23_10nM_AZD1390_1[13:16,3:3]))
H23_10nM_AZD1390_1_4Gy <- as.numeric(H23_10nM_AZD1390_1[17:20,3:3])
PE_H23_10nM_AZD1390_1<- mean(mean(H23_10nM_AZD1390_1_0Gy)/as.numeric(H23_10nM_AZD1390_1[1,2]))
mean_H23_10nM_AZD1390_1 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                                      "N0"=c(cells_H23_10nM_AZD1390_1),
                                      "N_mean"= c(mean(H23_10nM_AZD1390_1_0Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_1_0.5Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_1_1Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_1_2Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_1_4Gy,na.rm=TRUE)),
                                      "N_SD"= c(sd(H23_10nM_AZD1390_1_0Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_1_0.5Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_1_1Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_1_2Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_1_4Gy,na.rm=TRUE)))
mean_H23_10nM_AZD1390_1["SEM"] <- mean_H23_10nM_AZD1390_1[,4]/sqrt(4)
mean_H23_10nM_AZD1390_1["N/N0"] <- mean_H23_10nM_AZD1390_1[,3]/mean_H23_10nM_AZD1390_1[,2]
mean_H23_10nM_AZD1390_1["PE"] <- mean_H23_10nM_AZD1390_1[1,3]/mean_H23_10nM_AZD1390_1[1,2]
mean_H23_10nM_AZD1390_1["SF"] <- mean_H23_10nM_AZD1390_1["N_mean"]/(mean_H23_10nM_AZD1390_1["PE"]*mean_H23_10nM_AZD1390_1["N0"])

H23_10nM_AZD1390_1["SF"]<- as.numeric(unlist(H23_10nM_AZD1390_1["N"]))/(mean_H23_10nM_AZD1390_1[1,7]*as.numeric(unlist(H23_10nM_AZD1390_1["N0"])))
#####Caclulation H23 DMSO II####
H23_DMSO_2 <- H23_DMSO_2[2:21,1:3]
colnames(H23_DMSO_2)<- c("Dose","N0","N")
H23_DMSO_2["N0/N1"]<-as.numeric(H23_DMSO_2[,3])/as.numeric(H23_DMSO_2[,2])
H23_DMSO_2["Plating effeciency"]<- (H23_DMSO_2[,4]/mean(as.numeric(H23_DMSO_2[1:4,4])))
cells_H23_DMSO_2 <- data.frame("N0"=as.numeric(H23_DMSO_2[c(1,5,9,13,17),2]))

H23_DMSO_2_0Gy <- (as.numeric(H23_DMSO_2[1:4,3:3]))
H23_DMSO_2_0.5Gy <- (as.numeric(H23_DMSO_2[5:8,3:3]))
H23_DMSO_2_1Gy <- (as.numeric(H23_DMSO_2[9:12,3:3]))
H23_DMSO_2_2Gy <- (as.numeric(H23_DMSO_2[13:16,3:3]))
H23_DMSO_2_4Gy <- as.numeric(H23_DMSO_2[17:20,3:3])
PE_H23_DMSO_2<- mean(mean(H23_DMSO_2_0Gy)/as.numeric(H23_DMSO_2[1,2]))
mean_H23_DMSO_2 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                              "N0"=c(cells_H23_DMSO_2),
                              "N_mean"= c(mean(H23_DMSO_2_0Gy,na.rm=TRUE),mean(H23_DMSO_2_0.5Gy,na.rm=TRUE),mean(H23_DMSO_2_1Gy,na.rm=TRUE),mean(H23_DMSO_2_2Gy,na.rm=TRUE),mean(H23_DMSO_2_4Gy,na.rm=TRUE)),
                              "N_SD"= c(sd(H23_DMSO_2_0Gy,na.rm=TRUE),sd(H23_DMSO_2_0.5Gy,na.rm=TRUE),sd(H23_DMSO_2_1Gy,na.rm=TRUE),sd(H23_DMSO_2_2Gy,na.rm=TRUE),sd(H23_DMSO_2_4Gy,na.rm=TRUE)))
mean_H23_DMSO_2["SEM"] <- mean_H23_DMSO_2[,4]/sqrt(4)
mean_H23_DMSO_2["N/N0"] <- mean_H23_DMSO_2[,3]/mean_H23_DMSO_2[,2]
mean_H23_DMSO_2["PE"] <- mean_H23_DMSO_2[1,3]/mean_H23_DMSO_2[1,2]
mean_H23_DMSO_2["SF"] <- mean_H23_DMSO_2["N_mean"]/(mean_H23_DMSO_2["PE"]*mean_H23_DMSO_2["N0"])

H23_DMSO_2["SF"]<- as.numeric(unlist(H23_DMSO_2["N"]))/(mean_H23_DMSO_2[1,7]*as.numeric(unlist(H23_DMSO_2["N0"])))
#####Caclulation H23 1nM AZD1390 II####
H23_1nM_AZD1390_2 <- H23_1nM_AZD1390_2[2:21,1:3]
colnames(H23_1nM_AZD1390_2)<- c("Dose","N0","N")
H23_1nM_AZD1390_2["N0/N1"]<-as.numeric(H23_1nM_AZD1390_2[,3])/as.numeric(H23_1nM_AZD1390_2[,2])
H23_1nM_AZD1390_2["Plating effeciency"]<- (H23_1nM_AZD1390_2[,4]/mean(as.numeric(H23_1nM_AZD1390_2[1:4,4])))
cells_H23_1nM_AZD1390_2 <- data.frame("N0"=as.numeric(H23_1nM_AZD1390_2[c(1,5,9,13,17),2]))

H23_1nM_AZD1390_2_0Gy <- (as.numeric(H23_1nM_AZD1390_2[1:4,3:3]))
H23_1nM_AZD1390_2_0.5Gy <- (as.numeric(H23_1nM_AZD1390_2[5:8,3:3]))
H23_1nM_AZD1390_2_1Gy <- (as.numeric(H23_1nM_AZD1390_2[9:12,3:3]))
H23_1nM_AZD1390_2_2Gy <- (as.numeric(H23_1nM_AZD1390_2[13:16,3:3]))
H23_1nM_AZD1390_2_4Gy <- as.numeric(H23_1nM_AZD1390_2[17:20,3:3])
PE_H23_1nM_AZD1390_2<- mean(mean(H23_1nM_AZD1390_2_0Gy)/as.numeric(H23_1nM_AZD1390_2[1,2]))
mean_H23_1nM_AZD1390_2 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                                     "N0"=c(cells_H23_1nM_AZD1390_2),
                                     "N_mean"= c(mean(H23_1nM_AZD1390_2_0Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_2_0.5Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_2_1Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_2_2Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_2_4Gy,na.rm=TRUE)),
                                     "N_SD"= c(sd(H23_1nM_AZD1390_2_0Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_2_0.5Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_2_1Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_2_2Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_2_4Gy,na.rm=TRUE)))
mean_H23_1nM_AZD1390_2["SEM"] <- mean_H23_1nM_AZD1390_2[,4]/sqrt(4)
mean_H23_1nM_AZD1390_2["N/N0"] <- mean_H23_1nM_AZD1390_2[,3]/mean_H23_1nM_AZD1390_2[,2]
mean_H23_1nM_AZD1390_2["PE"] <- mean_H23_1nM_AZD1390_2[1,3]/mean_H23_1nM_AZD1390_2[1,2]
mean_H23_1nM_AZD1390_2["SF"] <- mean_H23_1nM_AZD1390_2["N_mean"]/(mean_H23_1nM_AZD1390_2["PE"]*mean_H23_1nM_AZD1390_2["N0"])

H23_1nM_AZD1390_2["SF"]<- as.numeric(unlist(H23_1nM_AZD1390_2["N"]))/(mean_H23_1nM_AZD1390_2[1,7]*as.numeric(unlist(H23_1nM_AZD1390_2["N0"])))
#####Caclulation H23 10nM AZD1390 II####
H23_10nM_AZD1390_2 <- H23_10nM_AZD1390_2[2:21,1:3]
colnames(H23_10nM_AZD1390_2)<- c("Dose","N0","N")
H23_10nM_AZD1390_2["N0/N1"]<-as.numeric(H23_10nM_AZD1390_2[,3])/as.numeric(H23_10nM_AZD1390_2[,2])
H23_10nM_AZD1390_2["Plating effeciency"]<- (H23_10nM_AZD1390_2[,4]/mean(as.numeric(H23_10nM_AZD1390_2[1:4,4])))
cells_H23_10nM_AZD1390_2 <- data.frame("N0"=as.numeric(H23_10nM_AZD1390_2[c(1,5,9,13,17),2]))

H23_10nM_AZD1390_2_0Gy <- (as.numeric(H23_10nM_AZD1390_2[1:4,3:3]))
H23_10nM_AZD1390_2_0.5Gy <- (as.numeric(H23_10nM_AZD1390_2[5:8,3:3]))
H23_10nM_AZD1390_2_1Gy <- (as.numeric(H23_10nM_AZD1390_2[9:12,3:3]))
H23_10nM_AZD1390_2_2Gy <- (as.numeric(H23_10nM_AZD1390_2[13:16,3:3]))
H23_10nM_AZD1390_2_4Gy <- as.numeric(H23_10nM_AZD1390_2[17:20,3:3])
PE_H23_10nM_AZD1390_2<- mean(mean(H23_10nM_AZD1390_2_0Gy)/as.numeric(H23_10nM_AZD1390_2[1,2]))
mean_H23_10nM_AZD1390_2 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                                      "N0"=c(cells_H23_10nM_AZD1390_2),
                                      "N_mean"= c(mean(H23_10nM_AZD1390_2_0Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_2_0.5Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_2_1Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_2_2Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_2_4Gy,na.rm=TRUE)),
                                      "N_SD"= c(sd(H23_10nM_AZD1390_2_0Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_2_0.5Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_2_1Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_2_2Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_1_4Gy,na.rm=TRUE)))
mean_H23_10nM_AZD1390_2["SEM"] <- mean_H23_10nM_AZD1390_2[,4]/sqrt(4)
mean_H23_10nM_AZD1390_2["N/N0"] <- mean_H23_10nM_AZD1390_2[,3]/mean_H23_10nM_AZD1390_2[,2]
mean_H23_10nM_AZD1390_2["PE"] <- mean_H23_10nM_AZD1390_2[1,3]/mean_H23_10nM_AZD1390_2[1,2]
mean_H23_10nM_AZD1390_2["SF"] <- mean_H23_10nM_AZD1390_2["N_mean"]/(mean_H23_10nM_AZD1390_2["PE"]*mean_H23_10nM_AZD1390_2["N0"])

H23_10nM_AZD1390_2["SF"]<- as.numeric(unlist(H23_10nM_AZD1390_2["N"]))/(mean_H23_10nM_AZD1390_2[1,7]*as.numeric(unlist(H23_10nM_AZD1390_2["N0"])))
#####Caclulation H23 DMSO III####
H23_DMSO_3 <- H23_DMSO_3[2:21,1:3]
colnames(H23_DMSO_3)<- c("Dose","N0","N")
H23_DMSO_3["N0/N1"]<-as.numeric(H23_DMSO_3[,3])/as.numeric(H23_DMSO_3[,2])
H23_DMSO_3["Plating effeciency"]<- (H23_DMSO_3[,4]/mean(as.numeric(H23_DMSO_3[1:4,4])))
cells_H23_DMSO_3 <- data.frame("N0"=as.numeric(H23_DMSO_3[c(1,5,9,13,17),2]))

H23_DMSO_3_0Gy <- (as.numeric(H23_DMSO_3[1:4,3:3]))
H23_DMSO_3_0.5Gy <- (as.numeric(H23_DMSO_3[5:8,3:3]))
H23_DMSO_3_1Gy <- (as.numeric(H23_DMSO_3[9:12,3:3]))
H23_DMSO_3_2Gy <- (as.numeric(H23_DMSO_3[13:16,3:3]))
H23_DMSO_3_4Gy <- as.numeric(H23_DMSO_3[17:20,3:3])
PE_H23_DMSO_3<- mean(mean(H23_DMSO_3_0Gy)/as.numeric(H23_DMSO_3[1,2]))
mean_H23_DMSO_3 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                              "N0"=c(cells_H23_DMSO_3),
                              "N_mean"= c(mean(H23_DMSO_3_0Gy,na.rm=TRUE),mean(H23_DMSO_3_0.5Gy,na.rm=TRUE),mean(H23_DMSO_3_1Gy,na.rm=TRUE),mean(H23_DMSO_3_2Gy,na.rm=TRUE),mean(H23_DMSO_3_4Gy,na.rm=TRUE)),
                              "N_SD"= c(sd(H23_DMSO_3_0Gy,na.rm=TRUE),sd(H23_DMSO_3_0.5Gy,na.rm=TRUE),sd(H23_DMSO_3_1Gy,na.rm=TRUE),sd(H23_DMSO_3_2Gy,na.rm=TRUE),sd(H23_DMSO_3_4Gy,na.rm=TRUE)))
mean_H23_DMSO_3["SEM"] <- mean_H23_DMSO_3[,4]/sqrt(4)
mean_H23_DMSO_3["N/N0"] <- mean_H23_DMSO_3[,3]/mean_H23_DMSO_3[,2]
mean_H23_DMSO_3["PE"] <- mean_H23_DMSO_3[1,3]/mean_H23_DMSO_3[1,2]
mean_H23_DMSO_3["SF"] <- mean_H23_DMSO_3["N_mean"]/(mean_H23_DMSO_3["PE"]*mean_H23_DMSO_3["N0"])

H23_DMSO_3["SF"]<- as.numeric(unlist(H23_DMSO_3["N"]))/(mean_H23_DMSO_3[1,7]*as.numeric(unlist(H23_DMSO_3["N0"])))
#####Caclulation H23 1nM AZD1390 III####
H23_1nM_AZD1390_3 <- H23_1nM_AZD1390_3[2:21,1:3]
colnames(H23_1nM_AZD1390_3)<- c("Dose","N0","N")
H23_1nM_AZD1390_3["N0/N1"]<-as.numeric(H23_1nM_AZD1390_3[,3])/as.numeric(H23_1nM_AZD1390_3[,2])
H23_1nM_AZD1390_3["Plating effeciency"]<- (H23_1nM_AZD1390_3[,4]/mean(as.numeric(H23_1nM_AZD1390_3[1:4,4])))
cells_H23_1nM_AZD1390_3 <- data.frame("N0"=as.numeric(H23_1nM_AZD1390_3[c(1,5,9,13,17),2]))

H23_1nM_AZD1390_3_0Gy <- (as.numeric(H23_1nM_AZD1390_3[1:4,3:3]))
H23_1nM_AZD1390_3_0.5Gy <- (as.numeric(H23_1nM_AZD1390_3[5:8,3:3]))
H23_1nM_AZD1390_3_1Gy <- (as.numeric(H23_1nM_AZD1390_3[9:12,3:3]))
H23_1nM_AZD1390_3_2Gy <- (as.numeric(H23_1nM_AZD1390_3[13:16,3:3]))
H23_1nM_AZD1390_3_4Gy <- as.numeric(H23_1nM_AZD1390_3[17:20,3:3])
PE_H23_1nM_AZD1390_3<- mean(mean(H23_1nM_AZD1390_3_0Gy)/as.numeric(H23_1nM_AZD1390_3[1,2]))
mean_H23_1nM_AZD1390_3 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                                     "N0"=c(cells_H23_1nM_AZD1390_3),
                                     "N_mean"= c(mean(H23_1nM_AZD1390_3_0Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_3_0.5Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_3_1Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_3_2Gy,na.rm=TRUE),mean(H23_1nM_AZD1390_3_4Gy,na.rm=TRUE)),
                                     "N_SD"= c(sd(H23_1nM_AZD1390_3_0Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_3_0.5Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_3_1Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_3_2Gy,na.rm=TRUE),sd(H23_1nM_AZD1390_3_4Gy,na.rm=TRUE)))
mean_H23_1nM_AZD1390_3["SEM"] <- mean_H23_1nM_AZD1390_3[,4]/sqrt(4)
mean_H23_1nM_AZD1390_3["N/N0"] <- mean_H23_1nM_AZD1390_3[,3]/mean_H23_1nM_AZD1390_3[,2]
mean_H23_1nM_AZD1390_3["PE"] <- mean_H23_1nM_AZD1390_3[1,3]/mean_H23_1nM_AZD1390_3[1,2]
mean_H23_1nM_AZD1390_3["SF"] <- mean_H23_1nM_AZD1390_3["N_mean"]/(mean_H23_1nM_AZD1390_3["PE"]*mean_H23_1nM_AZD1390_3["N0"])

H23_1nM_AZD1390_3["SF"]<- as.numeric(unlist(H23_1nM_AZD1390_3["N"]))/(mean_H23_1nM_AZD1390_3[1,7]*as.numeric(unlist(H23_1nM_AZD1390_3["N0"])))
#####Caclulation H23 10nM AZD1390 III####
H23_10nM_AZD1390_3 <- H23_10nM_AZD1390_3[2:21,1:3]
colnames(H23_10nM_AZD1390_3)<- c("Dose","N0","N")
H23_10nM_AZD1390_3["N0/N1"]<-as.numeric(H23_10nM_AZD1390_3[,3])/as.numeric(H23_10nM_AZD1390_3[,2])
H23_10nM_AZD1390_3["Plating effeciency"]<- (H23_10nM_AZD1390_3[,4]/mean(as.numeric(H23_10nM_AZD1390_3[1:4,4])))
cells_H23_10nM_AZD1390_3 <- data.frame("N0"=as.numeric(H23_10nM_AZD1390_3[c(1,5,9,13,17),2]))

H23_10nM_AZD1390_3_0Gy <- (as.numeric(H23_10nM_AZD1390_3[1:4,3:3]))
H23_10nM_AZD1390_3_0.5Gy <- (as.numeric(H23_10nM_AZD1390_3[5:8,3:3]))
H23_10nM_AZD1390_3_1Gy <- (as.numeric(H23_10nM_AZD1390_3[9:12,3:3]))
H23_10nM_AZD1390_3_2Gy <- (as.numeric(H23_10nM_AZD1390_3[13:16,3:3]))
H23_10nM_AZD1390_3_4Gy <- as.numeric(H23_10nM_AZD1390_3[17:20,3:3])
PE_H23_10nM_AZD1390_3<- mean(mean(H23_10nM_AZD1390_3_0Gy)/as.numeric(H23_10nM_AZD1390_3[1,2]))
mean_H23_10nM_AZD1390_3 <- data.frame("dose_Gy"=c(0,0.5,1,2,4),
                                      "N0"=c(cells_H23_10nM_AZD1390_3),
                                      "N_mean"= c(mean(H23_10nM_AZD1390_3_0Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_3_0.5Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_3_1Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_3_2Gy,na.rm=TRUE),mean(H23_10nM_AZD1390_3_4Gy,na.rm=TRUE)),
                                      "N_SD"= c(sd(H23_10nM_AZD1390_3_0Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_3_0.5Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_3_1Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_3_2Gy,na.rm=TRUE),sd(H23_10nM_AZD1390_1_4Gy,na.rm=TRUE)))
mean_H23_10nM_AZD1390_3["SEM"] <- mean_H23_10nM_AZD1390_3[,4]/sqrt(4)
mean_H23_10nM_AZD1390_3["N/N0"] <- mean_H23_10nM_AZD1390_3[,3]/mean_H23_10nM_AZD1390_3[,2]
mean_H23_10nM_AZD1390_3["PE"] <- mean_H23_10nM_AZD1390_3[1,3]/mean_H23_10nM_AZD1390_3[1,2]
mean_H23_10nM_AZD1390_3["SF"] <- mean_H23_10nM_AZD1390_3["N_mean"]/(mean_H23_10nM_AZD1390_3["PE"]*mean_H23_10nM_AZD1390_3["N0"])

H23_10nM_AZD1390_3["SF"]<- as.numeric(unlist(H23_10nM_AZD1390_3["N"]))/(mean_H23_10nM_AZD1390_3[1,7]*as.numeric(unlist(H23_10nM_AZD1390_3["N0"])))
#####binding H23 results####
all_H23_DMSO<- cbind(mean_H23_DMSO_1[,c(1:3,8)],mean_H23_DMSO_2[,c(3,8)],mean_H23_DMSO_3[,c(3,8)])
colnames(all_H23_DMSO)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3")
all_H23_DMSO["SF_mean"]<- rowMeans(all_H23_DMSO[,c(4,6,8)])
all_H23_DMSO["SF_sd"]<- apply(all_H23_DMSO[,c(4,6,8)], 1, sd)
all_H23_DMSO["SF_SEM"]<- all_H23_DMSO["SF_sd"]/sqrt(2)

all_H23_1nM_AZD1390<- cbind(mean_H23_1nM_AZD1390_1[,c(1:3,8)],mean_H23_1nM_AZD1390_2[,c(3,8)],mean_H23_1nM_AZD1390_3[,c(3,8)])
colnames(all_H23_1nM_AZD1390)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3")
all_H23_1nM_AZD1390["SF_mean"]<- rowMeans(all_H23_1nM_AZD1390[,c(4,6,8)])
all_H23_1nM_AZD1390["SF_sd"]<- apply(all_H23_1nM_AZD1390[,c(4,6,8)], 1, sd)
all_H23_1nM_AZD1390["SF_SEM"]<- all_H23_1nM_AZD1390["SF_sd"]/sqrt(2)

all_H23_10nM_AZD1390<- cbind(mean_H23_10nM_AZD1390_1[,c(1:3,8)],mean_H23_10nM_AZD1390_2[,c(3,8)],mean_H23_10nM_AZD1390_3[,c(3,8)])
colnames(all_H23_10nM_AZD1390)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3")
all_H23_10nM_AZD1390["SF_mean"]<- rowMeans(all_H23_10nM_AZD1390[,c(4,6,8)])
all_H23_10nM_AZD1390["SF_sd"]<- apply(all_H23_10nM_AZD1390[,c(4,6,8)], 1, sd)
all_H23_10nM_AZD1390["SF_SEM"]<- all_H23_10nM_AZD1390["SF_sd"]/sqrt(2)
#####Graphs####
setwd(save_wd)
jpeg(filename = "CFA_AZD1390_1_10_nM.jpeg", width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1,2,2,2, 1,1,1,2,2,2, 1,1,1,2,2,2, 1,1,1,2,2,2),ncol=6,nrow = 4,byrow = TRUE))
par(mar = c(10, 14, 6, 10),xpd=TRUE)
leg.txt <- c("DMSO","1nM AZD1390","10nM AZD1390")
plot(x=dose_HCC44, y=all_HCC44_1nM_AZD1390[,7],xlim=c(0,8),mgp=c(6, 0,1), type="l",cex.lab=4,lwd=3,ylim=c(0.0001,1000),axes=FALSE, cex=2, 
     pch=22,log="y", col="black", bg="#F46D43", lty=3, cex.main=4, ann=FALSE)
lines(x=dose_HCC44, y=all_HCC44_DMSO[,7], lwd=3,ylim=c(0.0001,100), cex=3, pch=22, col="black", bg="#ABDDA4", lty=3)
lines(x=dose_HCC44[1:3], y=all_HCC44_10nM_AZD1390[1:3,7], lwd=3,ylim=c(0.0001,10), cex=3, pch=21, col="black", lty=3)
lines(x=c(2,4), y=c(all_HCC44_10nM_AZD1390[3,7],0.0001), lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black", bg="#9E0142", lty=3)
points(x=dose_HCC44, y=all_HCC44_1nM_AZD1390[,7], lwd=3,ylim=c(0.001,10), cex=4, pch=23, col="black", bg="#F46D43", lty=2)
points(x=dose_HCC44, y=all_HCC44_10nM_AZD1390[,7], lwd=3,ylim=c(0.001,10), cex=4, pch=22, col="black", bg="#9E0142", lty=2)
points(x=dose_HCC44, y=all_HCC44_DMSO[,7], lwd=3,ylim=c(0.001,10), cex=4, pch=21, col="black", bg="#66C2A5", lty=2)
minor.tick(nx =2/1,ny=1,tick.ratio =1,  y.args = list( lwd = 3,tcl = -1),  x.args = list( lwd =3,tcl = -1)) 
title( ylab="Surviving Fraction [%]", cex.lab=3.5,cex.main=6,line=9)
title( xlab="Dose [Gy]", cex.lab=3.5,cex.main=6,line=7)
axis(1,lwd=3,at = c(0,4,8),mgp=c(5, 3.5,0),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5)
axis(2,lwd=3,at = c(1000,100,10,1, 0.1, 0.01,0.001,0.0001),label= c("","","","1", "", "",parse(text="10^-3"),""),las=2,mgp=c(7, 2,0),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5)
box(lwd=3,lty="solid")
title(main="HCC44", cex.lab=3,cex.main=3,line=2)
legend(x=c(3.5,8),y=c(2,1000),leg.txt,pch=c(21,22,23), box.col="white",y.intersp = 1,adj = c(0, 0.3), x.intersp = 0.7,ncol=1,pt.bg=c("#66C2A5","#F46D43","#9E0142"),inset=c(-0.0,0.0),pt.cex=3,pt.lwd = 3, cex =2.75,box.lwd=2)

leg.txt <- c("DMSO","1nM AZD1390","10nM AZD1390")
plot(x=dose_H23, y=all_H23_1nM_AZD1390[,9],xlim=c(0,8),mgp=c(6, 0,1), type="l",cex.lab=4,lwd=3,ylim=c(0.0001,1000),axes=FALSE, cex=4, 
     pch=23,log="y", col="black", lty=3, cex.main=4,ann=FALSE)
lines(x=dose_H23, y=all_H23_DMSO[,9], lwd=3,ylim=c(0.0001,750), cex=2, pch=21, col="black",  lty=3)
lines(x=dose_H23[1:4], y=all_H23_10nM_AZD1390[1:4,9], lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black", lty=3)
lines(x=c(2,3), y=c(all_H23_10nM_AZD1390[4,9],0.0001), lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black",lty=3)
lines(x=c(4,5), y=c(all_H23_1nM_AZD1390[5,9],0.0001), lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black", lty=3)
lines(x=c(2,3), y=c(all_H23_DMSO[4,9],0.0001), lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black", lty=3)
points(x=dose_H23, y=all_H23_1nM_AZD1390[,9], lwd=3,ylim=c(0.001,10), cex=4, pch=23, col="black", bg="#F46D43", lty=2)
points(x=dose_H23, y=all_H23_10nM_AZD1390[,9], lwd=3,ylim=c(0.001,10), cex=4, pch=22, col="black", bg="#9E0142", lty=2)
points(x=dose_H23, y=all_H23_DMSO[,9], lwd=3,ylim=c(0.001,10), cex=4, pch=21, col="black", bg="#66C2A5", lty=2)
axis(1,lwd=3,at = c(0,4,8),mgp=c(5, 4,0),lwd.ticks = 3,cex=3, cex.axis=3,tcl = -1.5)
axis(2,lwd=3,at = c(1000,100,10,1, 0.1, 0.01,0.001,0.0001),las=2,label= c("","","","1", "", "",parse(text="10^-3"),""),mgp=c(7, 2,0),lwd.ticks = 3,cex=3.5, cex.axis=3,tcl = -1.5)
box(lwd=3,lty="solid")
title(main="NCI-H23", cex.lab=3,cex.main=3,line=2)
title( ylab="Surviving Fraction [%]", cex.lab=3.5,cex.main=6,line=9)
title( xlab="Dose [Gy]", cex.lab=3.5,cex.main=6,line=7)
minor.tick(nx =2/1,ny=1,tick.ratio =1, y.args = list( lwd = 3,tcl = -1),  x.args = list( lwd = 3,tcl = -1)) 
legend(x=c(3.5,8),y=c(2,1000),leg.txt,pch=c(21,22,23),  box.col="white",y.intersp = 1,adj = c(0, 0.3), x.intersp = 0.7,ncol=1,pt.bg=c("#66C2A5","#F46D43","#9E0142"),inset=c(-0.0,0.0),pt.cex=3,pt.lwd = 3, cex =2.75,box.lwd=2)
dev.off()
setwd(source_wd)