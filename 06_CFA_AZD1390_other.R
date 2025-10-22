#####set directory#####
setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/CFA/Exp4_ATMinhibition/raw")
dose <-c(0,2,4,8)
A549<- read.csv("CFA_AZD1390_A549.csv",sep=",")
HCC827<- read.csv("CFA_AZD1390_HCC827.csv",sep=",")
HCC78<- read.csv("CFA_AZD1390_HCC78.csv",sep=",")
#####Caclulation A549 DMSO I####
A549_DMSO_1 <- A549[19:30,2:4]
colnames(A549_DMSO_1)<- c("Dose","N0","N")
A549_DMSO_1["N0/N1"]<-as.numeric(A549_DMSO_1[,3])/as.numeric(A549_DMSO_1[,2])
A549_DMSO_1["Plating effeciency"]<- (A549_DMSO_1[,4]/mean(as.numeric(A549_DMSO_1[1:3,4])))
cells_A549_DMSO_1 <- data.frame("N0"=as.numeric(A549_DMSO_1[c(1,4,7,10),2]))

A549_DMSO_1_0Gy <- (as.numeric(A549_DMSO_1[1:3,3:3]))
A549_DMSO_1_2Gy <- (as.numeric(A549_DMSO_1[4:6,3:3]))
A549_DMSO_1_4Gy <- (as.numeric(A549_DMSO_1[7:9,3:3]))
A549_DMSO_1_8Gy <- as.numeric(A549_DMSO_1[10:12,3:3])
PE_A549_DMSO_1<- mean(mean(A549_DMSO_1_0Gy)/as.numeric(A549_DMSO_1[1,2]))
mean_A549_DMSO_1 <- data.frame("dose_Gy"=c(0,2,4,8),
                               "N0"=c(cells_A549_DMSO_1),
                               "N_mean"= c(mean(A549_DMSO_1_0Gy,na.rm=TRUE),mean(A549_DMSO_1_2Gy,na.rm=TRUE),mean(A549_DMSO_1_4Gy,na.rm=TRUE),mean(A549_DMSO_1_8Gy,na.rm=TRUE)),
                               "N_SD"= c(sd(A549_DMSO_1_0Gy,na.rm=TRUE),sd(A549_DMSO_1_2Gy,na.rm=TRUE),sd(A549_DMSO_1_4Gy,na.rm=TRUE),sd(A549_DMSO_1_8Gy,na.rm=TRUE)))
mean_A549_DMSO_1["SEM"] <- mean_A549_DMSO_1[,4]/sqrt(3)
mean_A549_DMSO_1["N/N0"] <- mean_A549_DMSO_1[,3]/mean_A549_DMSO_1[,2]
mean_A549_DMSO_1["PE"] <- mean_A549_DMSO_1[1,3]/mean_A549_DMSO_1[1,2]
mean_A549_DMSO_1["SF"] <- mean_A549_DMSO_1["N_mean"]/(mean_A549_DMSO_1["PE"]*mean_A549_DMSO_1["N0"])

A549_DMSO_1["SF"]<- as.numeric(unlist(A549_DMSO_1["N"]))/(mean_A549_DMSO_1[1,7]*as.numeric(unlist(A549_DMSO_1["N0"])))
#####Caclulation A549 10nM AZD1390 I####
A549_10nM_AZD1390_1 <- A549[3:14,2:4]
colnames(A549_10nM_AZD1390_1)<- c("Dose","N0","N")
A549_10nM_AZD1390_1["N0/N1"]<-as.numeric(A549_10nM_AZD1390_1[,3])/as.numeric(A549_10nM_AZD1390_1[,2])
A549_10nM_AZD1390_1["Plating effeciency"]<- (A549_10nM_AZD1390_1[,4]/mean(as.numeric(A549_10nM_AZD1390_1[1:3,4])))
cells_A549_10nM_AZD1390_1 <- data.frame("N0"=as.numeric(A549_10nM_AZD1390_1[c(1,4,7,10),2]))

A549_10nM_AZD1390_1_0Gy <- (as.numeric(A549_10nM_AZD1390_1[1:3,3:3]))
A549_10nM_AZD1390_1_2Gy <- (as.numeric(A549_10nM_AZD1390_1[4:6,3:3]))
A549_10nM_AZD1390_1_4Gy <- (as.numeric(A549_10nM_AZD1390_1[7:9,3:3]))
A549_10nM_AZD1390_1_8Gy <- as.numeric(A549_10nM_AZD1390_1[10:12,3:3])
PE_A549_10nM_AZD1390_1<- mean(mean(A549_10nM_AZD1390_1_0Gy)/as.numeric(A549_10nM_AZD1390_1[1,2]))
mean_A549_10nM_AZD1390_1 <- data.frame("dose_Gy"=c(0,2,4,8),
                                       "N0"=c(cells_A549_10nM_AZD1390_1),
                                       "N_mean"= c(mean(A549_10nM_AZD1390_1_0Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_1_2Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_1_4Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_1_8Gy,na.rm=TRUE)),
                                       "N_SD"= c(sd(A549_10nM_AZD1390_1_0Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_1_2Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_1_4Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_1_8Gy,na.rm=TRUE)))
mean_A549_10nM_AZD1390_1["SEM"] <- mean_A549_10nM_AZD1390_1[,4]/sqrt(3)
mean_A549_10nM_AZD1390_1["N/N0"] <- mean_A549_10nM_AZD1390_1[,3]/mean_A549_10nM_AZD1390_1[,2]
mean_A549_10nM_AZD1390_1["PE"] <- mean_A549_10nM_AZD1390_1[1,3]/mean_A549_10nM_AZD1390_1[1,2]
mean_A549_10nM_AZD1390_1["SF"] <- mean_A549_10nM_AZD1390_1["N_mean"]/(mean_A549_10nM_AZD1390_1["PE"]*mean_A549_10nM_AZD1390_1["N0"])

A549_10nM_AZD1390_1["SF"]<- as.numeric(unlist(A549_10nM_AZD1390_1["N"]))/(mean_A549_10nM_AZD1390_1[1,7]*as.numeric(unlist(A549_10nM_AZD1390_1["N0"])))
#####Caclulation A549 DMSO II####
A549_DMSO_2 <- A549[19:30,7:9]
colnames(A549_DMSO_2)<- c("Dose","N0","N")
A549_DMSO_2["N0/N1"]<-as.numeric(A549_DMSO_2[,3])/as.numeric(A549_DMSO_2[,2])
A549_DMSO_2["Plating effeciency"]<- (A549_DMSO_2[,4]/mean(as.numeric(A549_DMSO_2[1:3,4])))
cells_A549_DMSO_2 <- data.frame("N0"=as.numeric(A549_DMSO_2[c(1,4,7,10),2]))

A549_DMSO_2_0Gy <- (as.numeric(A549_DMSO_2[1:3,3:3]))
A549_DMSO_2_2Gy <- (as.numeric(A549_DMSO_2[4:6,3:3]))
A549_DMSO_2_4Gy <- (as.numeric(A549_DMSO_2[7:9,3:3]))
A549_DMSO_2_8Gy <- as.numeric(A549_DMSO_2[10:12,3:3])
PE_A549_DMSO_2<- mean(mean(A549_DMSO_2_0Gy)/as.numeric(A549_DMSO_2[1,2]))
mean_A549_DMSO_2 <- data.frame("dose_Gy"=c(0,2,4,8),
                               "N0"=c(cells_A549_DMSO_2),
                               "N_mean"= c(mean(A549_DMSO_2_0Gy,na.rm=TRUE),mean(A549_DMSO_2_2Gy,na.rm=TRUE),mean(A549_DMSO_2_4Gy,na.rm=TRUE),mean(A549_DMSO_2_8Gy,na.rm=TRUE)),
                               "N_SD"= c(sd(A549_DMSO_2_0Gy,na.rm=TRUE),sd(A549_DMSO_2_2Gy,na.rm=TRUE),sd(A549_DMSO_2_4Gy,na.rm=TRUE),sd(A549_DMSO_2_8Gy,na.rm=TRUE)))
mean_A549_DMSO_2["SEM"] <- mean_A549_DMSO_2[,4]/sqrt(3)
mean_A549_DMSO_2["N/N0"] <- mean_A549_DMSO_2[,3]/mean_A549_DMSO_2[,2]
mean_A549_DMSO_2["PE"] <- mean_A549_DMSO_2[1,3]/mean_A549_DMSO_2[1,2]
mean_A549_DMSO_2["SF"] <- mean_A549_DMSO_2["N_mean"]/(mean_A549_DMSO_2["PE"]*mean_A549_DMSO_2["N0"])

A549_DMSO_2["SF"]<- as.numeric(unlist(A549_DMSO_2["N"]))/(mean_A549_DMSO_2[1,7]*as.numeric(unlist(A549_DMSO_2["N0"])))
#####Caclulation A549 10nM AZD1390 II####
A549_10nM_AZD1390_2 <- A549[3:14,7:9]
colnames(A549_10nM_AZD1390_2)<- c("Dose","N0","N")
A549_10nM_AZD1390_2["N0/N1"]<-as.numeric(A549_10nM_AZD1390_2[,3])/as.numeric(A549_10nM_AZD1390_2[,2])
A549_10nM_AZD1390_2["Plating effeciency"]<- (A549_10nM_AZD1390_2[,4]/mean(as.numeric(A549_10nM_AZD1390_2[1:3,4])))
cells_A549_10nM_AZD1390_2 <- data.frame("N0"=as.numeric(A549_10nM_AZD1390_2[c(1,4,7,10),2]))

A549_10nM_AZD1390_2_0Gy <- (as.numeric(A549_10nM_AZD1390_2[1:3,3:3]))
A549_10nM_AZD1390_2_2Gy <- (as.numeric(A549_10nM_AZD1390_2[4:6,3:3]))
A549_10nM_AZD1390_2_4Gy <- (as.numeric(A549_10nM_AZD1390_2[7:9,3:3]))
A549_10nM_AZD1390_2_8Gy <- as.numeric(A549_10nM_AZD1390_2[10:12,3:3])
PE_A549_10nM_AZD1390_2<- mean(mean(A549_10nM_AZD1390_2_0Gy)/as.numeric(A549_10nM_AZD1390_2[1,2]))
mean_A549_10nM_AZD1390_2 <- data.frame("dose_Gy"=c(0,2,4,8),
                                       "N0"=c(cells_A549_10nM_AZD1390_2),
                                       "N_mean"= c(mean(A549_10nM_AZD1390_2_0Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_2_2Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_2_4Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_2_8Gy,na.rm=TRUE)),
                                       "N_SD"= c(sd(A549_10nM_AZD1390_2_0Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_2_2Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_2_4Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_2_8Gy,na.rm=TRUE)))
mean_A549_10nM_AZD1390_2["SEM"] <- mean_A549_10nM_AZD1390_2[,4]/sqrt(3)
mean_A549_10nM_AZD1390_2["N/N0"] <- mean_A549_10nM_AZD1390_2[,3]/mean_A549_10nM_AZD1390_2[,2]
mean_A549_10nM_AZD1390_2["PE"] <- mean_A549_10nM_AZD1390_2[1,3]/mean_A549_10nM_AZD1390_2[1,2]
mean_A549_10nM_AZD1390_2["SF"] <- mean_A549_10nM_AZD1390_2["N_mean"]/(mean_A549_10nM_AZD1390_2["PE"]*mean_A549_10nM_AZD1390_2["N0"])

A549_10nM_AZD1390_2["SF"]<- as.numeric(unlist(A549_10nM_AZD1390_2["N"]))/(mean_A549_10nM_AZD1390_2[1,7]*as.numeric(unlist(A549_10nM_AZD1390_2["N0"])))
#####Caclulation A549 DMSO III####
A549_DMSO_3 <- A549[19:30,12:14]
colnames(A549_DMSO_3)<- c("Dose","N0","N")
A549_DMSO_3["N0/N1"]<-as.numeric(A549_DMSO_3[,3])/as.numeric(A549_DMSO_3[,2])
A549_DMSO_3["Plating effeciency"]<- (A549_DMSO_3[,4]/mean(as.numeric(A549_DMSO_3[1:3,4])))
cells_A549_DMSO_3 <- data.frame("N0"=as.numeric(A549_DMSO_3[c(1,4,7,10),2]))

A549_DMSO_3_0Gy <- (as.numeric(A549_DMSO_3[1:3,3:3]))
A549_DMSO_3_2Gy <- (as.numeric(A549_DMSO_3[4:6,3:3]))
A549_DMSO_3_4Gy <- (as.numeric(A549_DMSO_3[7:9,3:3]))
A549_DMSO_3_8Gy <- as.numeric(A549_DMSO_3[10:12,3:3])
PE_A549_DMSO_3<- mean(mean(A549_DMSO_3_0Gy)/as.numeric(A549_DMSO_3[1,2]))
mean_A549_DMSO_3 <- data.frame("dose_Gy"=c(0,2,4,8),
                               "N0"=c(cells_A549_DMSO_3),
                               "N_mean"= c(mean(A549_DMSO_3_0Gy,na.rm=TRUE),mean(A549_DMSO_3_2Gy,na.rm=TRUE),mean(A549_DMSO_3_4Gy,na.rm=TRUE),mean(A549_DMSO_3_8Gy,na.rm=TRUE)),
                               "N_SD"= c(sd(A549_DMSO_3_0Gy,na.rm=TRUE),sd(A549_DMSO_3_2Gy,na.rm=TRUE),sd(A549_DMSO_3_4Gy,na.rm=TRUE),sd(A549_DMSO_3_8Gy,na.rm=TRUE)))
mean_A549_DMSO_3["SEM"] <- mean_A549_DMSO_3[,4]/sqrt(3)
mean_A549_DMSO_3["N/N0"] <- mean_A549_DMSO_3[,3]/mean_A549_DMSO_3[,2]
mean_A549_DMSO_3["PE"] <- mean_A549_DMSO_3[1,3]/mean_A549_DMSO_3[1,2]
mean_A549_DMSO_3["SF"] <- mean_A549_DMSO_3["N_mean"]/(mean_A549_DMSO_3["PE"]*mean_A549_DMSO_3["N0"])

A549_DMSO_3["SF"]<- as.numeric(unlist(A549_DMSO_3["N"]))/(mean_A549_DMSO_3[1,7]*as.numeric(unlist(A549_DMSO_3["N0"])))
#####Caclulation A549 10nM AZD1390 III####
A549_10nM_AZD1390_3 <- A549[3:14,12:14]
colnames(A549_10nM_AZD1390_3)<- c("Dose","N0","N")
A549_10nM_AZD1390_3["N0/N1"]<-as.numeric(A549_10nM_AZD1390_3[,3])/as.numeric(A549_10nM_AZD1390_3[,2])
A549_10nM_AZD1390_3["Plating effeciency"]<- (A549_10nM_AZD1390_3[,4]/mean(as.numeric(A549_10nM_AZD1390_3[1:3,4])))
cells_A549_10nM_AZD1390_3 <- data.frame("N0"=as.numeric(A549_10nM_AZD1390_3[c(1,4,7,10),2]))

A549_10nM_AZD1390_3_0Gy <- (as.numeric(A549_10nM_AZD1390_3[1:3,3:3]))
A549_10nM_AZD1390_3_2Gy <- (as.numeric(A549_10nM_AZD1390_3[4:6,3:3]))
A549_10nM_AZD1390_3_4Gy <- (as.numeric(A549_10nM_AZD1390_3[7:9,3:3]))
A549_10nM_AZD1390_3_8Gy <- as.numeric(A549_10nM_AZD1390_3[10:12,3:3])
PE_A549_10nM_AZD1390_3<- mean(mean(A549_10nM_AZD1390_3_0Gy)/as.numeric(A549_10nM_AZD1390_3[1,2]))
mean_A549_10nM_AZD1390_3 <- data.frame("dose_Gy"=c(0,2,4,8),
                                       "N0"=c(cells_A549_10nM_AZD1390_3),
                                       "N_mean"= c(mean(A549_10nM_AZD1390_3_0Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_3_2Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_3_4Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_3_8Gy,na.rm=TRUE)),
                                       "N_SD"= c(sd(A549_10nM_AZD1390_3_0Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_3_2Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_3_4Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_3_8Gy,na.rm=TRUE)))
mean_A549_10nM_AZD1390_3["SEM"] <- mean_A549_10nM_AZD1390_3[,4]/sqrt(3)
mean_A549_10nM_AZD1390_3["N/N0"] <- mean_A549_10nM_AZD1390_3[,3]/mean_A549_10nM_AZD1390_3[,2]
mean_A549_10nM_AZD1390_3["PE"] <- mean_A549_10nM_AZD1390_3[1,3]/mean_A549_10nM_AZD1390_3[1,2]
mean_A549_10nM_AZD1390_3["SF"] <- mean_A549_10nM_AZD1390_3["N_mean"]/(mean_A549_10nM_AZD1390_3["PE"]*mean_A549_10nM_AZD1390_3["N0"])

A549_10nM_AZD1390_3["SF"]<- as.numeric(unlist(A549_10nM_AZD1390_3["N"]))/(mean_A549_10nM_AZD1390_3[1,7]*as.numeric(unlist(A549_10nM_AZD1390_3["N0"])))
#####Caclulation A549 DMSO IV####
A549_DMSO_4 <- A549[19:30,17:19]
colnames(A549_DMSO_4)<- c("Dose","N0","N")
A549_DMSO_4["N0/N1"]<-as.numeric(A549_DMSO_4[,3])/as.numeric(A549_DMSO_4[,2])
A549_DMSO_4["Plating effeciency"]<- (A549_DMSO_4[,4]/mean(as.numeric(A549_DMSO_4[1:3,4])))
cells_A549_DMSO_4 <- data.frame("N0"=as.numeric(A549_DMSO_4[c(1,4,7,10),2]))

A549_DMSO_4_0Gy <- (as.numeric(A549_DMSO_4[1:3,3:3]))
A549_DMSO_4_2Gy <- (as.numeric(A549_DMSO_4[4:6,3:3]))
A549_DMSO_4_4Gy <- (as.numeric(A549_DMSO_4[7:9,3:3]))
A549_DMSO_4_8Gy <- as.numeric(A549_DMSO_4[10:12,3:3])
PE_A549_DMSO_4<- mean(mean(A549_DMSO_4_0Gy)/as.numeric(A549_DMSO_4[1,2]))
mean_A549_DMSO_4 <- data.frame("dose_Gy"=c(0,2,4,8),
                               "N0"=c(cells_A549_DMSO_4),
                               "N_mean"= c(mean(A549_DMSO_4_0Gy,na.rm=TRUE),mean(A549_DMSO_4_2Gy,na.rm=TRUE),mean(A549_DMSO_4_4Gy,na.rm=TRUE),mean(A549_DMSO_4_8Gy,na.rm=TRUE)),
                               "N_SD"= c(sd(A549_DMSO_4_0Gy,na.rm=TRUE),sd(A549_DMSO_4_2Gy,na.rm=TRUE),sd(A549_DMSO_4_4Gy,na.rm=TRUE),sd(A549_DMSO_4_8Gy,na.rm=TRUE)))
mean_A549_DMSO_4["SEM"] <- mean_A549_DMSO_4[,4]/sqrt(3)
mean_A549_DMSO_4["N/N0"] <- mean_A549_DMSO_4[,3]/mean_A549_DMSO_4[,2]
mean_A549_DMSO_4["PE"] <- mean_A549_DMSO_4[1,3]/mean_A549_DMSO_4[1,2]
mean_A549_DMSO_4["SF"] <- mean_A549_DMSO_4["N_mean"]/(mean_A549_DMSO_4["PE"]*mean_A549_DMSO_4["N0"])

A549_DMSO_4["SF"]<- as.numeric(unlist(A549_DMSO_4["N"]))/(mean_A549_DMSO_4[1,7]*as.numeric(unlist(A549_DMSO_4["N0"])))
#####Caclulation A549 10nM AZD1390 IV####
A549_10nM_AZD1390_4 <- A549[3:14,17:19]
colnames(A549_10nM_AZD1390_4)<- c("Dose","N0","N")
A549_10nM_AZD1390_4["N0/N1"]<-as.numeric(A549_10nM_AZD1390_4[,3])/as.numeric(A549_10nM_AZD1390_4[,2])
A549_10nM_AZD1390_4["Plating effeciency"]<- (A549_10nM_AZD1390_4[,4]/mean(as.numeric(A549_10nM_AZD1390_4[1:3,4])))
cells_A549_10nM_AZD1390_4 <- data.frame("N0"=as.numeric(A549_10nM_AZD1390_4[c(1,4,7,10),2]))

A549_10nM_AZD1390_4_0Gy <- (as.numeric(A549_10nM_AZD1390_4[1:3,3:3]))
A549_10nM_AZD1390_4_2Gy <- (as.numeric(A549_10nM_AZD1390_4[4:6,3:3]))
A549_10nM_AZD1390_4_4Gy <- (as.numeric(A549_10nM_AZD1390_4[7:9,3:3]))
A549_10nM_AZD1390_4_8Gy <- as.numeric(A549_10nM_AZD1390_4[10:12,3:3])
PE_A549_10nM_AZD1390_4<- mean(mean(A549_10nM_AZD1390_4_0Gy)/as.numeric(A549_10nM_AZD1390_4[1,2]))
mean_A549_10nM_AZD1390_4 <- data.frame("dose_Gy"=c(0,2,4,8),
                                       "N0"=c(cells_A549_10nM_AZD1390_4),
                                       "N_mean"= c(mean(A549_10nM_AZD1390_4_0Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_4_2Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_4_4Gy,na.rm=TRUE),mean(A549_10nM_AZD1390_4_8Gy,na.rm=TRUE)),
                                       "N_SD"= c(sd(A549_10nM_AZD1390_4_0Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_4_2Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_4_4Gy,na.rm=TRUE),sd(A549_10nM_AZD1390_4_8Gy,na.rm=TRUE)))
mean_A549_10nM_AZD1390_4["SEM"] <- mean_A549_10nM_AZD1390_4[,4]/sqrt(3)
mean_A549_10nM_AZD1390_4["N/N0"] <- mean_A549_10nM_AZD1390_4[,3]/mean_A549_10nM_AZD1390_4[,2]
mean_A549_10nM_AZD1390_4["PE"] <- mean_A549_10nM_AZD1390_4[1,3]/mean_A549_10nM_AZD1390_4[1,2]
mean_A549_10nM_AZD1390_4["SF"] <- mean_A549_10nM_AZD1390_4["N_mean"]/(mean_A549_10nM_AZD1390_4["PE"]*mean_A549_10nM_AZD1390_4["N0"])

A549_10nM_AZD1390_4["SF"]<- as.numeric(unlist(A549_10nM_AZD1390_4["N"]))/(mean_A549_10nM_AZD1390_4[1,7]*as.numeric(unlist(A549_10nM_AZD1390_4["N0"])))
#####binding A549 results####
all_A549_DMSO <- cbind(mean_A549_DMSO_1[,c(1:3,8)],mean_A549_DMSO_2[,c(3,8)],mean_A549_DMSO_3[,c(3,8)],mean_A549_DMSO_4[,c(3,8)])
colnames(all_A549_DMSO)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_A549_DMSO["SF_mean"]<- rowMeans(all_A549_DMSO[,c(4,6,8,10)])
all_A549_DMSO["SF_sd"]<- apply(all_A549_DMSO[,c(4,6,8,10)], 1, sd)
all_A549_DMSO["SF_SEM"]<- all_A549_DMSO["SF_sd"]/sqrt(4)

all_A549_10nM_AZD1390 <- cbind(mean_A549_10nM_AZD1390_1[,c(1:3,8)],mean_A549_10nM_AZD1390_2[,c(3,8)],mean_A549_10nM_AZD1390_3[,c(3,8)],mean_A549_10nM_AZD1390_4[,c(3,8)])
colnames(all_A549_10nM_AZD1390)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_A549_10nM_AZD1390["SF_mean"]<- rowMeans(all_A549_10nM_AZD1390[,c(4,6,8,10)])
all_A549_10nM_AZD1390["SF_sd"]<- apply(all_A549_10nM_AZD1390[,c(4,6,8,10)], 1, sd)
all_A549_10nM_AZD1390["SF_SEM"]<- all_A549_10nM_AZD1390["SF_sd"]/sqrt(4)

#####Caclulation HCC827 DMSO I####
HCC827_DMSO_1 <- HCC827[19:30,2:4]
colnames(HCC827_DMSO_1)<- c("Dose","N0","N")
HCC827_DMSO_1["N0/N1"]<-as.numeric(HCC827_DMSO_1[,3])/as.numeric(HCC827_DMSO_1[,2])
HCC827_DMSO_1["Plating effeciency"]<- (HCC827_DMSO_1[,4]/mean(as.numeric(HCC827_DMSO_1[1:3,4])))
cells_HCC827_DMSO_1 <- data.frame("N0"=as.numeric(HCC827_DMSO_1[c(1,4,7,10),2]))

HCC827_DMSO_1_0Gy <- (as.numeric(HCC827_DMSO_1[1:3,3:3]))
HCC827_DMSO_1_2Gy <- (as.numeric(HCC827_DMSO_1[4:6,3:3]))
HCC827_DMSO_1_4Gy <- (as.numeric(HCC827_DMSO_1[7:9,3:3]))
HCC827_DMSO_1_8Gy <- as.numeric(HCC827_DMSO_1[10:12,3:3])
PE_HCC827_DMSO_1<- mean(mean(HCC827_DMSO_1_0Gy)/as.numeric(HCC827_DMSO_1[1,2]))
mean_HCC827_DMSO_1 <- data.frame("dose_Gy"=c(0,2,4,8),
                                 "N0"=c(cells_HCC827_DMSO_1),
                                 "N_mean"= c(mean(HCC827_DMSO_1_0Gy,na.rm=TRUE),mean(HCC827_DMSO_1_2Gy,na.rm=TRUE),mean(HCC827_DMSO_1_4Gy,na.rm=TRUE),mean(HCC827_DMSO_1_8Gy,na.rm=TRUE)),
                                 "N_SD"= c(sd(HCC827_DMSO_1_0Gy,na.rm=TRUE),sd(HCC827_DMSO_1_2Gy,na.rm=TRUE),sd(HCC827_DMSO_1_4Gy,na.rm=TRUE),sd(HCC827_DMSO_1_8Gy,na.rm=TRUE)))
mean_HCC827_DMSO_1["SEM"] <- mean_HCC827_DMSO_1[,4]/sqrt(3)
mean_HCC827_DMSO_1["N/N0"] <- mean_HCC827_DMSO_1[,3]/mean_HCC827_DMSO_1[,2]
mean_HCC827_DMSO_1["PE"] <- mean_HCC827_DMSO_1[1,3]/mean_HCC827_DMSO_1[1,2]
mean_HCC827_DMSO_1["SF"] <- mean_HCC827_DMSO_1["N_mean"]/(mean_HCC827_DMSO_1["PE"]*mean_HCC827_DMSO_1["N0"])

HCC827_DMSO_1["SF"]<- as.numeric(unlist(HCC827_DMSO_1["N"]))/(mean_HCC827_DMSO_1[1,7]*as.numeric(unlist(HCC827_DMSO_1["N0"])))
#####Caclulation HCC827 10nM AZD1390 I####
HCC827_10nM_AZD1390_1 <- HCC827[3:14,2:4]
colnames(HCC827_10nM_AZD1390_1)<- c("Dose","N0","N")
HCC827_10nM_AZD1390_1["N0/N1"]<-as.numeric(HCC827_10nM_AZD1390_1[,3])/as.numeric(HCC827_10nM_AZD1390_1[,2])
HCC827_10nM_AZD1390_1["Plating effeciency"]<- (HCC827_10nM_AZD1390_1[,4]/mean(as.numeric(HCC827_10nM_AZD1390_1[1:3,4])))
cells_HCC827_10nM_AZD1390_1 <- data.frame("N0"=as.numeric(HCC827_10nM_AZD1390_1[c(1,4,7,10),2]))

HCC827_10nM_AZD1390_1_0Gy <- (as.numeric(HCC827_10nM_AZD1390_1[1:3,3:3]))
HCC827_10nM_AZD1390_1_2Gy <- (as.numeric(HCC827_10nM_AZD1390_1[4:6,3:3]))
HCC827_10nM_AZD1390_1_4Gy <- (as.numeric(HCC827_10nM_AZD1390_1[7:9,3:3]))
HCC827_10nM_AZD1390_1_8Gy <- as.numeric(HCC827_10nM_AZD1390_1[10:12,3:3])
PE_HCC827_10nM_AZD1390_1<- mean(mean(HCC827_10nM_AZD1390_1_0Gy)/as.numeric(HCC827_10nM_AZD1390_1[1,2]))
mean_HCC827_10nM_AZD1390_1 <- data.frame("dose_Gy"=c(0,2,4,8),
                                         "N0"=c(cells_HCC827_10nM_AZD1390_1),
                                         "N_mean"= c(mean(HCC827_10nM_AZD1390_1_0Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_1_2Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_1_4Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_1_8Gy,na.rm=TRUE)),
                                         "N_SD"= c(sd(HCC827_10nM_AZD1390_1_0Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_1_2Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_1_4Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_1_8Gy,na.rm=TRUE)))
mean_HCC827_10nM_AZD1390_1["SEM"] <- mean_HCC827_10nM_AZD1390_1[,4]/sqrt(3)
mean_HCC827_10nM_AZD1390_1["N/N0"] <- mean_HCC827_10nM_AZD1390_1[,3]/mean_HCC827_10nM_AZD1390_1[,2]
mean_HCC827_10nM_AZD1390_1["PE"] <- mean_HCC827_10nM_AZD1390_1[1,3]/mean_HCC827_10nM_AZD1390_1[1,2]
mean_HCC827_10nM_AZD1390_1["SF"] <- mean_HCC827_10nM_AZD1390_1["N_mean"]/(mean_HCC827_10nM_AZD1390_1["PE"]*mean_HCC827_10nM_AZD1390_1["N0"])

HCC827_10nM_AZD1390_1["SF"]<- as.numeric(unlist(HCC827_10nM_AZD1390_1["N"]))/(mean_HCC827_10nM_AZD1390_1[1,7]*as.numeric(unlist(HCC827_10nM_AZD1390_1["N0"])))
#####Caclulation HCC827 DMSO II####
HCC827_DMSO_2 <- HCC827[19:30,7:9]
colnames(HCC827_DMSO_2)<- c("Dose","N0","N")
HCC827_DMSO_2["N0/N1"]<-as.numeric(HCC827_DMSO_2[,3])/as.numeric(HCC827_DMSO_2[,2])
HCC827_DMSO_2["Plating effeciency"]<- (HCC827_DMSO_2[,4]/mean(as.numeric(HCC827_DMSO_2[1:3,4])))
cells_HCC827_DMSO_2 <- data.frame("N0"=as.numeric(HCC827_DMSO_2[c(1,4,7,10),2]))

HCC827_DMSO_2_0Gy <- (as.numeric(HCC827_DMSO_2[1:3,3:3]))
HCC827_DMSO_2_2Gy <- (as.numeric(HCC827_DMSO_2[4:6,3:3]))
HCC827_DMSO_2_4Gy <- (as.numeric(HCC827_DMSO_2[7:9,3:3]))
HCC827_DMSO_2_8Gy <- as.numeric(HCC827_DMSO_2[10:12,3:3])
PE_HCC827_DMSO_2<- mean(mean(HCC827_DMSO_2_0Gy)/as.numeric(HCC827_DMSO_2[1,2]))
mean_HCC827_DMSO_2 <- data.frame("dose_Gy"=c(0,2,4,8),
                                 "N0"=c(cells_HCC827_DMSO_2),
                                 "N_mean"= c(mean(HCC827_DMSO_2_0Gy,na.rm=TRUE),mean(HCC827_DMSO_2_2Gy,na.rm=TRUE),mean(HCC827_DMSO_2_4Gy,na.rm=TRUE),mean(HCC827_DMSO_2_8Gy,na.rm=TRUE)),
                                 "N_SD"= c(sd(HCC827_DMSO_2_0Gy,na.rm=TRUE),sd(HCC827_DMSO_2_2Gy,na.rm=TRUE),sd(HCC827_DMSO_2_4Gy,na.rm=TRUE),sd(HCC827_DMSO_2_8Gy,na.rm=TRUE)))
mean_HCC827_DMSO_2["SEM"] <- mean_HCC827_DMSO_2[,4]/sqrt(3)
mean_HCC827_DMSO_2["N/N0"] <- mean_HCC827_DMSO_2[,3]/mean_HCC827_DMSO_2[,2]
mean_HCC827_DMSO_2["PE"] <- mean_HCC827_DMSO_2[1,3]/mean_HCC827_DMSO_2[1,2]
mean_HCC827_DMSO_2["SF"] <- mean_HCC827_DMSO_2["N_mean"]/(mean_HCC827_DMSO_2["PE"]*mean_HCC827_DMSO_2["N0"])

HCC827_DMSO_2["SF"]<- as.numeric(unlist(HCC827_DMSO_2["N"]))/(mean_HCC827_DMSO_2[1,7]*as.numeric(unlist(HCC827_DMSO_2["N0"])))
#####Caclulation HCC827 10nM AZD1390 II####
HCC827_10nM_AZD1390_2 <- HCC827[3:14,7:9]
colnames(HCC827_10nM_AZD1390_2)<- c("Dose","N0","N")
HCC827_10nM_AZD1390_2["N0/N1"]<-as.numeric(HCC827_10nM_AZD1390_2[,3])/as.numeric(HCC827_10nM_AZD1390_2[,2])
HCC827_10nM_AZD1390_2["Plating effeciency"]<- (HCC827_10nM_AZD1390_2[,4]/mean(as.numeric(HCC827_10nM_AZD1390_2[1:3,4])))
cells_HCC827_10nM_AZD1390_2 <- data.frame("N0"=as.numeric(HCC827_10nM_AZD1390_2[c(1,4,7,10),2]))

HCC827_10nM_AZD1390_2_0Gy <- (as.numeric(HCC827_10nM_AZD1390_2[1:3,3:3]))
HCC827_10nM_AZD1390_2_2Gy <- (as.numeric(HCC827_10nM_AZD1390_2[4:6,3:3]))
HCC827_10nM_AZD1390_2_4Gy <- (as.numeric(HCC827_10nM_AZD1390_2[7:9,3:3]))
HCC827_10nM_AZD1390_2_8Gy <- as.numeric(HCC827_10nM_AZD1390_2[10:12,3:3])
PE_HCC827_10nM_AZD1390_2<- mean(mean(HCC827_10nM_AZD1390_2_0Gy)/as.numeric(HCC827_10nM_AZD1390_2[1,2]))
mean_HCC827_10nM_AZD1390_2 <- data.frame("dose_Gy"=c(0,2,4,8),
                                         "N0"=c(cells_HCC827_10nM_AZD1390_2),
                                         "N_mean"= c(mean(HCC827_10nM_AZD1390_2_0Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_2_2Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_2_4Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_2_8Gy,na.rm=TRUE)),
                                         "N_SD"= c(sd(HCC827_10nM_AZD1390_2_0Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_2_2Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_2_4Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_2_8Gy,na.rm=TRUE)))
mean_HCC827_10nM_AZD1390_2["SEM"] <- mean_HCC827_10nM_AZD1390_2[,4]/sqrt(3)
mean_HCC827_10nM_AZD1390_2["N/N0"] <- mean_HCC827_10nM_AZD1390_2[,3]/mean_HCC827_10nM_AZD1390_2[,2]
mean_HCC827_10nM_AZD1390_2["PE"] <- mean_HCC827_10nM_AZD1390_2[1,3]/mean_HCC827_10nM_AZD1390_2[1,2]
mean_HCC827_10nM_AZD1390_2["SF"] <- mean_HCC827_10nM_AZD1390_2["N_mean"]/(mean_HCC827_10nM_AZD1390_2["PE"]*mean_HCC827_10nM_AZD1390_2["N0"])

HCC827_10nM_AZD1390_2["SF"]<- as.numeric(unlist(HCC827_10nM_AZD1390_2["N"]))/(mean_HCC827_10nM_AZD1390_2[1,7]*as.numeric(unlist(HCC827_10nM_AZD1390_2["N0"])))
#####Caclulation HCC827 DMSO III####
HCC827_DMSO_3 <- HCC827[19:30,12:14]
colnames(HCC827_DMSO_3)<- c("Dose","N0","N")
HCC827_DMSO_3["N0/N1"]<-as.numeric(HCC827_DMSO_3[,3])/as.numeric(HCC827_DMSO_3[,2])
HCC827_DMSO_3["Plating effeciency"]<- (HCC827_DMSO_3[,4]/mean(as.numeric(HCC827_DMSO_3[1:3,4])))
cells_HCC827_DMSO_3 <- data.frame("N0"=as.numeric(HCC827_DMSO_3[c(1,4,7,10),2]))

HCC827_DMSO_3_0Gy <- (as.numeric(HCC827_DMSO_3[1:3,3:3]))
HCC827_DMSO_3_2Gy <- (as.numeric(HCC827_DMSO_3[4:6,3:3]))
HCC827_DMSO_3_4Gy <- (as.numeric(HCC827_DMSO_3[7:9,3:3]))
HCC827_DMSO_3_8Gy <- as.numeric(HCC827_DMSO_3[10:12,3:3])
PE_HCC827_DMSO_3<- mean(mean(HCC827_DMSO_3_0Gy)/as.numeric(HCC827_DMSO_3[1,2]))
mean_HCC827_DMSO_3 <- data.frame("dose_Gy"=c(0,2,4,8),
                                 "N0"=c(cells_HCC827_DMSO_3),
                                 "N_mean"= c(mean(HCC827_DMSO_3_0Gy,na.rm=TRUE),mean(HCC827_DMSO_3_2Gy,na.rm=TRUE),mean(HCC827_DMSO_3_4Gy,na.rm=TRUE),mean(HCC827_DMSO_3_8Gy,na.rm=TRUE)),
                                 "N_SD"= c(sd(HCC827_DMSO_3_0Gy,na.rm=TRUE),sd(HCC827_DMSO_3_2Gy,na.rm=TRUE),sd(HCC827_DMSO_3_4Gy,na.rm=TRUE),sd(HCC827_DMSO_3_8Gy,na.rm=TRUE)))
mean_HCC827_DMSO_3["SEM"] <- mean_HCC827_DMSO_3[,4]/sqrt(3)
mean_HCC827_DMSO_3["N/N0"] <- mean_HCC827_DMSO_3[,3]/mean_HCC827_DMSO_3[,2]
mean_HCC827_DMSO_3["PE"] <- mean_HCC827_DMSO_3[1,3]/mean_HCC827_DMSO_3[1,2]
mean_HCC827_DMSO_3["SF"] <- mean_HCC827_DMSO_3["N_mean"]/(mean_HCC827_DMSO_3["PE"]*mean_HCC827_DMSO_3["N0"])

HCC827_DMSO_3["SF"]<- as.numeric(unlist(HCC827_DMSO_3["N"]))/(mean_HCC827_DMSO_3[1,7]*as.numeric(unlist(HCC827_DMSO_3["N0"])))
#####Caclulation HCC827 10nM AZD1390 III####
HCC827_10nM_AZD1390_3 <- HCC827[3:14,12:14]
colnames(HCC827_10nM_AZD1390_3)<- c("Dose","N0","N")
HCC827_10nM_AZD1390_3["N0/N1"]<-as.numeric(HCC827_10nM_AZD1390_3[,3])/as.numeric(HCC827_10nM_AZD1390_3[,2])
HCC827_10nM_AZD1390_3["Plating effeciency"]<- (HCC827_10nM_AZD1390_3[,4]/mean(as.numeric(HCC827_10nM_AZD1390_3[1:3,4])))
cells_HCC827_10nM_AZD1390_3 <- data.frame("N0"=as.numeric(HCC827_10nM_AZD1390_3[c(1,4,7,10),2]))

HCC827_10nM_AZD1390_3_0Gy <- (as.numeric(HCC827_10nM_AZD1390_3[1:3,3:3]))
HCC827_10nM_AZD1390_3_2Gy <- (as.numeric(HCC827_10nM_AZD1390_3[4:6,3:3]))
HCC827_10nM_AZD1390_3_4Gy <- (as.numeric(HCC827_10nM_AZD1390_3[7:9,3:3]))
HCC827_10nM_AZD1390_3_8Gy <- as.numeric(HCC827_10nM_AZD1390_3[10:12,3:3])
PE_HCC827_10nM_AZD1390_3<- mean(mean(HCC827_10nM_AZD1390_3_0Gy)/as.numeric(HCC827_10nM_AZD1390_3[1,2]))
mean_HCC827_10nM_AZD1390_3 <- data.frame("dose_Gy"=c(0,2,4,8),
                                         "N0"=c(cells_HCC827_10nM_AZD1390_3),
                                         "N_mean"= c(mean(HCC827_10nM_AZD1390_3_0Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_3_2Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_3_4Gy,na.rm=TRUE),mean(HCC827_10nM_AZD1390_3_8Gy,na.rm=TRUE)),
                                         "N_SD"= c(sd(HCC827_10nM_AZD1390_3_0Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_3_2Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_3_4Gy,na.rm=TRUE),sd(HCC827_10nM_AZD1390_3_8Gy,na.rm=TRUE)))
mean_HCC827_10nM_AZD1390_3["SEM"] <- mean_HCC827_10nM_AZD1390_3[,4]/sqrt(3)
mean_HCC827_10nM_AZD1390_3["N/N0"] <- mean_HCC827_10nM_AZD1390_3[,3]/mean_HCC827_10nM_AZD1390_3[,2]
mean_HCC827_10nM_AZD1390_3["PE"] <- mean_HCC827_10nM_AZD1390_3[1,3]/mean_HCC827_10nM_AZD1390_3[1,2]
mean_HCC827_10nM_AZD1390_3["SF"] <- mean_HCC827_10nM_AZD1390_3["N_mean"]/(mean_HCC827_10nM_AZD1390_3["PE"]*mean_HCC827_10nM_AZD1390_3["N0"])

HCC827_10nM_AZD1390_3["SF"]<- as.numeric(unlist(HCC827_10nM_AZD1390_3["N"]))/(mean_HCC827_10nM_AZD1390_3[1,7]*as.numeric(unlist(HCC827_10nM_AZD1390_3["N0"])))
#####binding HCC827 results####
all_HCC827_DMSO <- cbind(mean_HCC827_DMSO_1[,c(1:3,8)],mean_HCC827_DMSO_2[,c(3,8)],mean_HCC827_DMSO_3[,c(3,8)])
colnames(all_HCC827_DMSO)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3")
all_HCC827_DMSO["SF_mean"]<- rowMeans(all_HCC827_DMSO[,c(4,6,8)])
all_HCC827_DMSO["SF_sd"]<- apply(all_HCC827_DMSO[,c(4,6,8)], 1, sd)
all_HCC827_DMSO["SF_SEM"]<- all_HCC827_DMSO["SF_sd"]/sqrt(3)

all_HCC827_10nM_AZD1390 <- cbind(mean_HCC827_10nM_AZD1390_1[,c(1:3,8)],mean_HCC827_10nM_AZD1390_2[,c(3,8)],mean_HCC827_10nM_AZD1390_3[,c(3,8)])
colnames(all_HCC827_10nM_AZD1390)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3")
all_HCC827_10nM_AZD1390["SF_mean"]<- rowMeans(all_HCC827_10nM_AZD1390[,c(4,6,8)])
all_HCC827_10nM_AZD1390["SF_sd"]<- apply(all_HCC827_10nM_AZD1390[,c(4,6,8)], 1, sd)
all_HCC827_10nM_AZD1390["SF_SEM"]<- all_HCC827_10nM_AZD1390["SF_sd"]/sqrt(3)
#####Caclulation HCC78 DMSO I####
HCC78_DMSO_1 <- HCC78[19:30,2:4]
colnames(HCC78_DMSO_1)<- c("Dose","N0","N")
HCC78_DMSO_1["N0/N1"]<-as.numeric(HCC78_DMSO_1[,3])/as.numeric(HCC78_DMSO_1[,2])
HCC78_DMSO_1["Plating effeciency"]<- (HCC78_DMSO_1[,4]/mean(as.numeric(HCC78_DMSO_1[1:3,4])))
cells_HCC78_DMSO_1 <- data.frame("N0"=as.numeric(HCC78_DMSO_1[c(1,4,7,10),2]))

HCC78_DMSO_1_0Gy <- (as.numeric(HCC78_DMSO_1[1:3,3:3]))
HCC78_DMSO_1_2Gy <- (as.numeric(HCC78_DMSO_1[4:6,3:3]))
HCC78_DMSO_1_4Gy <- (as.numeric(HCC78_DMSO_1[7:9,3:3]))
HCC78_DMSO_1_8Gy <- as.numeric(HCC78_DMSO_1[10:12,3:3])
PE_HCC78_DMSO_1<- mean(mean(HCC78_DMSO_1_0Gy)/as.numeric(HCC78_DMSO_1[1,2]))
mean_HCC78_DMSO_1 <- data.frame("dose_Gy"=c(0,2,4,8),
                                "N0"=c(cells_HCC78_DMSO_1),
                                "N_mean"= c(mean(HCC78_DMSO_1_0Gy,na.rm=TRUE),mean(HCC78_DMSO_1_2Gy,na.rm=TRUE),mean(HCC78_DMSO_1_4Gy,na.rm=TRUE),mean(HCC78_DMSO_1_8Gy,na.rm=TRUE)),
                                "N_SD"= c(sd(HCC78_DMSO_1_0Gy,na.rm=TRUE),sd(HCC78_DMSO_1_2Gy,na.rm=TRUE),sd(HCC78_DMSO_1_4Gy,na.rm=TRUE),sd(HCC78_DMSO_1_8Gy,na.rm=TRUE)))
mean_HCC78_DMSO_1["SEM"] <- mean_HCC78_DMSO_1[,4]/sqrt(3)
mean_HCC78_DMSO_1["N/N0"] <- mean_HCC78_DMSO_1[,3]/mean_HCC78_DMSO_1[,2]
mean_HCC78_DMSO_1["PE"] <- mean_HCC78_DMSO_1[1,3]/mean_HCC78_DMSO_1[1,2]
mean_HCC78_DMSO_1["SF"] <- mean_HCC78_DMSO_1["N_mean"]/(mean_HCC78_DMSO_1["PE"]*mean_HCC78_DMSO_1["N0"])

HCC78_DMSO_1["SF"]<- as.numeric(unlist(HCC78_DMSO_1["N"]))/(mean_HCC78_DMSO_1[1,7]*as.numeric(unlist(HCC78_DMSO_1["N0"])))
#####Caclulation HCC78 10nM AZD1390 I####
HCC78_10nM_AZD1390_1 <- HCC78[3:14,2:4]
colnames(HCC78_10nM_AZD1390_1)<- c("Dose","N0","N")
HCC78_10nM_AZD1390_1["N0/N1"]<-as.numeric(HCC78_10nM_AZD1390_1[,3])/as.numeric(HCC78_10nM_AZD1390_1[,2])
HCC78_10nM_AZD1390_1["Plating effeciency"]<- (HCC78_10nM_AZD1390_1[,4]/mean(as.numeric(HCC78_10nM_AZD1390_1[1:3,4])))
cells_HCC78_10nM_AZD1390_1 <- data.frame("N0"=as.numeric(HCC78_10nM_AZD1390_1[c(1,4,7,10),2]))

HCC78_10nM_AZD1390_1_0Gy <- (as.numeric(HCC78_10nM_AZD1390_1[1:3,3:3]))
HCC78_10nM_AZD1390_1_2Gy <- (as.numeric(HCC78_10nM_AZD1390_1[4:6,3:3]))
HCC78_10nM_AZD1390_1_4Gy <- (as.numeric(HCC78_10nM_AZD1390_1[7:9,3:3]))
HCC78_10nM_AZD1390_1_8Gy <- as.numeric(HCC78_10nM_AZD1390_1[10:12,3:3])
PE_HCC78_10nM_AZD1390_1<- mean(mean(HCC78_10nM_AZD1390_1_0Gy)/as.numeric(HCC78_10nM_AZD1390_1[1,2]))
mean_HCC78_10nM_AZD1390_1 <- data.frame("dose_Gy"=c(0,2,4,8),
                                        "N0"=c(cells_HCC78_10nM_AZD1390_1),
                                        "N_mean"= c(mean(HCC78_10nM_AZD1390_1_0Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_1_2Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_1_4Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_1_8Gy,na.rm=TRUE)),
                                        "N_SD"= c(sd(HCC78_10nM_AZD1390_1_0Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_1_2Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_1_4Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_1_8Gy,na.rm=TRUE)))
mean_HCC78_10nM_AZD1390_1["SEM"] <- mean_HCC78_10nM_AZD1390_1[,4]/sqrt(3)
mean_HCC78_10nM_AZD1390_1["N/N0"] <- mean_HCC78_10nM_AZD1390_1[,3]/mean_HCC78_10nM_AZD1390_1[,2]
mean_HCC78_10nM_AZD1390_1["PE"] <- mean_HCC78_10nM_AZD1390_1[1,3]/mean_HCC78_10nM_AZD1390_1[1,2]
mean_HCC78_10nM_AZD1390_1["SF"] <- mean_HCC78_10nM_AZD1390_1["N_mean"]/(mean_HCC78_10nM_AZD1390_1["PE"]*mean_HCC78_10nM_AZD1390_1["N0"])

HCC78_10nM_AZD1390_1["SF"]<- as.numeric(unlist(HCC78_10nM_AZD1390_1["N"]))/(mean_HCC78_10nM_AZD1390_1[1,7]*as.numeric(unlist(HCC78_10nM_AZD1390_1["N0"])))
#####Caclulation HCC78 DMSO II####
HCC78_DMSO_2 <- HCC78[19:30,7:9]
colnames(HCC78_DMSO_2)<- c("Dose","N0","N")
HCC78_DMSO_2["N0/N1"]<-as.numeric(HCC78_DMSO_2[,3])/as.numeric(HCC78_DMSO_2[,2])
HCC78_DMSO_2["Plating effeciency"]<- (HCC78_DMSO_2[,4]/mean(as.numeric(HCC78_DMSO_2[1:3,4])))
cells_HCC78_DMSO_2 <- data.frame("N0"=as.numeric(HCC78_DMSO_2[c(1,4,7,10),2]))

HCC78_DMSO_2_0Gy <- (as.numeric(HCC78_DMSO_2[1:3,3:3]))
HCC78_DMSO_2_2Gy <- (as.numeric(HCC78_DMSO_2[4:6,3:3]))
HCC78_DMSO_2_4Gy <- (as.numeric(HCC78_DMSO_2[7:9,3:3]))
HCC78_DMSO_2_8Gy <- as.numeric(HCC78_DMSO_2[10:12,3:3])
PE_HCC78_DMSO_2<- mean(mean(HCC78_DMSO_2_0Gy)/as.numeric(HCC78_DMSO_2[1,2]))
mean_HCC78_DMSO_2 <- data.frame("dose_Gy"=c(0,2,4,8),
                                "N0"=c(cells_HCC78_DMSO_2),
                                "N_mean"= c(mean(HCC78_DMSO_2_0Gy,na.rm=TRUE),mean(HCC78_DMSO_2_2Gy,na.rm=TRUE),mean(HCC78_DMSO_2_4Gy,na.rm=TRUE),mean(HCC78_DMSO_2_8Gy,na.rm=TRUE)),
                                "N_SD"= c(sd(HCC78_DMSO_2_0Gy,na.rm=TRUE),sd(HCC78_DMSO_2_2Gy,na.rm=TRUE),sd(HCC78_DMSO_2_4Gy,na.rm=TRUE),sd(HCC78_DMSO_2_8Gy,na.rm=TRUE)))
mean_HCC78_DMSO_2["SEM"] <- mean_HCC78_DMSO_2[,4]/sqrt(3)
mean_HCC78_DMSO_2["N/N0"] <- mean_HCC78_DMSO_2[,3]/mean_HCC78_DMSO_2[,2]
mean_HCC78_DMSO_2["PE"] <- mean_HCC78_DMSO_2[1,3]/mean_HCC78_DMSO_2[1,2]
mean_HCC78_DMSO_2["SF"] <- mean_HCC78_DMSO_2["N_mean"]/(mean_HCC78_DMSO_2["PE"]*mean_HCC78_DMSO_2["N0"])

HCC78_DMSO_2["SF"]<- as.numeric(unlist(HCC78_DMSO_2["N"]))/(mean_HCC78_DMSO_2[1,7]*as.numeric(unlist(HCC78_DMSO_2["N0"])))
#####Caclulation HCC78 10nM AZD1390 II####
HCC78_10nM_AZD1390_2 <- HCC78[3:14,7:9]
colnames(HCC78_10nM_AZD1390_2)<- c("Dose","N0","N")
HCC78_10nM_AZD1390_2["N0/N1"]<-as.numeric(HCC78_10nM_AZD1390_2[,3])/as.numeric(HCC78_10nM_AZD1390_2[,2])
HCC78_10nM_AZD1390_2["Plating effeciency"]<- (HCC78_10nM_AZD1390_2[,4]/mean(as.numeric(HCC78_10nM_AZD1390_2[1:3,4])))
cells_HCC78_10nM_AZD1390_2 <- data.frame("N0"=as.numeric(HCC78_10nM_AZD1390_2[c(1,4,7,10),2]))

HCC78_10nM_AZD1390_2_0Gy <- (as.numeric(HCC78_10nM_AZD1390_2[1:3,3:3]))
HCC78_10nM_AZD1390_2_2Gy <- (as.numeric(HCC78_10nM_AZD1390_2[4:6,3:3]))
HCC78_10nM_AZD1390_2_4Gy <- (as.numeric(HCC78_10nM_AZD1390_2[7:9,3:3]))
HCC78_10nM_AZD1390_2_8Gy <- as.numeric(HCC78_10nM_AZD1390_2[10:12,3:3])
PE_HCC78_10nM_AZD1390_2<- mean(mean(HCC78_10nM_AZD1390_2_0Gy)/as.numeric(HCC78_10nM_AZD1390_2[1,2]))
mean_HCC78_10nM_AZD1390_2 <- data.frame("dose_Gy"=c(0,2,4,8),
                                        "N0"=c(cells_HCC78_10nM_AZD1390_2),
                                        "N_mean"= c(mean(HCC78_10nM_AZD1390_2_0Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_2_2Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_2_4Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_2_8Gy,na.rm=TRUE)),
                                        "N_SD"= c(sd(HCC78_10nM_AZD1390_2_0Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_2_2Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_2_4Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_2_8Gy,na.rm=TRUE)))
mean_HCC78_10nM_AZD1390_2["SEM"] <- mean_HCC78_10nM_AZD1390_2[,4]/sqrt(3)
mean_HCC78_10nM_AZD1390_2["N/N0"] <- mean_HCC78_10nM_AZD1390_2[,3]/mean_HCC78_10nM_AZD1390_2[,2]
mean_HCC78_10nM_AZD1390_2["PE"] <- mean_HCC78_10nM_AZD1390_2[1,3]/mean_HCC78_10nM_AZD1390_2[1,2]
mean_HCC78_10nM_AZD1390_2["SF"] <- mean_HCC78_10nM_AZD1390_2["N_mean"]/(mean_HCC78_10nM_AZD1390_2["PE"]*mean_HCC78_10nM_AZD1390_2["N0"])

HCC78_10nM_AZD1390_2["SF"]<- as.numeric(unlist(HCC78_10nM_AZD1390_2["N"]))/(mean_HCC78_10nM_AZD1390_2[1,7]*as.numeric(unlist(HCC78_10nM_AZD1390_2["N0"])))
#####Caclulation HCC78 DMSO III####
HCC78_DMSO_3 <- HCC78[19:30,12:14]
colnames(HCC78_DMSO_3)<- c("Dose","N0","N")
HCC78_DMSO_3["N0/N1"]<-as.numeric(HCC78_DMSO_3[,3])/as.numeric(HCC78_DMSO_3[,2])
HCC78_DMSO_3["Plating effeciency"]<- (HCC78_DMSO_3[,4]/mean(as.numeric(HCC78_DMSO_3[1:3,4])))
cells_HCC78_DMSO_3 <- data.frame("N0"=as.numeric(HCC78_DMSO_3[c(1,4,7,10),2]))

HCC78_DMSO_3_0Gy <- (as.numeric(HCC78_DMSO_3[1:3,3:3]))
HCC78_DMSO_3_2Gy <- (as.numeric(HCC78_DMSO_3[4:6,3:3]))
HCC78_DMSO_3_4Gy <- (as.numeric(HCC78_DMSO_3[7:9,3:3]))
HCC78_DMSO_3_8Gy <- as.numeric(HCC78_DMSO_3[10:12,3:3])
PE_HCC78_DMSO_3<- mean(mean(HCC78_DMSO_3_0Gy)/as.numeric(HCC78_DMSO_3[1,2]))
mean_HCC78_DMSO_3 <- data.frame("dose_Gy"=c(0,2,4,8),
                                "N0"=c(cells_HCC78_DMSO_3),
                                "N_mean"= c(mean(HCC78_DMSO_3_0Gy,na.rm=TRUE),mean(HCC78_DMSO_3_2Gy,na.rm=TRUE),mean(HCC78_DMSO_3_4Gy,na.rm=TRUE),mean(HCC78_DMSO_3_8Gy,na.rm=TRUE)),
                                "N_SD"= c(sd(HCC78_DMSO_3_0Gy,na.rm=TRUE),sd(HCC78_DMSO_3_2Gy,na.rm=TRUE),sd(HCC78_DMSO_3_4Gy,na.rm=TRUE),sd(HCC78_DMSO_3_8Gy,na.rm=TRUE)))
mean_HCC78_DMSO_3["SEM"] <- mean_HCC78_DMSO_3[,4]/sqrt(3)
mean_HCC78_DMSO_3["N/N0"] <- mean_HCC78_DMSO_3[,3]/mean_HCC78_DMSO_3[,2]
mean_HCC78_DMSO_3["PE"] <- mean_HCC78_DMSO_3[1,3]/mean_HCC78_DMSO_3[1,2]
mean_HCC78_DMSO_3["SF"] <- mean_HCC78_DMSO_3["N_mean"]/(mean_HCC78_DMSO_3["PE"]*mean_HCC78_DMSO_3["N0"])

HCC78_DMSO_3["SF"]<- as.numeric(unlist(HCC78_DMSO_3["N"]))/(mean_HCC78_DMSO_3[1,7]*as.numeric(unlist(HCC78_DMSO_3["N0"])))
#####Caclulation HCC78 10nM AZD1390 III####
HCC78_10nM_AZD1390_3 <- HCC78[3:14,12:14]
colnames(HCC78_10nM_AZD1390_3)<- c("Dose","N0","N")
HCC78_10nM_AZD1390_3["N0/N1"]<-as.numeric(HCC78_10nM_AZD1390_3[,3])/as.numeric(HCC78_10nM_AZD1390_3[,2])
HCC78_10nM_AZD1390_3["Plating effeciency"]<- (HCC78_10nM_AZD1390_3[,4]/mean(as.numeric(HCC78_10nM_AZD1390_3[1:3,4])))
cells_HCC78_10nM_AZD1390_3 <- data.frame("N0"=as.numeric(HCC78_10nM_AZD1390_3[c(1,4,7,10),2]))

HCC78_10nM_AZD1390_3_0Gy <- (as.numeric(HCC78_10nM_AZD1390_3[1:3,3:3]))
HCC78_10nM_AZD1390_3_2Gy <- (as.numeric(HCC78_10nM_AZD1390_3[4:6,3:3]))
HCC78_10nM_AZD1390_3_4Gy <- (as.numeric(HCC78_10nM_AZD1390_3[7:9,3:3]))
HCC78_10nM_AZD1390_3_8Gy <- as.numeric(HCC78_10nM_AZD1390_3[10:12,3:3])
PE_HCC78_10nM_AZD1390_3<- mean(mean(HCC78_10nM_AZD1390_3_0Gy)/as.numeric(HCC78_10nM_AZD1390_3[1,2]))
mean_HCC78_10nM_AZD1390_3 <- data.frame("dose_Gy"=c(0,2,4,8),
                                        "N0"=c(cells_HCC78_10nM_AZD1390_3),
                                        "N_mean"= c(mean(HCC78_10nM_AZD1390_3_0Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_3_2Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_3_4Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_3_8Gy,na.rm=TRUE)),
                                        "N_SD"= c(sd(HCC78_10nM_AZD1390_3_0Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_3_2Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_3_4Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_3_8Gy,na.rm=TRUE)))
mean_HCC78_10nM_AZD1390_3["SEM"] <- mean_HCC78_10nM_AZD1390_3[,4]/sqrt(3)
mean_HCC78_10nM_AZD1390_3["N/N0"] <- mean_HCC78_10nM_AZD1390_3[,3]/mean_HCC78_10nM_AZD1390_3[,2]
mean_HCC78_10nM_AZD1390_3["PE"] <- mean_HCC78_10nM_AZD1390_3[1,3]/mean_HCC78_10nM_AZD1390_3[1,2]
mean_HCC78_10nM_AZD1390_3["SF"] <- mean_HCC78_10nM_AZD1390_3["N_mean"]/(mean_HCC78_10nM_AZD1390_3["PE"]*mean_HCC78_10nM_AZD1390_3["N0"])

HCC78_10nM_AZD1390_3["SF"]<- as.numeric(unlist(HCC78_10nM_AZD1390_3["N"]))/(mean_HCC78_10nM_AZD1390_3[1,7]*as.numeric(unlist(HCC78_10nM_AZD1390_3["N0"])))
#####Caclulation HCC78 DMSO IV####
HCC78_DMSO_4 <- HCC78[19:30,17:19]
colnames(HCC78_DMSO_4)<- c("Dose","N0","N")
HCC78_DMSO_4["N0/N1"]<-as.numeric(HCC78_DMSO_4[,3])/as.numeric(HCC78_DMSO_4[,2])
HCC78_DMSO_4["Plating effeciency"]<- (HCC78_DMSO_4[,4]/mean(as.numeric(HCC78_DMSO_4[1:3,4])))
cells_HCC78_DMSO_4 <- data.frame("N0"=as.numeric(HCC78_DMSO_4[c(1,4,7,10),2]))

HCC78_DMSO_4_0Gy <- (as.numeric(HCC78_DMSO_4[1:3,3:3]))
HCC78_DMSO_4_2Gy <- (as.numeric(HCC78_DMSO_4[4:6,3:3]))
HCC78_DMSO_4_4Gy <- (as.numeric(HCC78_DMSO_4[7:9,3:3]))
HCC78_DMSO_4_8Gy <- as.numeric(HCC78_DMSO_4[10:12,3:3])
PE_HCC78_DMSO_4<- mean(mean(HCC78_DMSO_4_0Gy)/as.numeric(HCC78_DMSO_4[1,2]))
mean_HCC78_DMSO_4 <- data.frame("dose_Gy"=c(0,2,4,8),
                                "N0"=c(cells_HCC78_DMSO_4),
                                "N_mean"= c(mean(HCC78_DMSO_4_0Gy,na.rm=TRUE),mean(HCC78_DMSO_4_2Gy,na.rm=TRUE),mean(HCC78_DMSO_4_4Gy,na.rm=TRUE),mean(HCC78_DMSO_4_8Gy,na.rm=TRUE)),
                                "N_SD"= c(sd(HCC78_DMSO_4_0Gy,na.rm=TRUE),sd(HCC78_DMSO_4_2Gy,na.rm=TRUE),sd(HCC78_DMSO_4_4Gy,na.rm=TRUE),sd(HCC78_DMSO_4_8Gy,na.rm=TRUE)))
mean_HCC78_DMSO_4["SEM"] <- mean_HCC78_DMSO_4[,4]/sqrt(3)
mean_HCC78_DMSO_4["N/N0"] <- mean_HCC78_DMSO_4[,3]/mean_HCC78_DMSO_4[,2]
mean_HCC78_DMSO_4["PE"] <- mean_HCC78_DMSO_4[1,3]/mean_HCC78_DMSO_4[1,2]
mean_HCC78_DMSO_4["SF"] <- mean_HCC78_DMSO_4["N_mean"]/(mean_HCC78_DMSO_4["PE"]*mean_HCC78_DMSO_4["N0"])

HCC78_DMSO_4["SF"]<- as.numeric(unlist(HCC78_DMSO_4["N"]))/(mean_HCC78_DMSO_4[1,7]*as.numeric(unlist(HCC78_DMSO_4["N0"])))
#####Caclulation HCC78 10nM AZD1390 IV####
HCC78_10nM_AZD1390_4 <- HCC78[3:14,17:19]
colnames(HCC78_10nM_AZD1390_4)<- c("Dose","N0","N")
HCC78_10nM_AZD1390_4["N0/N1"]<-as.numeric(HCC78_10nM_AZD1390_4[,3])/as.numeric(HCC78_10nM_AZD1390_4[,2])
HCC78_10nM_AZD1390_4["Plating effeciency"]<- (HCC78_10nM_AZD1390_4[,4]/mean(as.numeric(HCC78_10nM_AZD1390_4[1:3,4])))
cells_HCC78_10nM_AZD1390_4 <- data.frame("N0"=as.numeric(HCC78_10nM_AZD1390_4[c(1,4,7,10),2]))

HCC78_10nM_AZD1390_4_0Gy <- (as.numeric(HCC78_10nM_AZD1390_4[1:3,3:3]))
HCC78_10nM_AZD1390_4_2Gy <- (as.numeric(HCC78_10nM_AZD1390_4[4:6,3:3]))
HCC78_10nM_AZD1390_4_4Gy <- (as.numeric(HCC78_10nM_AZD1390_4[7:9,3:3]))
HCC78_10nM_AZD1390_4_8Gy <- as.numeric(HCC78_10nM_AZD1390_4[10:12,3:3])
PE_HCC78_10nM_AZD1390_4<- mean(mean(HCC78_10nM_AZD1390_4_0Gy)/as.numeric(HCC78_10nM_AZD1390_4[1,2]))
mean_HCC78_10nM_AZD1390_4 <- data.frame("dose_Gy"=c(0,2,4,8),
                                        "N0"=c(cells_HCC78_10nM_AZD1390_4),
                                        "N_mean"= c(mean(HCC78_10nM_AZD1390_4_0Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_4_2Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_4_4Gy,na.rm=TRUE),mean(HCC78_10nM_AZD1390_4_8Gy,na.rm=TRUE)),
                                        "N_SD"= c(sd(HCC78_10nM_AZD1390_4_0Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_4_2Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_4_4Gy,na.rm=TRUE),sd(HCC78_10nM_AZD1390_4_8Gy,na.rm=TRUE)))
mean_HCC78_10nM_AZD1390_4["SEM"] <- mean_HCC78_10nM_AZD1390_4[,4]/sqrt(3)
mean_HCC78_10nM_AZD1390_4["N/N0"] <- mean_HCC78_10nM_AZD1390_4[,3]/mean_HCC78_10nM_AZD1390_4[,2]
mean_HCC78_10nM_AZD1390_4["PE"] <- mean_HCC78_10nM_AZD1390_4[1,3]/mean_HCC78_10nM_AZD1390_4[1,2]
mean_HCC78_10nM_AZD1390_4["SF"] <- mean_HCC78_10nM_AZD1390_4["N_mean"]/(mean_HCC78_10nM_AZD1390_4["PE"]*mean_HCC78_10nM_AZD1390_4["N0"])

HCC78_10nM_AZD1390_4["SF"]<- as.numeric(unlist(HCC78_10nM_AZD1390_4["N"]))/(mean_HCC78_10nM_AZD1390_4[1,7]*as.numeric(unlist(HCC78_10nM_AZD1390_4["N0"])))
#####binding HCC78 results####
all_HCC78_DMSO <- cbind(mean_HCC78_DMSO_1[,c(1:3,8)],mean_HCC78_DMSO_2[,c(3,8)],mean_HCC78_DMSO_3[,c(3,8)],mean_HCC78_DMSO_4[,c(3,8)])
colnames(all_HCC78_DMSO)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_HCC78_DMSO["SF_mean"]<- rowMeans(all_HCC78_DMSO[,c(4,6,8,10)])
all_HCC78_DMSO["SF_sd"]<- apply(all_HCC78_DMSO[,c(4,6,8,10)], 1, sd)
all_HCC78_DMSO["SF_SEM"]<- all_HCC78_DMSO["SF_sd"]/sqrt(4)

all_HCC78_10nM_AZD1390 <- cbind(mean_HCC78_10nM_AZD1390_1[,c(1:3,8)],mean_HCC78_10nM_AZD1390_2[,c(3,8)],mean_HCC78_10nM_AZD1390_3[,c(3,8)],mean_HCC78_10nM_AZD1390_4[,c(3,8)])
colnames(all_HCC78_10nM_AZD1390)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_HCC78_10nM_AZD1390["SF_mean"]<- rowMeans(all_HCC78_10nM_AZD1390[,c(4,6,8,10)])
all_HCC78_10nM_AZD1390["SF_sd"]<- apply(all_HCC78_10nM_AZD1390[,c(4,6,8,10)], 1, sd)
all_HCC78_10nM_AZD1390["SF_SEM"]<- all_HCC78_10nM_AZD1390["SF_sd"]/sqrt(4)
#####Plotting####
setwd(save_wd)
jpeg(filename = "CFA_AZD1390_10nM.jpeg",  width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2,2,3,3),ncol=6,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt <- c("DMSO","10nM AZD1390")
plot(x=dose[1:2], y=all_A549_10nM_AZD1390[1:2,11],xlim=c(0,8),mgp=c(7, 0,1), lwd=3, type="l", cex.lab=4,ylim=c(0.0001,100),axes=FALSE, cex=6, 
     pch=23,log="y", col="black",  lty=3, ann=FALSE,cex.main=4, main="A594")
lines(x=dose, y=all_A549_DMSO[,11], lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black", lty=3)
lines(x=c(dose[2:3]), y=c(all_A549_10nM_AZD1390[2,11],0.0001), lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black",  lty=3)
points(x=dose, y=all_A549_DMSO[,11], lwd=3,ylim=c(0.001,10), cex=6, pch=21, col="black", bg="#66C2A5", lty=3)
points(x=dose, y=all_A549_10nM_AZD1390[,11], lwd=3,ylim=c(0.001,10), cex=6, pch=23, col="black", bg="#9E0142", lty=3)
axis(1,lwd=3,at = c(0,4,8),mgp=c(7, 4,0),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5)
axis(2,lwd=3,at = c(100,10,1, 0.1, 0.01,0.001,0.0001),label= c("","","1", "", "",parse(text="10^-3"),""),mgp=c(7, 2,0),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,las=2)
arrows(x0=dose[2:4], y0=all_A549_10nM_AZD1390[2:4,11]-all_A549_10nM_AZD1390[2:4,13], x1=dose[2:4], y1=all_A549_10nM_AZD1390[2:4,11]+all_A549_10nM_AZD1390[2:4,13], code=3, angle=90, length=0.1,lwd=3)
arrows(x0=dose[2:4], y0=all_A549_DMSO[2:4,11]-all_A549_DMSO[2:4,13], x1=dose[2:4], y1=all_A549_DMSO[2:4,11]+all_A549_DMSO[2:4,13], code=3, angle=90, length=0.1,lwd=3)
box(lwd=3,lty="solid")
title(main="A549", cex.lab=3,cex.main=3,line=2)
title( ylab="Surviving Fraction [%]", cex.lab=3.5,cex.main=6,line=9,adj=0.75)
title( xlab="Dose [Gy]", cex.lab=3.5,cex.main=6,line=7)
minor.tick(nx =2/1,ny=0,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl=-0.75)) 
legend(x=c(0.5,8),y=c(3,100),leg.txt,pch=c(21,23),adj=c(0,0.25),box.col="white",pt.bg=c("#66C2A5","#9E0142"),inset=c(-0.0,0.0),pt.cex=4, cex =2.5,box.lwd=2,pt.lwd=2)

leg.txt <- c("DMSO","10nM AZD1390")
plot(x=dose[1:2], y=all_HCC78_10nM_AZD1390[1:2,11],xlim=c(0,8),mgp=c(7, 0,1), type="l", lwd=3,cex.lab=4,ylim=c(0.0001,100),axes=FALSE, cex=6, 
     pch=23,log="y", col="black",  lty=3, ann=FALSE,cex.main=4, main="HCC78")
lines(x=dose, y=all_HCC78_DMSO[,11], lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black",  lty=3)
points(x=dose[1:3], y=all_HCC78_DMSO[1:3,11], lwd=3,ylim=c(0.001,10), cex=6, pch=21, col="black", bg="#66C2A5", lty=3)
points(x=dose[1:3], y=all_HCC78_10nM_AZD1390[1:3,11], lwd=3,ylim=c(0.001,10), cex=6, pch=23, col="black", bg="#9E0142", lty=3)
axis(1,lwd=3,at = c(0,4,8),mgp=c(7, 4,0),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5)
axis(2,lwd=3,at = c(100,10,1, 0.1, 0.01,0.001,0.0001),label= c("","","1", "", "",parse(text="10^-3"),""),mgp=c(7, 2,0),lwd.ticks = 3,cex=6, cex.axis=3,tcl = -1.5,las=2)
arrows(x0=dose[2:4], y0=all_HCC78_10nM_AZD1390[2:4,11]-all_HCC78_10nM_AZD1390[2:4,13], x1=dose[2:4], y1=all_HCC78_10nM_AZD1390[2:4,11]+all_HCC78_10nM_AZD1390[2:4,13], code=3, angle=90, length=0.1,lwd=3)
arrows(x0=dose[2:4], y0=all_HCC78_DMSO[2:4,11]-all_HCC78_DMSO[2:4,13], x1=dose[2:4], y1=all_HCC78_DMSO[2:4,11]+all_HCC78_DMSO[2:4,13], code=3, angle=90, length=0.1,lwd=3)
lines(x = c(2,3), y = c(all_HCC78_10nM_AZD1390[2,11] ,0.0001), lwd=3, cex=1.4, pch=21, col="black", lty=3)
lines(x = c(4,5), y = c(all_HCC78_DMSO[3,11] ,0.0001), lwd=3, cex=1.4, pch=21, col="black", lty=3)
box(lwd=3,lty="solid")
title(main="HCC78", cex.lab=3,cex.main=3,line=2)
title( ylab="Surviving Fraction [%]", cex.lab=3.5,cex.main=6,line=9,adj=0.75)
title( xlab="Dose [Gy]", cex.lab=3.5,cex.main=6,line=7)
minor.tick(nx =2/1,ny=0,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 

leg.txt <- c("DMSO","10nM AZD1390")
plot(x=dose[1:2], y=all_HCC827_10nM_AZD1390[1:2,4],xlim=c(0,8),mgp=c(7, 0,1), type="l", cex.lab=2,lwd=3,ylim=c(0.0001,100),axes=FALSE, cex=6,
     pch=23,log="y", col="black",  lty=3, ann=FALSE)
lines(x=dose, y=all_HCC827_DMSO[,4], lwd=3,ylim=c(0.0001,10), cex=2, pch=21, col="black", lty=3)
points(x=dose[1:3], y=all_HCC827_DMSO[1:3,4], lwd=3,ylim=c(0.001,10), cex=6, pch=21, col="black", bg="#66C2A5", lty=2)
points(x=dose[1:3], y=all_HCC827_10nM_AZD1390[1:3,4], lwd=3,ylim=c(0.001,10), cex=6, pch=23, col="black", bg="#9E0142", lty=2)
axis(1,lwd=3,at = c(0,4,8),mgp=c(7, 4,0),lwd.ticks = 3,cex=6, cex.axis=3.5,tcl = -1.5)
axis(2,lwd=3,at = c(100,10,1, 0.1, 0.01,0.001,0.0001),label= c("","","1", "", "",parse(text="10^-3"),""),mgp=c(7, 2,0),lwd.ticks = 3,cex=6, cex.axis=3.5,tcl = -1.5,las=2)
arrows(x0=dose[2:4], y0=all_HCC827_10nM_AZD1390[2:4,4]-all_HCC827_10nM_AZD1390[2:4,5], x1=dose[2:4], y1=all_HCC827_10nM_AZD1390[2:4,4]+all_HCC827_10nM_AZD1390[2:4,5], code=3, angle=90, length=0.1,lwd=3)
arrows(x0=dose[2:4], y0=all_HCC827_DMSO[2:4,4]-all_HCC827_DMSO[2:4,5], x1=dose[2:4], y1=all_HCC827_DMSO[2:4,4]+all_HCC827_DMSO[2:4,5], code=3, angle=90, length=0.1,lwd=3)
lines(x = c(2,3), y = c(all_HCC827_10nM_AZD1390[2,4] ,0.0001), lwd=3, cex=1.4, pch=21, col="black", lty=3)
lines(x = c(4,5), y = c(all_HCC827_DMSO[3,4] ,0.0001), lwd=3, cex=1.4, pch=21, col="black", lty=3)
box(lwd=3,lty="solid")
title(main="HCC827", cex.lab=3,cex.main=3,line=2)
title( ylab="Surviving Fraction [%]", cex.lab=3.5,cex.main=6,line=9,adj=0.75)
title( xlab="Dose [Gy]", cex.lab=3.5,cex.main=6,line=7)
minor.tick(nx =2/1,ny=0,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
dev.off()
setwd(source_wd)