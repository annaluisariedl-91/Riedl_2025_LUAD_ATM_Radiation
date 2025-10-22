####Set directory#####
setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/CFA/Exp1_DoseResponse8Gy/raw")
dose <-c(0,2,4,6,8)
####I H23 - Extracting data#####
H23<- read.csv("CFA_1_8Gy_H23.csv",sep=",")
H23_1 <- H23[2:31,1:3]
colnames(H23_1)<- c("Dose","N0","N")
H23_1["N0/N1"]<-as.numeric(H23_1[,3])/as.numeric(H23_1[,2])
H23_1["Plating effeciency"]<- (H23_1[,4]/mean(as.numeric(H23_1[1:6,4])))
cells_H23_1 <- data.frame(H23_1[c(1,7,13,19,25),2])

H23_1_0Gy <- (as.numeric(H23_1[1:6,3]))
H23_1_2Gy <- (as.numeric(H23_1[7:12,3]))
H23_1_4Gy <- (as.numeric(H23_1[13:18,3]))
H23_1_6Gy <- (as.numeric(H23_1[19:24,3]))
H23_1_8Gy <- (as.numeric(H23_1[25:30,3]))
mean_H23_1 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=as.numeric(unlist(cells_H23_1)),
                           "N_mean"= c(mean(H23_1_0Gy,na.rm=TRUE),mean(H23_1_2Gy,na.rm=TRUE),mean(H23_1_4Gy,na.rm=TRUE),mean(H23_1_6Gy,na.rm=TRUE),mean(H23_1_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(H23_1_0Gy,na.rm=TRUE),sd(H23_1_2Gy,na.rm=TRUE),sd(H23_1_4Gy,na.rm=TRUE),sd(H23_1_6Gy,na.rm=TRUE),sd(H23_1_8Gy,na.rm=TRUE)))
mean_H23_1["SEM"] <- mean_H23_1[,4]/sqrt(6)
mean_H23_1["N/N0"] <- mean_H23_1[,3]/(mean_H23_1[,2])
mean_H23_1["PE"] <- mean_H23_1[1,3]/(mean_H23_1[1,2])
mean_H23_1["SF"] <- mean_H23_1["N_mean"]/(mean_H23_1["PE"]*mean_H23_1["N0"])
CFA_mean_H23_1<-mean_H23_1
H23_1["SF"]<- as.numeric(unlist(H23_1["N"]))/(mean_H23_1[1,7]*as.numeric(unlist(H23_1["N0"])))
####I HCC44 - Extracting Data####
HCC44<- read.csv("CFA_1_8Gy_HCC44.csv")
HCC44_1 <- HCC44[2:31,1:3]
colnames(HCC44_1)<- c("Dose","N0","N")
HCC44_1["N0/N1"]<-as.numeric(HCC44_1[,3])/as.numeric(HCC44_1[,2])
HCC44_1["Plating effeciency"]<- (HCC44_1[,4]/mean(as.numeric(HCC44_1[1:6,4])))
cells_HCC44_1 <- data.frame("N0"=as.numeric(HCC44_1[c(1,7,13,19,25),2]))

HCC44_1_0Gy <- (as.numeric(HCC44_1[1:6,3]))
HCC44_1_2Gy <- (as.numeric(HCC44_1[7:12,3]))
HCC44_1_4Gy <- (as.numeric(HCC44_1[13:18,3]))
HCC44_1_6Gy <- (as.numeric(HCC44_1[19:24,3]))
HCC44_1_8Gy <- as.numeric(HCC44_1[25:30,3])
PE_HCC44_1<- mean(mean(HCC44_1_0Gy)/as.numeric(HCC44_1[1,2]))
mean_HCC44_1 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_HCC44_1),
                           "N_mean"= c(mean(HCC44_1_0Gy,na.rm=TRUE),mean(HCC44_1_2Gy,na.rm=TRUE),mean(HCC44_1_4Gy,na.rm=TRUE),mean(HCC44_1_6Gy,na.rm=TRUE),mean(HCC44_1_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC44_1_0Gy,na.rm=TRUE),sd(HCC44_1_2Gy,na.rm=TRUE),sd(HCC44_1_4Gy,na.rm=TRUE),sd(HCC44_1_6Gy,na.rm=TRUE),sd(HCC44_1_8Gy,na.rm=TRUE)))
mean_HCC44_1["SEM"] <- mean_HCC44_1[,4]/sqrt(6)
mean_HCC44_1["N/N0"] <- mean_HCC44_1[,3]/mean_HCC44_1[,2]
mean_HCC44_1["PE"] <- mean_HCC44_1[1,3]/mean_HCC44_1[1,2]
mean_HCC44_1["SF"] <- mean_HCC44_1["N_mean"]/(mean_HCC44_1["PE"]*mean_HCC44_1["N0"])

HCC44_1["SF"]<- as.numeric(unlist(HCC44_1["N"]))/(mean_HCC44_1[1,7]*as.numeric(unlist(HCC44_1["N0"])))
####II HCC44 - Extracting Data####
HCC44_2 <- HCC44[2:31,7:9]
colnames(HCC44_2)<- c("Dose","N0","N")
HCC44_2["N0/N1"]<-as.numeric(HCC44_2[,3])/as.numeric(HCC44_2[,2])
HCC44_2["Plating effeciency"]<- (HCC44_2[,4]/mean(as.numeric(HCC44_2[1:6,4])))
cells_HCC44_2 <- data.frame("N0"=as.numeric(HCC44_2[c(1,7,13,19,25),2]))

HCC44_2_0Gy <- (as.numeric(HCC44_2[1:6,3:3]))
HCC44_2_2Gy <- (as.numeric(HCC44_2[7:12,3:3]))
HCC44_2_4Gy <- (as.numeric(HCC44_2[13:18,3:3]))
HCC44_2_6Gy <- (as.numeric(HCC44_2[19:24,3:3]))
HCC44_2_8Gy <- as.numeric(HCC44_2[25:30,3:3])
PE_HCC44_2<- mean(mean(HCC44_2_0Gy)/as.numeric(HCC44_2[1,2]))
mean_HCC44_2 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_HCC44_2),
                           "N_mean"= c(mean(HCC44_2_0Gy,na.rm=TRUE),mean(HCC44_2_2Gy,na.rm=TRUE),mean(HCC44_2_4Gy,na.rm=TRUE),mean(HCC44_2_6Gy,na.rm=TRUE),mean(HCC44_2_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC44_2_0Gy,na.rm=TRUE),sd(HCC44_2_2Gy,na.rm=TRUE),sd(HCC44_2_4Gy,na.rm=TRUE),sd(HCC44_2_6Gy,na.rm=TRUE),sd(HCC44_2_8Gy,na.rm=TRUE)))
mean_HCC44_2["SEM"] <- mean_HCC44_2[,4]/sqrt(6)
mean_HCC44_2["N/N0"] <- mean_HCC44_2[,3]/mean_HCC44_2[,2]
mean_HCC44_2["PE"] <- mean_HCC44_2[1,3]/mean_HCC44_2[1,2]
mean_HCC44_2["SF"] <- mean_HCC44_2["N_mean"]/(mean_HCC44_2["PE"]*mean_HCC44_2["N0"])

HCC44_2["SF"]<- as.numeric(unlist(HCC44_2["N"]))/(mean_HCC44_2[1,7]*as.numeric(unlist(HCC44_2["N0"])))
####III HCC44 - Extracting Data####
HCC44_3 <- HCC44[2:31,13:15]
colnames(HCC44_3)<- c("Dose","N0","N")
HCC44_3["N0/N1"]<-as.numeric(HCC44_3[,3])/as.numeric(HCC44_3[,2])
HCC44_3["Plating effeciency"]<- (HCC44_3[,4]/mean(as.numeric(HCC44_3[1:6,4])))
cells_HCC44_3 <- data.frame("N0"=as.numeric(HCC44_3[c(1,7,13,19,25),2]))

HCC44_3_0Gy <- (as.numeric(HCC44_3[1:6,3:3]))
HCC44_3_2Gy <- (as.numeric(HCC44_3[7:12,3:3]))
HCC44_3_4Gy <- (as.numeric(HCC44_3[13:18,3:3]))
HCC44_3_6Gy <- (as.numeric(HCC44_3[19:24,3:3]))
HCC44_3_8Gy <- as.numeric(HCC44_3[25:30,3:3])
PE_HCC44_3<- mean(mean(HCC44_3_0Gy)/as.numeric(HCC44_3[1,2]))
mean_HCC44_3 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_HCC44_3),
                           "N_mean"= c(mean(HCC44_3_0Gy,na.rm=TRUE),mean(HCC44_3_2Gy,na.rm=TRUE),mean(HCC44_3_4Gy,na.rm=TRUE),mean(HCC44_3_6Gy,na.rm=TRUE),mean(HCC44_3_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC44_3_0Gy,na.rm=TRUE),sd(HCC44_3_2Gy,na.rm=TRUE),sd(HCC44_3_4Gy,na.rm=TRUE),sd(HCC44_3_6Gy,na.rm=TRUE),sd(HCC44_3_8Gy,na.rm=TRUE)))
mean_HCC44_3["SEM"] <- mean_HCC44_3[,4]/sqrt(6)
mean_HCC44_3["N/N0"] <- mean_HCC44_3[,3]/mean_HCC44_3[,2]
mean_HCC44_3["PE"] <- mean_HCC44_3[1,3]/mean_HCC44_3[1,2]
mean_HCC44_3["SF"] <- mean_HCC44_3["N_mean"]/(mean_HCC44_3["PE"]*mean_HCC44_3["N0"])

HCC44_3["SF"]<- as.numeric(unlist(HCC44_3["N"]))/(mean_HCC44_3[1,7]*as.numeric(unlist(HCC44_3["N0"])))
####IV HCC44 - Extracting Data####
HCC44_4 <- HCC44[2:31,13:15]
colnames(HCC44_4)<- c("Dose","N0","N")
HCC44_4["N0/N1"]<-as.numeric(HCC44_4[,3])/as.numeric(HCC44_4[,2])
HCC44_4["Plating effeciency"]<- (HCC44_4[,4]/mean(as.numeric(HCC44_4[1:6,4])))
cells_HCC44_4 <- data.frame("N0"=as.numeric(HCC44_4[c(1,7,13,19,25),2]))

HCC44_4_0Gy <- (as.numeric(HCC44_4[1:6,3]))
HCC44_4_2Gy <- (as.numeric(HCC44_4[7:12,3]))
HCC44_4_4Gy <- (as.numeric(HCC44_4[13:18,3]))
HCC44_4_6Gy <- (as.numeric(HCC44_4[19:24,3]))
HCC44_4_8Gy <- as.numeric(HCC44_4[25:30,3])
PE_HCC44_4<- mean(mean(HCC44_4_0Gy)/as.numeric(HCC44_4[1,2]))
mean_HCC44_4 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_HCC44_4),
                           "N_mean"= c(mean(HCC44_4_0Gy,na.rm=TRUE),mean(HCC44_4_2Gy,na.rm=TRUE),mean(HCC44_4_4Gy,na.rm=TRUE),mean(HCC44_4_6Gy,na.rm=TRUE),mean(HCC44_4_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC44_4_0Gy,na.rm=TRUE),sd(HCC44_4_2Gy,na.rm=TRUE),sd(HCC44_4_4Gy,na.rm=TRUE),sd(HCC44_4_6Gy,na.rm=TRUE),sd(HCC44_4_8Gy,na.rm=TRUE)))
mean_HCC44_4["SEM"] <- mean_HCC44_4[,4]/sqrt(6)
mean_HCC44_4["N/N0"] <- mean_HCC44_4[,3]/mean_HCC44_4[,2]
mean_HCC44_4["PE"] <- mean_HCC44_4[1,3]/mean_HCC44_4[1,2]
mean_HCC44_4["SF"] <- mean_HCC44_4["N_mean"]/(mean_HCC44_4["PE"]*mean_HCC44_4["N0"])

HCC44_4["SF"]<- as.numeric(unlist(HCC44_4["N"]))/(mean_HCC44_4[1,7]*as.numeric(unlist(HCC44_4["N0"])))
#####binding HCC44 results####
all_HCC44 <- cbind(mean_HCC44_1[,c(1:3,8)],mean_HCC44_2[,c(3,8)],mean_HCC44_3[,c(3,8)],mean_HCC44_4[,c(3,8)])
colnames(all_HCC44)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_HCC44["SF_mean"]<- rowMeans(all_HCC44[,c(4,6,8,10)],na.rm=TRUE)
all_HCC44["SF_sd"]<- apply(all_HCC44[,c(4,6,8,10)], 1, sd,na.rm=TRUE)
all_HCC44["SF_SEM"]<- all_HCC44["SF_sd"]/sqrt(4)
CFA_all_HCC44<-all_HCC44
####I A594 - Extracting Data####
A549<- read.csv("CFA_1_8Gy_A549.csv")
A549_1 <- A549[2:21,1:3]
colnames(A549_1)<- c("Dose","N0","N")
A549_1["N0/N1"]<-as.numeric(A549_1[,3])/as.numeric(A549_1[,2])
A549_1["Plating effeciency"]<- (A549_1[,4]/mean(as.numeric(A549_1[1:6,4])))
cells_A549_1 <- data.frame("N0"=as.numeric(A549_1[c(1,5,9,13,17),2]))

A549_1_0Gy <- (as.numeric(A549_1[1:4,3]))
A549_1_2Gy <- (as.numeric(A549_1[5:8,3]))
A549_1_4Gy <- (as.numeric(A549_1[9:12,3]))
A549_1_6Gy <- (as.numeric(A549_1[13:16,3]))
A549_1_8Gy <- as.numeric(A549_1[17:20,3])
PE_A549_1<- mean(mean(A549_1_0Gy)/as.numeric(A549_1[1,2]))
mean_A549_1 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_A549_1),
                           "N_mean"= c(mean(A549_1_0Gy,na.rm=TRUE),mean(A549_1_2Gy,na.rm=TRUE),mean(A549_1_4Gy,na.rm=TRUE),mean(A549_1_6Gy,na.rm=TRUE),mean(A549_1_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(A549_1_0Gy,na.rm=TRUE),sd(A549_1_2Gy,na.rm=TRUE),sd(A549_1_4Gy,na.rm=TRUE),sd(A549_1_6Gy,na.rm=TRUE),sd(A549_1_8Gy,na.rm=TRUE)))
mean_A549_1["SEM"] <- mean_A549_1[,4]/sqrt(6)
mean_A549_1["N/N0"] <- mean_A549_1[,3]/mean_A549_1[,2]
mean_A549_1["PE"] <- mean_A549_1[1,3]/mean_A549_1[1,2]
mean_A549_1["SF"] <- mean_A549_1["N_mean"]/(mean_A549_1["PE"]*mean_A549_1["N0"])

A549_1["SF"]<- as.numeric(unlist(A549_1["N"]))/(mean_A549_1[1,7]*as.numeric(unlist(A549_1["N0"])))
####II A549 - Extracting Data####
A549_2 <- A549[2:21,7:9]
colnames(A549_2)<- c("Dose","N0","N")
A549_2["N0/N1"]<-as.numeric(A549_2[,3])/as.numeric(A549_2[,2])
A549_2["Plating effeciency"]<- (A549_2[,4]/mean(as.numeric(A549_2[1:6,4])))
cells_A549_2 <- data.frame("N0"=as.numeric(A549_2[c(1,5,9,13,17),2]))

A549_2_0Gy <- (as.numeric(A549_2[1:4,3:3]))
A549_2_2Gy <- (as.numeric(A549_2[5:8,3:3]))
A549_2_4Gy <- (as.numeric(A549_2[9:12,3:3]))
A549_2_6Gy <- (as.numeric(A549_2[13:16,3:3]))
A549_2_8Gy <- as.numeric(A549_2[17:20,3:3])
PE_A549_2<- mean(mean(A549_2_0Gy)/as.numeric(A549_2[1,2]))
mean_A549_2 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_A549_2),
                           "N_mean"= c(mean(A549_2_0Gy,na.rm=TRUE),mean(A549_2_2Gy,na.rm=TRUE),mean(A549_2_4Gy,na.rm=TRUE),mean(A549_2_6Gy,na.rm=TRUE),mean(A549_2_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(A549_2_0Gy,na.rm=TRUE),sd(A549_2_2Gy,na.rm=TRUE),sd(A549_2_4Gy,na.rm=TRUE),sd(A549_2_6Gy,na.rm=TRUE),sd(A549_2_8Gy,na.rm=TRUE)))
mean_A549_2["SEM"] <- mean_A549_2[,4]/sqrt(6)
mean_A549_2["N/N0"] <- mean_A549_2[,3]/mean_A549_2[,2]
mean_A549_2["PE"] <- mean_A549_2[1,3]/mean_A549_2[1,2]
mean_A549_2["SF"] <- mean_A549_2["N_mean"]/(mean_A549_2["PE"]*mean_A549_2["N0"])

A549_2["SF"]<- as.numeric(unlist(A549_2["N"]))/(mean_A549_2[1,7]*as.numeric(unlist(A549_2["N0"])))
####III A549 - Extracting Data####
A549_3 <- A549[2:21,13:15]
colnames(A549_3)<- c("Dose","N0","N")
A549_3["N0/N1"]<-as.numeric(A549_3[,3])/as.numeric(A549_3[,2])
A549_3["Plating effeciency"]<- (A549_3[,4]/mean(as.numeric(A549_3[1:6,4])))
cells_A549_3 <- data.frame("N0"=as.numeric(A549_3[c(1,5,9,13,17),2]))

A549_3_0Gy <- (as.numeric(A549_3[1:4,3:3]))
A549_3_2Gy <- (as.numeric(A549_3[5:8,3:3]))
A549_3_4Gy <- (as.numeric(A549_3[9:12,3:3]))
A549_3_6Gy <- (as.numeric(A549_3[13:16,3:3]))
A549_3_8Gy <- as.numeric(A549_3[17:20,3:3])
PE_A549_3<- mean(mean(A549_3_0Gy)/as.numeric(A549_3[1,2]))
mean_A549_3 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_A549_3),
                           "N_mean"= c(mean(A549_3_0Gy,na.rm=TRUE),mean(A549_3_2Gy,na.rm=TRUE),mean(A549_3_4Gy,na.rm=TRUE),mean(A549_3_6Gy,na.rm=TRUE),mean(A549_3_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(A549_3_0Gy,na.rm=TRUE),sd(A549_3_2Gy,na.rm=TRUE),sd(A549_3_4Gy,na.rm=TRUE),sd(A549_3_6Gy,na.rm=TRUE),sd(A549_3_8Gy,na.rm=TRUE)))
mean_A549_3["SEM"] <- mean_A549_3[,4]/sqrt(6)
mean_A549_3["N/N0"] <- mean_A549_3[,3]/mean_A549_3[,2]
mean_A549_3["PE"] <- mean_A549_3[1,3]/mean_A549_3[1,2]
mean_A549_3["SF"] <- mean_A549_3["N_mean"]/(mean_A549_3["PE"]*mean_A549_3["N0"])

A549_3["SF"]<- as.numeric(unlist(A549_3["N"]))/(mean_A549_3[1,7]*as.numeric(unlist(A549_3["N0"])))
####IV A549 - Extracting Data####
A549_4 <- A549[2:21,13:15]
colnames(A549_4)<- c("Dose","N0","N")
A549_4["N0/N1"]<-as.numeric(A549_4[,3])/as.numeric(A549_4[,2])
A549_4["Plating effeciency"]<- (A549_4[,4]/mean(as.numeric(A549_4[1:6,4])))
cells_A549_4 <- data.frame("N0"=as.numeric(A549_4[c(1,5,9,13,17),2]))

A549_4_0Gy <- (as.numeric(A549_4[1:4,3]))
A549_4_2Gy <- (as.numeric(A549_4[5:8,3]))
A549_4_4Gy <- (as.numeric(A549_4[9:12,3]))
A549_4_6Gy <- (as.numeric(A549_4[13:16,3]))
A549_4_8Gy <- as.numeric(A549_4[17:20,3])
PE_A549_4<- mean(mean(A549_4_0Gy)/as.numeric(A549_4[1,2]))
mean_A549_4 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_A549_4),
                           "N_mean"= c(mean(A549_4_0Gy,na.rm=TRUE),mean(A549_4_2Gy,na.rm=TRUE),mean(A549_4_4Gy,na.rm=TRUE),mean(A549_4_6Gy,na.rm=TRUE),mean(A549_4_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(A549_4_0Gy,na.rm=TRUE),sd(A549_4_2Gy,na.rm=TRUE),sd(A549_4_4Gy,na.rm=TRUE),sd(A549_4_6Gy,na.rm=TRUE),sd(A549_4_8Gy,na.rm=TRUE)))
mean_A549_4["SEM"] <- mean_A549_4[,4]/sqrt(6)
mean_A549_4["N/N0"] <- mean_A549_4[,3]/mean_A549_4[,2]
mean_A549_4["PE"] <- mean_A549_4[1,3]/mean_A549_4[1,2]
mean_A549_4["SF"] <- mean_A549_4["N_mean"]/(mean_A549_4["PE"]*mean_A549_4["N0"])

A549_4["SF"]<- as.numeric(unlist(A549_4["N"]))/(mean_A549_4[1,7]*as.numeric(unlist(A549_4["N0"])))
#####binding A549 results####
all_A549 <- cbind(mean_A549_1[,c(1:3,8)],mean_A549_2[,c(3,8)],mean_A549_3[,c(3,8)],mean_A549_4[,c(3,8)])
colnames(all_A549)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_A549["SF_mean"]<- rowMeans(all_A549[,c(4,6,8,10)],na.rm=TRUE)
all_A549["SF_sd"]<- apply(all_A549[,c(4,6,8,10)], 1, sd,na.rm=TRUE)
all_A549["SF_SEM"]<- all_A549["SF_sd"]/sqrt(4)
CFA_all_A549<-all_A549
####I HCC827 - Extracting Data####
HCC827<- read.csv("CFA_1_8Gy_HCC827.csv")
HCC827_1 <- HCC827[2:31,1:3]
colnames(HCC827_1)<- c("Dose","N0","N")
HCC827_1["N0/N1"]<-as.numeric(HCC827_1[,3])/as.numeric(HCC827_1[,2])
HCC827_1["Plating effeciency"]<- (HCC827_1[,4]/mean(as.numeric(HCC827_1[1:6,4])))
cells_HCC827_1 <- data.frame("N0"=as.numeric(HCC827_1[c(1,7,13,19,25),2]))

HCC827_1_0Gy <- (as.numeric(HCC827_1[1:6,3]))
HCC827_1_2Gy <- (as.numeric(HCC827_1[7:12,3]))
HCC827_1_4Gy <- (as.numeric(HCC827_1[13:18,3]))
HCC827_1_6Gy <- (as.numeric(HCC827_1[19:24,3]))
HCC827_1_8Gy <- as.numeric(HCC827_1[25:30,3])
PE_HCC827_1<- mean(mean(HCC827_1_0Gy)/as.numeric(HCC827_1[1,2]))
mean_HCC827_1 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_HCC827_1),
                           "N_mean"= c(mean(HCC827_1_0Gy,na.rm=TRUE),mean(HCC827_1_2Gy,na.rm=TRUE),mean(HCC827_1_4Gy,na.rm=TRUE),mean(HCC827_1_6Gy,na.rm=TRUE),mean(HCC827_1_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC827_1_0Gy,na.rm=TRUE),sd(HCC827_1_2Gy,na.rm=TRUE),sd(HCC827_1_4Gy,na.rm=TRUE),sd(HCC827_1_6Gy,na.rm=TRUE),sd(HCC827_1_8Gy,na.rm=TRUE)))
mean_HCC827_1["SEM"] <- mean_HCC827_1[,4]/sqrt(6)
mean_HCC827_1["N/N0"] <- mean_HCC827_1[,3]/mean_HCC827_1[,2]
mean_HCC827_1["PE"] <- mean_HCC827_1[1,3]/mean_HCC827_1[1,2]
mean_HCC827_1["SF"] <- mean_HCC827_1["N_mean"]/(mean_HCC827_1["PE"]*mean_HCC827_1["N0"])

HCC827_1["SF"]<- as.numeric(unlist(HCC827_1["N"]))/(mean_HCC827_1[1,7]*as.numeric(unlist(HCC827_1["N0"])))
####II HCC827 - Extracting Data####
HCC827_2 <- HCC827[2:31,7:9]
colnames(HCC827_2)<- c("Dose","N0","N")
HCC827_2["N0/N1"]<-as.numeric(HCC827_2[,3])/as.numeric(HCC827_2[,2])
HCC827_2["Plating effeciency"]<- (HCC827_2[,4]/mean(as.numeric(HCC827_2[1:6,4])))
cells_HCC827_2 <- data.frame("N0"=as.numeric(HCC827_2[c(1,7,13,19,25),2]))

HCC827_2_0Gy <- (as.numeric(HCC827_2[1:6,3:3]))
HCC827_2_2Gy <- (as.numeric(HCC827_2[7:12,3:3]))
HCC827_2_4Gy <- (as.numeric(HCC827_2[13:18,3:3]))
HCC827_2_6Gy <- (as.numeric(HCC827_2[19:24,3:3]))
HCC827_2_8Gy <- as.numeric(HCC827_2[25:30,3:3])
PE_HCC827_2<- mean(mean(HCC827_2_0Gy)/as.numeric(HCC827_2[1,2]))
mean_HCC827_2 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_HCC827_2),
                           "N_mean"= c(mean(HCC827_2_0Gy,na.rm=TRUE),mean(HCC827_2_2Gy,na.rm=TRUE),mean(HCC827_2_4Gy,na.rm=TRUE),mean(HCC827_2_6Gy,na.rm=TRUE),mean(HCC827_2_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC827_2_0Gy,na.rm=TRUE),sd(HCC827_2_2Gy,na.rm=TRUE),sd(HCC827_2_4Gy,na.rm=TRUE),sd(HCC827_2_6Gy,na.rm=TRUE),sd(HCC827_2_8Gy,na.rm=TRUE)))
mean_HCC827_2["SEM"] <- mean_HCC827_2[,4]/sqrt(6)
mean_HCC827_2["N/N0"] <- mean_HCC827_2[,3]/mean_HCC827_2[,2]
mean_HCC827_2["PE"] <- mean_HCC827_2[1,3]/mean_HCC827_2[1,2]
mean_HCC827_2["SF"] <- mean_HCC827_2["N_mean"]/(mean_HCC827_2["PE"]*mean_HCC827_2["N0"])

HCC827_2["SF"]<- as.numeric(unlist(HCC827_2["N"]))/(mean_HCC827_2[1,7]*as.numeric(unlist(HCC827_2["N0"])))
####III HCC827 - Extracting Data####
HCC827_3 <- HCC827[2:31,13:15]
colnames(HCC827_3)<- c("Dose","N0","N")
HCC827_3["N0/N1"]<-as.numeric(HCC827_3[,3])/as.numeric(HCC827_3[,2])
HCC827_3["Plating effeciency"]<- (HCC827_3[,4]/mean(as.numeric(HCC827_3[1:6,4])))
cells_HCC827_3 <- data.frame("N0"=as.numeric(HCC827_3[c(1,7,13,19,25),2]))

HCC827_3_0Gy <- (as.numeric(HCC827_3[1:6,3:3]))
HCC827_3_2Gy <- (as.numeric(HCC827_3[7:12,3:3]))
HCC827_3_4Gy <- (as.numeric(HCC827_3[13:18,3:3]))
HCC827_3_6Gy <- (as.numeric(HCC827_3[19:24,3:3]))
HCC827_3_8Gy <- as.numeric(HCC827_3[25:30,3:3])
PE_HCC827_3<- mean(mean(HCC827_3_0Gy)/as.numeric(HCC827_3[1,2]))
mean_HCC827_3 <- data.frame("dose_Gy"=c(0,2,4,6,8),
                           "N0"=c(cells_HCC827_3),
                           "N_mean"= c(mean(HCC827_3_0Gy,na.rm=TRUE),mean(HCC827_3_2Gy,na.rm=TRUE),mean(HCC827_3_4Gy,na.rm=TRUE),mean(HCC827_3_6Gy,na.rm=TRUE),mean(HCC827_3_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC827_3_0Gy,na.rm=TRUE),sd(HCC827_3_2Gy,na.rm=TRUE),sd(HCC827_3_4Gy,na.rm=TRUE),sd(HCC827_3_6Gy,na.rm=TRUE),sd(HCC827_3_8Gy,na.rm=TRUE)))
mean_HCC827_3["SEM"] <- mean_HCC827_3[,4]/sqrt(6)
mean_HCC827_3["N/N0"] <- mean_HCC827_3[,3]/mean_HCC827_3[,2]
mean_HCC827_3["PE"] <- mean_HCC827_3[1,3]/mean_HCC827_3[1,2]
mean_HCC827_3["SF"] <- mean_HCC827_3["N_mean"]/(mean_HCC827_3["PE"]*mean_HCC827_3["N0"])

HCC827_3["SF"]<- as.numeric(unlist(HCC827_3["N"]))/(mean_HCC827_3[1,7]*as.numeric(unlist(HCC827_3["N0"])))
#####binding HCC827 results####
all_HCC827 <- cbind(mean_HCC827_1[,c(1:3,8)],mean_HCC827_2[,c(3,8)],mean_HCC827_3[,c(3,8)])
colnames(all_HCC827)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3")
all_HCC827["SF_mean"]<- rowMeans(all_HCC827[,c(4,6,8)],na.rm=TRUE)
all_HCC827["SF_sd"]<- apply(all_HCC827[,c(4,6,8)], 1, sd,na.rm=TRUE)
all_HCC827["SF_SEM"]<- all_HCC827["SF_sd"]/sqrt(3)
CFA_all_HCC827<-all_HCC827
####I HCC78 - Extracting Data####
HCC78<- read.csv("CFA_1_8Gy_HCC78.csv")
HCC78_1 <- HCC78[2:12,1:3]
colnames(HCC78_1)<- c("Dose","N0","N")
HCC78_1["N0/N1"]<-as.numeric(HCC78_1[,3])/as.numeric(HCC78_1[,2])
HCC78_1["Plating effeciency"]<- (HCC78_1[,4]/mean(as.numeric(HCC78_1[1:6,4])))
cells_HCC78_1 <- data.frame("N0"=as.numeric(HCC78_1[c(1,4,7,10),2]))

HCC78_1_0Gy <- (as.numeric(HCC78_1[1:3,3]))
HCC78_1_4Gy <- (as.numeric(HCC78_1[4:6,3]))
HCC78_1_6Gy <- (as.numeric(HCC78_1[7:9,3]))
HCC78_1_8Gy <- as.numeric(HCC78_1[10:12,3])
PE_HCC78_1<- mean(mean(HCC78_1_0Gy)/as.numeric(HCC78_1[1,2]))
mean_HCC78_1 <- data.frame("dose_Gy"=c(0,4,6,8),
                           "N0"=c(cells_HCC78_1),
                           "N_mean"= c(mean(HCC78_1_0Gy,na.rm=TRUE),mean(HCC78_1_4Gy,na.rm=TRUE),mean(HCC78_1_6Gy,na.rm=TRUE),mean(HCC78_1_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC78_1_0Gy,na.rm=TRUE),sd(HCC78_1_4Gy,na.rm=TRUE),sd(HCC78_1_6Gy,na.rm=TRUE),sd(HCC78_1_8Gy,na.rm=TRUE)))
mean_HCC78_1["SEM"] <- mean_HCC78_1[,4]/sqrt(3)
mean_HCC78_1["N/N0"] <- mean_HCC78_1[,3]/mean_HCC78_1[,2]
mean_HCC78_1["PE"] <- mean_HCC78_1[1,3]/mean_HCC78_1[1,2]
mean_HCC78_1["SF"] <- mean_HCC78_1["N_mean"]/(mean_HCC78_1["PE"]*mean_HCC78_1["N0"])

HCC78_1["SF"]<- as.numeric(unlist(HCC78_1["N"]))/(mean_HCC78_1[1,7]*as.numeric(unlist(HCC78_1["N0"])))
####II HCC78 - Extracting Data####
HCC78_2 <- HCC78[2:31,7:9]
colnames(HCC78_2)<- c("Dose","N0","N")
HCC78_2["N0/N1"]<-as.numeric(HCC78_2[,3])/as.numeric(HCC78_2[,2])
HCC78_2["Plating effeciency"]<- (HCC78_2[,4]/mean(as.numeric(HCC78_2[1:6,4])))
cells_HCC78_2 <- data.frame("N0"=as.numeric(HCC78_2[c(1,4,7,10),2]))

HCC78_2_0Gy <- (as.numeric(HCC78_2[1:3,3:3]))
HCC78_2_4Gy <- (as.numeric(HCC78_2[4:6,3:3]))
HCC78_2_6Gy <- (as.numeric(HCC78_2[7:9,3:3]))
HCC78_2_8Gy <- as.numeric(HCC78_2[10:12,3:3])
PE_HCC78_2<- mean(mean(HCC78_2_0Gy)/as.numeric(HCC78_2[1,2]))
mean_HCC78_2 <- data.frame("dose_Gy"=c(0,4,6,8),
                           "N0"=c(cells_HCC78_2),
                           "N_mean"= c(mean(HCC78_2_0Gy,na.rm=TRUE),mean(HCC78_2_4Gy,na.rm=TRUE),mean(HCC78_2_6Gy,na.rm=TRUE),mean(HCC78_2_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC78_2_0Gy,na.rm=TRUE),sd(HCC78_2_4Gy,na.rm=TRUE),sd(HCC78_2_6Gy,na.rm=TRUE),sd(HCC78_2_8Gy,na.rm=TRUE)))
mean_HCC78_2["SEM"] <- mean_HCC78_2[,4]/sqrt(3)
mean_HCC78_2["N/N0"] <- mean_HCC78_2[,3]/mean_HCC78_2[,2]
mean_HCC78_2["PE"] <- mean_HCC78_2[1,3]/mean_HCC78_2[1,2]
mean_HCC78_2["SF"] <- mean_HCC78_2["N_mean"]/(mean_HCC78_2["PE"]*mean_HCC78_2["N0"])

HCC78_2["SF"]<- as.numeric(unlist(HCC78_2["N"]))/(mean_HCC78_2[1,7]*as.numeric(unlist(HCC78_2["N0"])))
####III HCC78 - Extracting Data####
HCC78_3 <- HCC78[2:31,13:15]
colnames(HCC78_3)<- c("Dose","N0","N")
HCC78_3["N0/N1"]<-as.numeric(HCC78_3[,3])/as.numeric(HCC78_3[,2])
HCC78_3["Plating effeciency"]<- (HCC78_3[,4]/mean(as.numeric(HCC78_3[1:6,4])))
cells_HCC78_3 <- data.frame("N0"=as.numeric(HCC78_3[c(1,4,7,10),2]))

HCC78_3_0Gy <- (as.numeric(HCC78_3[1:3,3:3]))
HCC78_3_4Gy <- (as.numeric(HCC78_3[4:6,3:3]))
HCC78_3_6Gy <- (as.numeric(HCC78_3[7:9,3:3]))
HCC78_3_8Gy <- as.numeric(HCC78_3[10:12,3:3])
PE_HCC78_3<- mean(mean(HCC78_3_0Gy)/as.numeric(HCC78_3[1,2]))
mean_HCC78_3 <- data.frame("dose_Gy"=c(0,4,6,8),
                           "N0"=c(cells_HCC78_3),
                           "N_mean"= c(mean(HCC78_3_0Gy,na.rm=TRUE),mean(HCC78_3_4Gy,na.rm=TRUE),mean(HCC78_3_6Gy,na.rm=TRUE),mean(HCC78_3_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC78_3_0Gy,na.rm=TRUE),sd(HCC78_3_4Gy,na.rm=TRUE),sd(HCC78_3_6Gy,na.rm=TRUE),sd(HCC78_3_8Gy,na.rm=TRUE)))
mean_HCC78_3["SEM"] <- mean_HCC78_3[,4]/sqrt(3)
mean_HCC78_3["N/N0"] <- mean_HCC78_3[,3]/mean_HCC78_3[,2]
mean_HCC78_3["PE"] <- mean_HCC78_3[1,3]/mean_HCC78_3[1,2]
mean_HCC78_3["SF"] <- mean_HCC78_3["N_mean"]/(mean_HCC78_3["PE"]*mean_HCC78_3["N0"])

HCC78_3["SF"]<- as.numeric(unlist(HCC78_3["N"]))/(mean_HCC78_3[1,7]*as.numeric(unlist(HCC78_3["N0"])))
####IV HCC78 - Extracting Data####
HCC78_4 <- HCC78[2:31,13:15]
colnames(HCC78_4)<- c("Dose","N0","N")
HCC78_4["N0/N1"]<-as.numeric(HCC78_4[,3])/as.numeric(HCC78_4[,2])
HCC78_4["Plating effeciency"]<- (HCC78_4[,4]/mean(as.numeric(HCC78_4[1:6,4])))
cells_HCC78_4 <- data.frame("N0"=as.numeric(HCC78_4[c(1,4,7,10),2]))

HCC78_4_0Gy <- (as.numeric(HCC78_4[1:3,3]))
HCC78_4_4Gy <- (as.numeric(HCC78_4[4:6,3]))
HCC78_4_6Gy <- (as.numeric(HCC78_4[7:9,3]))
HCC78_4_8Gy <- as.numeric(HCC78_4[10:12,3])
PE_HCC78_4<- mean(mean(HCC78_4_0Gy)/as.numeric(HCC78_4[1,2]))
mean_HCC78_4 <- data.frame("dose_Gy"=c(0,4,6,8),
                           "N0"=c(cells_HCC78_4),
                           "N_mean"= c(mean(HCC78_4_0Gy,na.rm=TRUE),mean(HCC78_4_4Gy,na.rm=TRUE),mean(HCC78_4_6Gy,na.rm=TRUE),mean(HCC78_4_8Gy,na.rm=TRUE)),
                           "N_SD"= c(sd(HCC78_4_0Gy,na.rm=TRUE),sd(HCC78_4_4Gy,na.rm=TRUE),sd(HCC78_4_6Gy,na.rm=TRUE),sd(HCC78_4_8Gy,na.rm=TRUE)))
mean_HCC78_4["SEM"] <- mean_HCC78_4[,4]/sqrt(3)
mean_HCC78_4["N/N0"] <- mean_HCC78_4[,3]/mean_HCC78_4[,2]
mean_HCC78_4["PE"] <- mean_HCC78_4[1,3]/mean_HCC78_4[1,2]
mean_HCC78_4["SF"] <- mean_HCC78_4["N_mean"]/(mean_HCC78_4["PE"]*mean_HCC78_4["N0"])

HCC78_4["SF"]<- as.numeric(unlist(HCC78_4["N"]))/(mean_HCC78_4[1,7]*as.numeric(unlist(HCC78_4["N0"])))
####binding HCC78 results####
all_HCC78 <- cbind(mean_HCC78_1[,c(1:3,8)],mean_HCC78_2[,c(3,8)],mean_HCC78_3[,c(3,8)],mean_HCC78_4[,c(3,8)])
colnames(all_HCC78)<- c("Dose_Gy","N0","N_mean_rep1","SF_rep1","N_mean_rep2","SF_rep2","N_mean_rep3","SF_rep3","N_mean_rep4","SF_rep4")
all_HCC78["SF_mean"]<- rowMeans(all_HCC78[,c(4,6,8,10)],na.rm=TRUE)
all_HCC78["SF_sd"]<- apply(all_HCC78[,c(4,6,8,10)], 1, sd,na.rm=TRUE)
all_HCC78["SF_SEM"]<- all_HCC78["SF_sd"]/sqrt(4)
CFA_all_HCC78<-all_HCC78
###Plotting data#####
setwd(save_wd)
jpeg(filename = "CFA_8Gy_base.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1),ncol=3,nrow = 1,byrow = TRUE))
leg.txt <- c("HCC44","NCI-H23","HCC827","A549","HCC78")
par(mar = c(10, 14, 6, 4),xpd=TRUE)
plot(x=dose, y=(all_HCC44[,10]),cex=4,axes=FALSE,cex.axis=2, lwd=3,pch=22,log="y",bg="#3288BD", 
     ylab=" ", xlab="", lty=1.5,cex.main=4,cex.lab=4, #main="CFA of HCC44 and H23",
     cex.main=2,mgp=c(6, 0,0),ylim=c(0.001,50),xlim=c(0,9))
points(x=dose[1:2], y=c(mean_H23_1[1,8],0.001),lwd=3,ylim=c(0.01,10), cex=4, pch=23, bg="#66C2A5")
points(x=dose, y=(all_HCC827[,8]), lwd=3,ylim=c(0.01,10), cex=4, pch=24, bg="#FDAE61")
points(x=dose, y=(all_A549[,10]), lwd=3,ylim=c(0.01,10), cex=4, pch=25, bg="#D53E4F")
points(x=c(0,2,4,8), y=all_HCC78[,10], lwd=3,ylim=c(0.01,10), cex=4, pch=25, bg="#5E4FA2")
lines(x=dose, y=(all_HCC44[,10]),col="black",lwd=3,lty=3)
lines(x=dose[1:2], y=c(mean_H23_1[1,8],0.001),col="black",lwd=3,lty=3)
lines(x=dose, y=(all_HCC827[,8]),col="black",lwd=3,lty=3)
lines(x=c(0,2,4,8), y=all_HCC78[,10],col="black",lwd=3,lty=3)
lines(x=dose[1:4], y=(all_A549[1:4,10]),col="black",lwd=3,lty=3)
lines(x=c(6,7), y=c(mean(as.numeric(all_A549[4,2:4])),0),col="black",lwd=3)
box(which="plot",lty="solid",lwd=3)
title( ylab="Surviving Fraction", cex.lab=3,cex.main=6,line=9,adj=0.75)
title( xlab="Dose [Gy]", cex.lab=3,cex.main=6,line=7)
axis(1,lwd=3,at = c(0,4,8),mgp=c(4, 3.5,0),lwd.ticks = 3,cex=2, cex.axis=3,tcl = -1.5)
axis(2,lwd=3,at = c(1, 0.1, 0.01,0.001),label=c("1","","",parse(text="10^-3")),las=2,mgp=c(4, 2,0),lwd.ticks =3,cex=6, cex.axis=3,tcl = -1.5)
legend(c(3.5,9),c(40,60),leg.txt, ncol=2, pch=c(22,23,24,25,24,25), bg = NA , box.col="white",adj=c(0,0.25),pt.cex=3,pt.lwd=3, cex = 2.75,pt.bg=c("#3288BD","#66C2A5","#FDAE61","#D53E4F","#5E4FA2"),box.lwd=2,y.intersp =1.1)
minor.tick(nx =2/1,ny=0,tick.ratio =1.5, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3)) 
par(mar = c(10, 2, 6, 2),xpd=TRUE)
dev.off()
setwd(source_wd)



