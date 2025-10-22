#####H23 - importing Data####
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
names <- c("0h","24h","48h","72h","96h","120h")
col_names <- c("control_1","control_2","control_3","2Gy_1","2Gy_2","2Gy_2")
conc <- c(0,0.001,0.1, 1, 10, 100)
setwd("C:/Users/annar/Desktop/Codes+Data/functional_assays/MTS assays/Exp2_IC50_AZD1390/processed")
#####loading raw data####
H23_1 <- read.csv("IC50_AZD1390_NCI-H23_rep1.csv",sep="\t",dec=",")
H23_1 <- data.frame(lapply(H23_1, as.numeric))
H23_2 <- read.csv("IC50_AZD1390_NCI-H23_rep2.csv",sep="\t",dec=",")
H23_2 <- data.frame(lapply(H23_2, as.numeric))
H23_3 <- read.csv("IC50_AZD1390_NCI-H23_rep3.csv",sep="\t",dec=",")
H23_3 <- data.frame(lapply(H23_3, as.numeric))

HCC44_1 <- read.csv("IC50_AZD1390_HCC44_rep1.csv",sep="\t",dec=",")
HCC44_1 <- data.frame(lapply(HCC44_1, as.numeric))
HCC44_2 <- read.csv("IC50_AZD1390_HCC44_rep2.csv",sep="\t",dec=",")
HCC44_1 <- data.frame(lapply(HCC44_1, as.numeric))
HCC44_3 <- read.csv("IC50_AZD1390_HCC44_rep3.csv",sep="\t",dec=",")
HCC44_3 <- data.frame(lapply(HCC44_3, as.numeric))
#####Splitting data - HCC44 biological replicate 1#####
HCC44_c_1.1 <- t(HCC44_1[2:7,2:11])
HCC44_c_1.1 <- range01(HCC44_c_1.1)
colnames(HCC44_c_1.1)<-conc

HCC44_c_1.2 <- t(HCC44_1[11:16,2:11])
HCC44_c_1.2 <- range01(HCC44_c_1.2)
colnames(HCC44_c_1.2)<-conc

HCC44_c_1.3 <- t(HCC44_1[20:25,2:11])
HCC44_c_1.3 <- range01(HCC44_c_1.3)
colnames(HCC44_c_1.3)<-conc
#####binding HCC44 replicate 1#####
calc_HCC44_c_1.1 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_1.1[,1],HCC44_c_1.1[,2],HCC44_c_1.1[,3],HCC44_c_1.1[,4],HCC44_c_1.1[,5],HCC44_c_1.1[,6]),
                               replicate= rep("Rep1", 60) )
calc_HCC44_c_1.2 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_1.2[,1],HCC44_c_1.2[,2],HCC44_c_1.2[,3],HCC44_c_1.2[,4],HCC44_c_1.2[,5],HCC44_c_1.2[,6]),
                               replicate= rep("Rep2",60) )
calc_HCC44_c_1.3 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_1.3[,1],HCC44_c_1.3[,2],HCC44_c_1.3[,3],HCC44_c_1.3[,4],HCC44_c_1.3[,5],HCC44_c_1.3[,6]),
                               replicate= rep("Rep3",60) )                               
calc_HCC44_c_1<- rbind(calc_HCC44_c_1.1,calc_HCC44_c_1.2,calc_HCC44_c_1.3)
#####modelling HCC44 replicate 1####
modell_HCC44_c_1 <- drm(value ~ conc, data = calc_HCC44_c_1, curveid = replicate, fct = LL.4())
ceofficient_HCC44_c_1 <-summary(modell_HCC44_c_1)
ic50_HCC44_c_1 <- rbind(ED(modell_HCC44_c_1, 50, interval = "delta"),ED(modell_HCC44_c_1, 10, interval = "delta"),ED(modell_HCC44_c_1, 90, interval = "delta"))
AIC_HCC44_c_1 <- AIC(modell_HCC44_c_1)
predict_HCC44_c_1 <- predict(modell_HCC44_c_1)
residuals_HCC44_c_1 <- calc_HCC44_c_1$value - predict_HCC44_c_1
ss_res <- sum(residuals_HCC44_c_1^2)
ss_tot <- sum((calc_HCC44_c_1$value - mean(calc_HCC44_c_1$value))^2)
r_squared_HCC44_c_1 <- 1 - (ss_res / ss_tot)
#####Splitting data - HCC44 biological replicate 2#####
HCC44_c_2.1 <- t(HCC44_2[2:7,2:11])
HCC44_c_2.1 <- range01(HCC44_c_2.1)
colnames(HCC44_c_2.1)<-conc

HCC44_c_2.2 <- t(HCC44_2[11:16,2:11])
HCC44_c_2.2 <- range01(HCC44_c_2.2)
colnames(HCC44_c_2.2)<-conc

HCC44_c_2.3 <- t(HCC44_2[20:25,2:11])
HCC44_c_2.3 <- range01(HCC44_c_2.3)
colnames(HCC44_c_2.3)<-conc
#####binding HCC44 replicate 2#####
calc_HCC44_c_2.1 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_2.1[,1],HCC44_c_2.1[,2],HCC44_c_2.1[,3],HCC44_c_2.1[,4],HCC44_c_2.1[,5],HCC44_c_2.1[,6]),
                               replicate= rep("Rep1", 60) )
calc_HCC44_c_2.2 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_2.2[,1],HCC44_c_2.2[,2],HCC44_c_2.2[,3],HCC44_c_2.2[,4],HCC44_c_2.2[,5],HCC44_c_2.2[,6]),
                               replicate= rep("Rep2",60) )
calc_HCC44_c_2.3 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_2.3[,1],HCC44_c_2.3[,2],HCC44_c_2.3[,3],HCC44_c_2.3[,4],HCC44_c_2.3[,5],HCC44_c_2.3[,6]),
                               replicate= rep("Rep3",60) )                               
calc_HCC44_c_2<- rbind(calc_HCC44_c_2.1,calc_HCC44_c_2.2,calc_HCC44_c_2.3)
#####modelling HCC44 replicate 2####
modell_HCC44_c_2 <- drm(value ~ conc, data = calc_HCC44_c_2, curveid = replicate, fct = LL.4())
ceofficient_HCC44_c_2 <-summary(modell_HCC44_c_2)
ic50_HCC44_c_2 <- rbind(ED(modell_HCC44_c_2, 50, interval = "delta"),ED(modell_HCC44_c_2, 10, interval = "delta"),ED(modell_HCC44_c_2, 90, interval = "delta"))
AIC_HCC44_c_2 <- AIC(modell_HCC44_c_2)
predict_HCC44_c_2 <- predict(modell_HCC44_c_2)
residuals_HCC44_c_2 <- calc_HCC44_c_2$value - predict_HCC44_c_2
ss_res <- sum(residuals_HCC44_c_2^2)
ss_tot <- sum((calc_HCC44_c_2$value - mean(calc_HCC44_c_2$value))^2)
r_squared_HCC44_c_2 <- 1 - (ss_res / ss_tot)
#####Splitting data - HCC44 biological replicate 3#####
HCC44_c_3.1 <- t(HCC44_3[2:7,2:11])
HCC44_c_3.1 <- range01(HCC44_c_3.1)
colnames(HCC44_c_3.1)<-conc

HCC44_c_3.2 <- t(HCC44_3[11:16,2:11])
HCC44_c_3.2 <- range01(HCC44_c_3.2)
colnames(HCC44_c_3.2)<-conc

HCC44_c_3.3 <- t(HCC44_3[20:25,2:11])
HCC44_c_3.3 <- range01(HCC44_c_3.3)
colnames(HCC44_c_3.3)<-conc
calc_HCC44_c_3<- rbind(calc_HCC44_c_3.1,calc_HCC44_c_3.2,calc_HCC44_c_3.3)
#####binding HCC44 replicate 2#####
calc_HCC44_c_3.1 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_3.1[,1],HCC44_c_3.1[,2],HCC44_c_3.1[,3],HCC44_c_3.1[,4],HCC44_c_3.1[,5],HCC44_c_3.1[,6]),
                               replicate= rep("Rep1", 60) )
calc_HCC44_c_3.2 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_3.2[,1],HCC44_c_3.2[,2],HCC44_c_3.2[,3],HCC44_c_3.2[,4],HCC44_c_3.2[,5],HCC44_c_3.2[,6]),
                               replicate= rep("Rep2",60) )
calc_HCC44_c_3.3 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                               value=c(HCC44_c_3.3[,1],HCC44_c_3.3[,2],HCC44_c_3.3[,3],HCC44_c_3.3[,4],HCC44_c_3.3[,5],HCC44_c_3.3[,6]),
                               replicate= rep("Rep3",60) )                               
calc_HCC44_c_3<- rbind(calc_HCC44_c_3.1,calc_HCC44_c_3.2,calc_HCC44_c_3.3)

#####modelling HCC44 replicate 3####
modell_HCC44_c_3 <- drm(value ~ conc, data = calc_HCC44_c_3, curveid = replicate, fct = LL.4())
ceofficient_HCC44_c_3 <-summary(modell_HCC44_c_3)
ic50_HCC44_c_3 <- rbind(ED(modell_HCC44_c_3, 50, interval = "delta"),ED(modell_HCC44_c_3, 10, interval = "delta"),ED(modell_HCC44_c_3, 90, interval = "delta"))
AIC_HCC44_c_3 <- AIC(modell_HCC44_c_3)
predict_HCC44_c_3 <- predict(modell_HCC44_c_3)
residuals_HCC44_c_3 <- calc_HCC44_c_3$value - predict_HCC44_c_3
ss_res <- sum(residuals_HCC44_c_3^2)
ss_tot <- sum((calc_HCC44_c_3$value - mean(calc_HCC44_c_3$value))^2)
r_squared_HCC44_c_3 <- 1 - (ss_res / ss_tot)
#####binding HCC44 all#####
calc_HCC44_c_all<- data.frame(conc=c(rep(c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),3)),
                              value=c(rowMeans(cbind(calc_HCC44_c_1.1[,2],calc_HCC44_c_1.2[,2],calc_HCC44_c_1.3[,2])),rowMeans(cbind(calc_HCC44_c_2.1[,2],calc_HCC44_c_2.2[,2],calc_HCC44_c_2.3[,2])),
                                      rowMeans(cbind(calc_HCC44_c_1.1[,2],calc_HCC44_c_1.2[,2],calc_HCC44_c_1.3[,2]))),
                              replicate= c(rep("Rep1",60),rep("Rep2",60),rep("Rep3",60) ))   
modell_HCC44_all <- drm(value ~ conc, data = calc_HCC44_c_all, curveid = replicate, fct = LL.4())
######Binding parameters HCC44####
coefficient_HCC44 <- cbind(ceofficient_HCC44_c_1 [["coefficients"]],rep("HCC44_rep1",12),ceofficient_HCC44_c_2 [["coefficients"]],rep("HCC44_rep2",12),ceofficient_HCC44_c_3 [["coefficients"]],rep("HCC44_rep3",12))
AIC_HCC44 <- rbind(c(AIC_HCC44_c_1,"HCC44_rep1"),c(AIC_HCC44_c_2,"HCC44_rep2"),c(AIC_HCC44_c_3,"HCC44_rep3"))
ic50_HCC44<- rbind(cbind(ic50_HCC44_c_1,rep("HCC44_rep1",9)),cbind(ic50_HCC44_c_2,rep("HCC44_rep2",9)),cbind(ic50_HCC44_c_3,rep("HCC44_rep3",9)))
r_squared_HCC44<- rbind(c(r_squared_HCC44_c_1,"HCC44_rep1"),c(r_squared_HCC44_c_2,"HCC44_rep2"),c(r_squared_HCC44_c_3,"HCC44_rep3"))
#####Splitting data - HCC44 biological replicate 1#####
H23_c_1.1 <- t(H23_1[2:7,2:11])
H23_c_1.1 <- range01(H23_c_1.1)
colnames(H23_c_1.1)<-conc

H23_c_1.2 <- t(H23_1[11:16,2:11])
H23_c_1.2 <- range01(H23_c_1.2)
colnames(H23_c_1.2)<-conc

H23_c_1.3 <- t(H23_1[20:25,2:11])
H23_c_1.3 <- range01(H23_c_1.3)
colnames(H23_c_1.3)<-conc

#####binding H23 replicate 1#####
calc_H23_c_1.1 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_1.1[,1],H23_c_1.1[,2],H23_c_1.1[,3],H23_c_1.1[,4],H23_c_1.1[,5],H23_c_1.1[,6]),
                             replicate= rep("Rep1", 60) )
calc_H23_c_1.2 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_1.2[,1],H23_c_1.2[,2],H23_c_1.2[,3],H23_c_1.2[,4],H23_c_1.2[,5],H23_c_1.2[,6]),
                             replicate= rep("Rep2",60) )
calc_H23_c_1.3 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_1.3[,1],H23_c_1.3[,2],H23_c_1.3[,3],H23_c_1.3[,4],H23_c_1.3[,5],H23_c_1.3[,6]),
                             replicate= rep("Rep3",60) )                               
calc_H23_c_1<- rbind(calc_H23_c_1.1,calc_H23_c_1.2,calc_H23_c_1.3)
#####modelling H23 replicate 1####
modell_H23_c_1 <- drm(value ~ conc, data = calc_H23_c_1, curveid = replicate, fct = LL.4())
ceofficient_H23_c_1 <-summary(modell_H23_c_1)
ic50_H23_c_1 <- (rbind(ED(modell_H23_c_1, 50, interval = "delta"),ED(modell_H23_c_1, 10, interval = "delta"),ED(modell_H23_c_1, 90, interval = "delta")))
AIC_H23_c_1 <- AIC(modell_H23_c_1)
predict_H23_c_1 <- predict(modell_H23_c_1)
residuals_H23_c_1 <- calc_H23_c_1$value - predict_H23_c_1
ss_res <- sum(residuals_H23_c_1^2)
ss_tot <- sum((calc_H23_c_1$value - mean(calc_H23_c_1$value))^2)
r_squared_H23_c_1 <- 1 - (ss_res / ss_tot)
#####Splitting data - H23 biological replicate 2#####
H23_c_2.1 <- t(H23_2[2:7,2:11])
H23_c_2.1 <- range01(H23_c_2.1)
colnames(H23_c_2.1)<-conc

H23_c_2.2 <- t(H23_2[11:16,2:11])
H23_c_2.2 <- range01(H23_c_2.2)
colnames(H23_c_2.2)<-conc

H23_c_2.3 <- t(H23_2[20:25,2:11])
H23_c_2.3 <- range01(H23_c_2.3)
colnames(H23_c_2.3)<-conc
#####binding H23 replicate 1#####
calc_H23_c_2.1 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_2.1[,1],H23_c_2.1[,2],H23_c_2.1[,3],H23_c_2.1[,4],H23_c_2.1[,5],H23_c_2.1[,6]),
                             replicate= rep("Rep1", 60) )
calc_H23_c_2.2 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_2.2[,1],H23_c_2.2[,2],H23_c_2.2[,3],H23_c_2.2[,4],H23_c_2.2[,5],H23_c_2.2[,6]),
                             replicate= rep("Rep2",60) )
calc_H23_c_2.3 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_2.3[,1],H23_c_2.3[,2],H23_c_2.3[,3],H23_c_2.3[,4],H23_c_2.3[,5],H23_c_2.3[,6]),
                             replicate= rep("Rep3",60) )                               
calc_H23_c_2<- rbind(calc_H23_c_2.1,calc_H23_c_2.2,calc_H23_c_2.3)

#####modelling H23 replicate 2####
modell_H23_c_2 <- drm(value ~ conc, data = calc_H23_c_2, curveid = replicate, fct = LL.4())
ceofficient_H23_c_2 <-summary(modell_H23_c_2)
ic50_H23_c_2 <- (rbind(ED(modell_H23_c_2, 50, interval = "delta",na.rm=TRUE),ED(modell_H23_c_2, 10, interval = "delta",na.rm=TRUE),ED(modell_H23_c_2, 90, interval = "delta",na.rm=TRUE)))
AIC_H23_c_2 <- AIC(modell_H23_c_2)
predict_H23_c_2 <- predict(modell_H23_c_2)
residuals_H23_c_2 <- calc_H23_c_2$value - predict_H23_c_2
ss_res <- sum(residuals_H23_c_2^2)
ss_tot <- sum((calc_H23_c_2$value - mean(calc_H23_c_2$value))^2)
r_squared_H23_c_2 <- 1 - (ss_res / ss_tot)
#####Splitting data - H23 biological replicate 3#####
H23_c_3.1 <- t(H23_3[2:7,2:11])
H23_c_3.1 <- range01(H23_c_3.1)
colnames(H23_c_3.1)<-conc

H23_c_3.2 <- t(H23_3[11:16,2:11])
H23_c_3.2 <- range01(H23_c_3.2)
colnames(H23_c_3.2)<-conc

H23_c_3.3 <- t(H23_3[20:25,2:11])
H23_c_3.3 <- range01(H23_c_3.3)
colnames(H23_c_3.3)<-conc
#####binding H23 replicate 3#####
calc_H23_c_3.1 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_3.1[,1],H23_c_3.1[,2],H23_c_3.1[,3],H23_c_3.1[,4],H23_c_3.1[,5],H23_c_3.1[,6]),
                             replicate= rep("Rep1", 60) )
calc_H23_c_3.2 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_3.2[,1],H23_c_3.2[,2],H23_c_3.2[,3],H23_c_3.2[,4],H23_c_3.2[,5],H23_c_3.2[,6]),
                             replicate= rep("Rep2",60) )
calc_H23_c_3.3 <- data.frame(conc=c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),
                             value=c(H23_c_3.3[,1],H23_c_3.3[,2],H23_c_3.3[,3],H23_c_3.3[,4],H23_c_3.3[,5],H23_c_3.3[,6]),
                             replicate= rep("Rep3",60) )                               
calc_H23_c_3<- rbind(calc_H23_c_3.1,calc_H23_c_3.2,calc_H23_c_3.3)
#####modelling H23 replicate 3####
#modell_H23_c_3 <- drm(value ~ conc, data = calc_H23_c_3, curveid = replicate, fct = LL.4())
#ceofficient_H23_c_3 <-summary(modell_H23_c_3)
#ic50_H23_c_3 <- rbind(ED(modell_H23_c_3, 50, interval = "delta"),ED(modell_H23_c_3, 10, interval = "delta"),ED(H23_H23_c_3, 90, interval = "delta"))
#AIC_H23_c_3 <- AIC(H23_H23_c_3)
#predict_H23_c_3 <- predict(modell_H23_c_3)
#residuals_H23_c_3 <- calc_H23_c_3$value - predict_H23_c_3
#ss_res <- sum(residuals_H23_c_3^2)
#ss_tot <- sum((calc_H23_c_3$value - mean(calc_H23_c_3$value))^2)
#r_squared_H23_c_3 <- 1 - (ss_res / ss_tot)
#####binding HCC44 all#####
calc_H23_c_all<- data.frame(conc=c(rep(c(rep(0,10),rep(0.01,10),rep(0.1,10),rep(1,10),rep(10,10),rep(100,10)),3)),
                            value=c(rowMeans(cbind(calc_H23_c_1.1[,2],calc_H23_c_1.2[,2],calc_H23_c_1.3[,2])),rowMeans(cbind(calc_H23_c_2.1[,2],calc_H23_c_2.2[,2],calc_H23_c_2.3[,2])),
                                    rowMeans(cbind(calc_H23_c_1.1[,2],calc_H23_c_1.2[,2],calc_H23_c_1.3[,2]))),
                            replicate= c(rep("Rep1",60),rep("Rep2",60),rep("Rep3",60) ))   
modell_H23_all <- drm(value ~ conc, data = calc_H23_c_all, curveid = replicate, fct = LL.4())
#####
IC50_mean<- data.frame(HCC44_mean=colMeans(rbind(HCC44_c_1.1,HCC44_c_1.2,HCC44_c_1.3,HCC44_c_2.1,HCC44_c_2.2,HCC44_c_2.3,HCC44_c_3.1,HCC44_c_3.2,HCC44_c_3.3)),
                       HCC44_sd=apply(rbind(HCC44_c_1.1,HCC44_c_1.2,HCC44_c_1.3,HCC44_c_2.1,HCC44_c_2.2,HCC44_c_2.3,HCC44_c_3.1,HCC44_c_3.2,HCC44_c_3.3), 2, sd, na.rm=TRUE),
                       H23=colMeans(rbind(H23_c_1.1,H23_c_1.2,H23_c_1.3,H23_c_2.1,H23_c_2.2,H23_c_2.3,H23_c_3.1,H23_c_3.2,H23_c_3.3)),
                       H23_sd=apply(rbind(H23_c_1.1,H23_c_1.2,H23_c_1.3,H23_c_2.1,H23_c_2.2,H23_c_2.3,H23_c_3.1,H23_c_3.2,H23_c_3.3), 2, sd, na.rm=TRUE))
#####plotting models#####
setwd(save_wd)
jpeg(filename = "MTS_IC50_AZD1390.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,1),ncol=3,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 6, 4),xpd=TRUE)
leg.txt<-c("HCC44","NCI-H23")
plot(c(0:5),IC50_mean[,3],pch=21,col="black",lwd=2,xlim=c(0,5.5), bg="#66C2A5", axes=FALSE,ylim = c(0,1.4),
     ann=FALSE,cex=4,cex.axis = 3,cex.lab=4,mgp=c(7, 2,0))
lines(c(0:5),IC50_mean[,3],col="black",lwd=3,lty=3)
points(c(0:5),IC50_mean[,1],pch=21,col="black",lwd=2, bg="#3288BD",cex=4)
lines(c(0:5),IC50_mean[,1],col="black",lwd=3,lty=3)
box(type="solid",lwd=3)
title( xlab="AZD1390[µm]", cex.lab=3.5,cex.main=6,line=7)
title( ylab="Viabilty [a.u.]", cex.lab=3.5,cex.main=6,line=9)
axis(1,lwd=3,at=c(0,2,4),labels  = c(0,1,100),lwd.ticks = 3,cex=2, cex.axis=3,tcl = -1.5,mgp=c(3, 4,0))
axis(2,lwd=3,at=c(0,0.5,1),labels  = c(0,50,100),las=2,lwd.ticks = 3,cex=2, cex.axis=3,tcl = -1.5,mgp=c(3, 3,0))
legend(x=c(4,5.5), y=c(1.1,1.4),ncol=1, box.col="white",bty = "n", bg = NA,pt.lwd=2,adj=c(0,0.25),leg.txt,cex=2.5, pch=c(21), pt.bg=c("#3288BD","#66C2A5"),box.lwd = 2,pt.cex=4,yjust=0.5,xjust=0,y.intersp = 0.9,x.intersp = 0.9)
minor.tick(nx =2/1,ny=2/1,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
dev.off()
setwd(source_wd)