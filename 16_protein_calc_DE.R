#####Proteome - Binding - HCC44 0h####
raw_HCC44_0h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[17],proteome_data[40],proteome_data[67],proteome_data[56])
raw_HCC44_0h_2Gy$mean<-rowMeans(raw_HCC44_0h_2Gy[,4:7])
raw_HCC44_0h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[18],proteome_data[41],proteome_data[68],proteome_data[57])
raw_HCC44_0h_c$mean<-rowMeans(raw_HCC44_0h_c[,4:7])

HCC44_0h <- cbind(raw_HCC44_0h_c,raw_HCC44_0h_2Gy[,4:7])
#####Proteome - Calculating Diff and pvalue - HCC44 0h####
HCC44_0h <- raw_HCC44_0h_2Gy[,1:3]
HCC44_0h$mean_c <- rowMeans(log2(raw_HCC44_0h_c[4:7]), na.rm=TRUE)
HCC44_0h$SD_c <- apply(log2(raw_HCC44_0h_c[4:7]), 1, sd, na.rm=TRUE)
HCC44_0h$mean_2Gy <- rowMeans(log2(raw_HCC44_0h_2Gy[4:7]), na.rm=TRUE)
HCC44_0h$SD_2Gy <- apply(log2(raw_HCC44_0h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pHCC44_0h<- ttest_protti(mean1 = HCC44_0h$mean_2Gy,
                                   mean2 = HCC44_0h$mean_c,
                                   sd1 = HCC44_0h$SD_2Gy,
                                   sd2 = HCC44_0h$SD_c,
                                   n1 = 4,
                                   n2 = 4,
                                   log_values=TRUE)

fdrs_pHCC44_0h<-data.frame(p.adjust(res.ttest_pHCC44_0h$pval, method="BH"))

ready_pHCC44_0h<-cbind(HCC44_0h[,1:3],res.ttest_pHCC44_0h,fdrs_pHCC44_0h,proteome_data[,c(2,95:96)])
colnames(ready_pHCC44_0h)[8] <- 'adj_pVal'

#####Proteome - Binding - HCC44 2h####
raw_HCC44_2h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[21],proteome_data[31],proteome_data[69],proteome_data[58])
raw_HCC44_2h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[22],proteome_data[32],proteome_data[70],proteome_data[59])

HCC44_2h <- cbind(raw_HCC44_2h_c,raw_HCC44_2h_2Gy[,4:7])
#####Proteome - Calculating Diff and pvalue - HCC44 2h####
HCC44_2h <- raw_HCC44_2h_2Gy[,1:3]
HCC44_2h$mean_c <- rowMeans(log2(raw_HCC44_2h_c[4:7]), na.rm=TRUE)
HCC44_2h$SD_c <- apply(log2(raw_HCC44_2h_c[4:7]), 1, sd, na.rm=TRUE)
HCC44_2h$mean_2Gy <- rowMeans(log2(raw_HCC44_2h_2Gy[4:7]), na.rm=TRUE)
HCC44_2h$SD_2Gy <- apply(log2(raw_HCC44_2h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pHCC44_2h<- ttest_protti(mean1 = HCC44_2h$mean_2Gy,
                                   mean2 = HCC44_2h$mean_c,
                                   sd1 = HCC44_2h$SD_2Gy,
                                   sd2 = HCC44_2h$SD_c,
                                   n1 = 4,
                                   n2 = 4,
                                   log_values=TRUE)

fdrs_pHCC44_2h<-data.frame(p.adjust(res.ttest_pHCC44_2h$pval, method="BH"))

ready_pHCC44_2h<-cbind(HCC44_2h[,1:3],res.ttest_pHCC44_2h,fdrs_pHCC44_2h,proteome_data[,c(2,95:96)])
colnames(ready_pHCC44_2h)[8] <- 'adj_pVal'
#####Proteome - Binding - HCC44 6h#####
raw_HCC44_6h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[34],proteome_data[33],proteome_data[71],proteome_data[60])
raw_HCC44_6h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[35],proteome_data[42],proteome_data[79],proteome_data[61])

HCC44_6h <- cbind(raw_HCC44_6h_c,raw_HCC44_6h_2Gy[,4:7])
#####Proteome - Calculating Diff and pvalue - HCC44 6h####
HCC44_6h <- raw_HCC44_6h_2Gy[,1:3]
HCC44_6h$mean_c <- rowMeans(log2(raw_HCC44_6h_c[4:7]), na.rm=TRUE)
HCC44_6h$SD_c <- apply(log2(raw_HCC44_6h_c[4:7]), 1, sd, na.rm=TRUE)
HCC44_6h$mean_2Gy <- rowMeans(log2(raw_HCC44_6h_2Gy[4:7]), na.rm=TRUE)
HCC44_6h$SD_2Gy <- apply(log2(raw_HCC44_6h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pHCC44_6h<- ttest_protti(mean1 = HCC44_6h$mean_2Gy,
                                   mean2 = HCC44_6h$mean_c,
                                   sd1 = HCC44_6h$SD_2Gy,
                                   sd2 = HCC44_6h$SD_c,
                                   n1 = 4,
                                   n2 = 4,
                                   log_values=TRUE)

fdrs_pHCC44_6h<-data.frame(p.adjust(res.ttest_pHCC44_6h$pval, method="BH"))

ready_pHCC44_6h<-cbind(HCC44_6h[,1:3],res.ttest_pHCC44_6h,fdrs_pHCC44_6h,proteome_data[,c(2,95:96)])
colnames(ready_pHCC44_6h)[8] <- 'adj_pVal'

#####Proteome - Binding - HCC44 24h####
raw_HCC44_24h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[36],proteome_data[43],proteome_data[80],proteome_data[62])
raw_HCC44_24h_2Gy$mean<-rowMeans(raw_HCC44_24h_2Gy[,4:7])

raw_HCC44_24h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[37],proteome_data[44],proteome_data[81],proteome_data[63])
raw_HCC44_24h_c$mean<-rowMeans(raw_HCC44_24h_c[,4:7])

HCC44_24h <- cbind(raw_HCC44_24h_c,raw_HCC44_24h_2Gy[,4:7])
#####Proteome - Calculating Diff and pvalue - HCC44 24h#####
HCC44_24h <- raw_HCC44_24h_2Gy[,1:3]
HCC44_24h$mean_c <- rowMeans(log2(raw_HCC44_24h_c[4:7]), na.rm=TRUE)
HCC44_24h$SD_c <- apply(log2(raw_HCC44_24h_c[4:7]), 1, sd, na.rm=TRUE)
HCC44_24h$mean_2Gy <- rowMeans(log2(raw_HCC44_24h_2Gy[4:7]), na.rm=TRUE)
HCC44_24h$SD_2Gy <- apply(log2(raw_HCC44_24h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pHCC44_24h<- ttest_protti(mean1 = HCC44_24h$mean_2Gy,
                                    mean2 = HCC44_24h$mean_c,
                                    sd1 = HCC44_24h$SD_2Gy,
                                    sd2 = HCC44_24h$SD_c,
                                    n1 = 4,
                                    n2 = 4,
                                    log_values=TRUE)

fdrs_pHCC44_24h<-data.frame(p.adjust(res.ttest_pHCC44_24h$pval, method="BH"))

ready_pHCC44_24h<-cbind(HCC44_24h[,1:3],res.ttest_pHCC44_24h,fdrs_pHCC44_24h,proteome_data[,c(2,95:96)])
colnames(ready_pHCC44_24h)[8] <- 'adj_pVal'
#####Proteome - Binding - HCC44 48h####
raw_HCC44_48h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[38],proteome_data[45],proteome_data[82],proteome_data[64])
raw_HCC44_48h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[39],proteome_data[46],proteome_data[83],proteome_data[65])

HCC44_48h <- cbind(raw_HCC44_48h_c,raw_HCC44_48h_2Gy[,4:7])
#####Proteome - calculating Diff and pvalue - HCC44 48h####
HCC44_48h <- raw_HCC44_48h_2Gy[,1:3]
HCC44_48h$mean_c <- rowMeans(log2(raw_HCC44_48h_c[4:7]), na.rm=TRUE)
HCC44_48h$SD_c <- apply(log2(raw_HCC44_48h_c[4:7]), 1, sd, na.rm=TRUE)
HCC44_48h$mean_2Gy <- rowMeans(log2(raw_HCC44_48h_2Gy[4:7]), na.rm=TRUE)
HCC44_48h$SD_2Gy <- apply(log2(raw_HCC44_48h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pHCC44_48h<- ttest_protti(mean1 = HCC44_48h$mean_2Gy,
                                    mean2 = HCC44_48h$mean_c,
                                    sd1 = HCC44_48h$SD_2Gy,
                                    sd2 = HCC44_48h$SD_c,
                                    n1 = 4,
                                    n2 = 4,
                                    log_values=TRUE)

fdrs_pHCC44_48h<-data.frame(p.adjust(res.ttest_pHCC44_48h$pval, method="BH"))

ready_pHCC44_48h<-cbind(HCC44_48h[,1:3],res.ttest_pHCC44_48h,fdrs_pHCC44_48h,proteome_data[,c(2,95:96)])
colnames(ready_pHCC44_48h)[8] <- 'adj_pVal'

#####Proteome - Binding - H23 0h#####
raw_H23_0h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[7],proteome_data[19],proteome_data[47],proteome_data[72])
raw_H23_0h_2Gy$mean<-rowMeans(raw_H23_0h_2Gy[,4:7])

raw_H23_0h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[8],proteome_data[20],proteome_data[48],proteome_data[73])
raw_H23_0h_c$mean<-rowMeans(raw_H23_0h_c[,4:7])

H23_0h <- cbind(raw_H23_0h_c,raw_H23_0h_2Gy[,4:7])
#####Proteome - calculating Diff and pvalue - HCC23 0h####
H23_0h <- raw_H23_0h_2Gy[,1:3]
H23_0h$mean_c <- rowMeans(log2(raw_H23_0h_c[4:7]), na.rm=TRUE)
H23_0h$SD_c <- apply(log2(raw_H23_0h_c[4:7]), 1, sd, na.rm=TRUE)
H23_0h$mean_2Gy <- rowMeans(log2(raw_H23_0h_2Gy[4:7]), na.rm=TRUE)
H23_0h$SD_2Gy <- apply(log2(raw_H23_0h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pH23_0h<- ttest_protti(mean1 = H23_0h$mean_2Gy,
                                 mean2 = H23_0h$mean_c,
                                 sd1 = H23_0h$SD_2Gy,
                                 sd2 = H23_0h$SD_c,
                                 n1 = 4,
                                 n2 = 4,
                                 log_values=TRUE)

fdrs_pH23_0h<-data.frame(p.adjust(res.ttest_pH23_0h$pval, method="BH"))

ready_pH23_0h<-cbind(H23_0h[,1:3],res.ttest_pH23_0h,fdrs_pH23_0h,proteome_data[,c(2,95:96)])
colnames(ready_pH23_0h)[8] <- 'adj_pVal'

#####Proteome - Binding - H23 2h#####
raw_H23_2h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[9],proteome_data[23],proteome_data[49],proteome_data[74])
raw_H23_2h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[10],proteome_data[24],proteome_data[50],proteome_data[75])
#proteome_data[54])

H23_2h <- cbind(raw_H23_2h_c,raw_H23_2h_2Gy[,4:6])
#####Proteome - Calculating Diff and pvalue - HCC23 2h#####
H23_2h <- raw_H23_2h_2Gy[,1:3]
H23_2h$mean_c <- rowMeans(log2(raw_H23_2h_c[4:7]), na.rm=TRUE)
H23_2h$SD_c <- apply(log2(raw_H23_2h_c[4:7]), 1, sd, na.rm=TRUE)
H23_2h$mean_2Gy <- rowMeans(log2(raw_H23_2h_2Gy[4:7]), na.rm=TRUE)
H23_2h$SD_2Gy <- apply(log2(raw_H23_2h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pH23_2h<- ttest_protti(mean1 = H23_2h$mean_2Gy,
                                 mean2 = H23_2h$mean_c,
                                 sd1 = H23_2h$SD_2Gy,
                                 sd2 = H23_2h$SD_c,
                                 n1 = 4,
                                 n2 = 4,
                                 log_values=TRUE)

fdrs_pH23_2h<-data.frame(p.adjust(res.ttest_pH23_2h$pval, method="BH"))

ready_pH23_2h<-cbind(H23_2h[,1:3],res.ttest_pH23_2h,fdrs_pH23_2h,proteome_data[,c(2,95:96)])
colnames(ready_pH23_2h)[8] <- 'adj_pVal'

#####Proteome - Binding - H23 6h#####
raw_H23_6h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[11],proteome_data[25],proteome_data[84],proteome_data[76])
raw_H23_6h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[12],proteome_data[26],proteome_data[90],proteome_data[85])

H23_6h <- cbind(raw_H23_6h_c,raw_H23_6h_2Gy[,4:7])
#####Proteome - Calculating Diff and pvalue - HCC23 6h#####
H23_6h <- raw_H23_6h_2Gy[,1:3]
H23_6h$mean_c <- rowMeans(log2(raw_H23_6h_c[4:7]), na.rm=TRUE)
H23_6h$SD_c <- apply(log2(raw_H23_6h_c[4:7]), 1, sd, na.rm=TRUE)
H23_6h$mean_2Gy <- rowMeans(log2(raw_H23_6h_2Gy[4:7]), na.rm=TRUE)
H23_6h$SD_2Gy <- apply(log2(raw_H23_6h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pH23_6h<- ttest_protti(mean1 = H23_6h$mean_2Gy,
                                 mean2 = H23_6h$mean_c,
                                 sd1 = H23_6h$SD_2Gy,
                                 sd2 = H23_6h$SD_c,
                                 n1 = 4,
                                 n2 = 4,
                                 log_values=TRUE)

fdrs_pH23_6h<-data.frame(p.adjust(res.ttest_pH23_6h$pval, method="BH"))

ready_pH23_6h<-cbind(H23_6h[,1:3],res.ttest_pH23_6h,fdrs_pH23_6h,proteome_data[,c(2,95:96)])
colnames(ready_pH23_6h)[8] <- 'adj_pVal'

#####Proteome - Binding - H23 24h#####
raw_H23_24h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[13],proteome_data[27],proteome_data[55],proteome_data[86])
raw_H23_24h_2Gy$mean<-rowMeans(raw_H23_24h_2Gy[,4:7])

raw_H23_24h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[14],proteome_data[28],proteome_data[66],proteome_data[87])
raw_H23_24h_c$mean<-rowMeans(raw_H23_24h_c[,4:7])

H23_24h <- cbind(raw_H23_24h_c[,1:7],raw_H23_24h_2Gy[,4:7])
#####Proteome - Calculating Diff and pvalue - HCC23 24h####
H23_24h <- raw_H23_24h_2Gy[,1:3]
H23_24h$mean_c <- rowMeans(log2(raw_H23_24h_c[4:7]), na.rm=TRUE)
H23_24h$SD_c <- apply(log2(raw_H23_24h_c[4:7]), 1, sd, na.rm=TRUE)
H23_24h$mean_2Gy <- rowMeans(log2(raw_H23_24h_2Gy[4:7]), na.rm=TRUE)
H23_24h$SD_2Gy <- apply(log2(raw_H23_24h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pH23_24h<- ttest_protti(mean1 = H23_24h$mean_2Gy,
                                  mean2 = H23_24h$mean_c,
                                  sd1 = H23_24h$SD_2Gy,
                                  sd2 = H23_24h$SD_c,
                                  n1 = 4,
                                  n2 = 4,
                                  log_values=TRUE)

fdrs_pH23_24h<-data.frame(p.adjust(res.ttest_pH23_24h$pval, method="BH"))

ready_pH23_24h<-cbind(H23_24h[,1:3],res.ttest_pH23_24h,fdrs_pH23_24h,proteome_data[,c(2,95:96)])
colnames(ready_pH23_24h)[8] <- 'adj_pVal'

#####Proteome - Binding - H23 48h#####
raw_H23_48h_2Gy <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[15],proteome_data[29],proteome_data[77],proteome_data[88])
raw_H23_48h_c <-  cbind(proteome_data[92],proteome_data[91],proteome_data[94],proteome_data[16],proteome_data[30],proteome_data[78],proteome_data[89])

H23_48h <- cbind(raw_H23_48h_c,raw_H23_48h_2Gy[,4:7])
#####Proteome - Calculating Diff and pvalue - H23 48h####
H23_48h <- raw_H23_48h_2Gy[,1:3]
H23_48h$mean_c <- rowMeans(log2(raw_H23_48h_c[4:7]), na.rm=TRUE)
H23_48h$SD_c <- apply(log2(raw_H23_48h_c[4:7]), 1, sd, na.rm=TRUE)
H23_48h$mean_2Gy <- rowMeans(log2(raw_H23_48h_2Gy[4:7]), na.rm=TRUE)
H23_48h$SD_2Gy <- apply(log2(raw_H23_48h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pH23_48h<- ttest_protti(mean1 = H23_48h$mean_2Gy,
                                  mean2 = H23_48h$mean_c,
                                  sd1 = H23_48h$SD_2Gy,
                                  sd2 = H23_48h$SD_c,
                                  n1 = 4,
                                  n2 = 4,
                                  log_values=TRUE)

fdrs_pH23_48h<-data.frame(p.adjust(res.ttest_pH23_48h$pval, method="BH"))

ready_pH23_48h<-cbind(H23_48h[,1:3],res.ttest_pH23_48h,fdrs_pH23_48h,proteome_data[,c(2,95:96)])
colnames(ready_pH23_24h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Radiation 2h####
rad_raw_2h <- raw_HCC44_2h_2Gy[,2:4]
rad_raw_2h$mean_HCC44 <- rowMeans(log2(raw_HCC44_2h_2Gy[4:7]), na.rm=TRUE)
rad_raw_2h$SD_HCC44 <- apply(log2(raw_HCC44_2h_2Gy[4:7]), 1, sd, na.rm=TRUE)
rad_raw_2h$mean_H23 <- rowMeans(log2(raw_H23_2h_2Gy[4:7]), na.rm=TRUE)
rad_raw_2h$SD_H23 <- apply(log2(raw_H23_2h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pRad_2h<- ttest_protti(mean1 = rad_raw_2h$mean_HCC44,
                                 mean2 = rad_raw_2h$mean_H23,
                                 sd1 = rad_raw_2h$SD_HCC44 ,
                                 sd2 = rad_raw_2h$SD_H23,
                                 n1 = 4,
                                 n2 = 4,
                                 log_values=TRUE)

fdrs_pRad_2h<-data.frame(p.adjust(res.ttest_pRad_2h$pval, method="BH"))

ready_pRAD_2h<-cbind(raw_HCC44_2h_2Gy[,1:3],res.ttest_pRad_2h,fdrs_pRad_2h,proteome_data[,c(2,95:96)])
colnames(ready_pRAD_2h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Radiation 6h####
rad_raw_6h <- raw_HCC44_6h_2Gy[,2:4]
rad_raw_6h$mean_HCC44 <- rowMeans(log2(raw_HCC44_6h_2Gy[4:7]), na.rm=TRUE)
rad_raw_6h$SD_HCC44 <- apply(log2(raw_HCC44_6h_2Gy[4:7]), 1, sd, na.rm=TRUE)
rad_raw_6h$mean_H23 <- rowMeans(log2(raw_H23_6h_2Gy[4:7]), na.rm=TRUE)
rad_raw_6h$SD_H23 <- apply(log2(raw_H23_6h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pRad_6h<- ttest_protti(mean1 = rad_raw_6h$mean_HCC44,
                                 mean2 = rad_raw_6h$mean_H23,
                                 sd1 = rad_raw_6h$SD_HCC44 ,
                                 sd2 = rad_raw_6h$SD_H23,
                                 n1 = 4,
                                 n2 = 4,
                                 log_values=TRUE)

fdrs_pRad_6h<-data.frame(p.adjust(res.ttest_pRad_6h$pval, method="BH"))

ready_pRAD_6h<-cbind(raw_HCC44_6h_2Gy[,1:3],res.ttest_pRad_6h,fdrs_pRad_6h,proteome_data[,c(2,95:96)])
colnames(ready_pRAD_6h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Radiation 24h####
rad_raw_24h <- raw_HCC44_24h_2Gy[,2:4]
rad_raw_24h$mean_HCC44 <- rowMeans(log2(raw_HCC44_24h_2Gy[4:7]), na.rm=TRUE)
rad_raw_24h$SD_HCC44 <- apply(log2(raw_HCC44_24h_2Gy[4:7]), 1, sd, na.rm=TRUE)
rad_raw_24h$mean_H23 <- rowMeans(log2(raw_H23_24h_2Gy[4:7]), na.rm=TRUE)
rad_raw_24h$SD_H23 <- apply(log2(raw_H23_24h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pRad_24h<- ttest_protti(mean1 = rad_raw_24h$mean_HCC44,
                                  mean2 = rad_raw_24h$mean_H23,
                                  sd1 = rad_raw_24h$SD_HCC44 ,
                                  sd2 = rad_raw_24h$SD_H23,
                                  n1 = 4,
                                  n2 = 4,
                                  log_values=TRUE)

fdrs_pRad_24h<-data.frame(p.adjust(res.ttest_pRad_24h$pval, method="BH"))

ready_pRAD_24h<-cbind(raw_HCC44_24h_2Gy[,1:3],res.ttest_pRad_24h,fdrs_pRad_24h,proteome_data[,c(2,95:96)])
colnames(ready_pRAD_24h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Radiation 48h####
rad_raw_48h <- raw_HCC44_48h_2Gy[,2:4]
rad_raw_48h$mean_HCC44 <- rowMeans(log2(raw_HCC44_48h_2Gy[4:7]), na.rm=TRUE)
rad_raw_48h$SD_HCC44 <- apply(log2(raw_HCC44_48h_2Gy[4:7]), 1, sd, na.rm=TRUE)
rad_raw_48h$mean_H23 <- rowMeans(log2(raw_H23_48h_2Gy[4:7]), na.rm=TRUE)
rad_raw_48h$SD_H23 <- apply(log2(raw_H23_48h_2Gy[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pRad_48h<- ttest_protti(mean1 = rad_raw_48h$mean_HCC44,
                                  mean2 = rad_raw_48h$mean_H23,
                                  sd1 = rad_raw_48h$SD_HCC44 ,
                                  sd2 = rad_raw_48h$SD_H23,
                                  n1 = 4,
                                  n2 = 4,
                                  log_values=TRUE)

fdrs_pRad_48h<-data.frame(p.adjust(res.ttest_pRad_48h$pval, method="BH"))

ready_pRAD_48h<-cbind(raw_HCC44_48h_2Gy[,1:3],res.ttest_pRad_48h,fdrs_pRad_48h,proteome_data[,c(2,95:96)])
colnames(ready_pRAD_48h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Control 0h####
base_raw_0h <- raw_HCC44_0h_2Gy[,2:4]
mean_Rad <- cbind(log2(raw_HCC44_0h_c[4:7]),log2(raw_HCC44_0h_2Gy[4:7]))
base_raw_0h$mean_Rad <- rowMeans(mean_Rad, na.rm=TRUE)
base_raw_0h$SD_Rad <- apply(mean_Rad[4:7], 1, sd, na.rm=TRUE)

mean_H23 <- cbind(log2(raw_H23_0h_c[4:7]),log2(raw_H23_0h_2Gy[4:7]))
base_raw_0h$mean_H23 <- rowMeans(mean_H23, na.rm=TRUE)
base_raw_0h$SD_H23 <- apply(mean_H23, 1, sd, na.rm=TRUE)

res.ttest_pbase_0h<- ttest_protti(mean1 = base_raw_0h$mean_Rad,
                                  mean2 = base_raw_0h$mean_H23,
                                  sd1 = base_raw_0h$SD_Rad ,
                                  sd2 = base_raw_0h$SD_H23,
                                  n1 = 8,
                                  n2 = 8,
                                  log_values=TRUE)

fdrs_pbase_0h<-data.frame(p.adjust(res.ttest_pbase_0h$pval, method="BH"))

ready_pbase_0h<-cbind(raw_HCC44_0h_c[,1:3],res.ttest_pbase_0h,fdrs_pbase_0h,proteome_data[,c(2,95:96)])
colnames(ready_pbase_0h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Control 2h####
Con_raw_2h <- raw_HCC44_2h_2Gy[,2:4]
Con_raw_2h$mean_HCC44 <- rowMeans(log2(raw_HCC44_2h_c[4:7]), na.rm=TRUE)
Con_raw_2h$SD_HCC44 <- apply(log2(raw_HCC44_2h_c[4:7]), 1, sd, na.rm=TRUE)
Con_raw_2h$mean_H23 <- rowMeans(log2(raw_H23_2h_c[4:7]), na.rm=TRUE)
Con_raw_2h$SD_H23 <- apply(log2(raw_H23_2h_c[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pCon_2h<- ttest_protti(mean1 = Con_raw_2h$mean_HCC44,
                                 mean2 = Con_raw_2h$mean_H23,
                                 sd1 = Con_raw_2h$SD_HCC44 ,
                                 sd2 = Con_raw_2h$SD_H23,
                                 n1 = 4,
                                 n2 = 4,
                                 log_values=TRUE)

fdrs_pCon_2h<-data.frame(p.adjust(res.ttest_pCon_2h$pval, method="BH"))

ready_pCon_2h<-cbind(raw_HCC44_2h_c[,1:3],res.ttest_pCon_2h,fdrs_pCon_2h,proteome_data[,c(2,95:96)])
colnames(ready_pCon_2h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Control 6h####
Con_raw_6h <- raw_HCC44_6h_c[,2:4]
Con_raw_6h$mean_HCC44 <- rowMeans(log2(raw_HCC44_6h_c[4:7]), na.rm=TRUE)
Con_raw_6h$SD_HCC44 <- apply(log2(raw_HCC44_6h_c[4:7]), 1, sd, na.rm=TRUE)
Con_raw_6h$mean_H23 <- rowMeans(log2(raw_H23_6h_c[4:7]), na.rm=TRUE)
Con_raw_6h$SD_H23 <- apply(log2(raw_H23_6h_c[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pCon_6h<- ttest_protti(mean1 = Con_raw_6h$mean_HCC44,
                                 mean2 = Con_raw_6h$mean_H23,
                                 sd1 = Con_raw_6h$SD_HCC44 ,
                                 sd2 = Con_raw_6h$SD_H23,
                                 n1 = 4,
                                 n2 = 4,
                                 log_values=TRUE)

fdrs_pCon_6h<-data.frame(p.adjust(res.ttest_pCon_6h$pval, method="BH"))

ready_pCon_6h<-cbind(raw_HCC44_6h_c[,1:3],res.ttest_pCon_6h,fdrs_pCon_6h,proteome_data[,c(2,95:96)])
colnames(ready_pCon_6h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Control 24h####
Con_raw_24h <- raw_HCC44_24h_c[,2:4]
Con_raw_24h$mean_HCC44 <- rowMeans(log2(raw_HCC44_24h_c[4:7]), na.rm=TRUE)
Con_raw_24h$SD_HCC44 <- apply(log2(raw_HCC44_24h_c[4:7]), 1, sd, na.rm=TRUE)
Con_raw_24h$mean_H23 <- rowMeans(log2(raw_H23_24h_c[4:7]), na.rm=TRUE)
Con_raw_24h$SD_H23 <- apply(log2(raw_H23_24h_c[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pCon_24h<- ttest_protti(mean1 = Con_raw_24h$mean_HCC44,
                                  mean2 = Con_raw_24h$mean_H23,
                                  sd1 = Con_raw_24h$SD_HCC44 ,
                                  sd2 = Con_raw_24h$SD_H23,
                                  n1 = 4,
                                  n2 = 4,
                                  log_values=TRUE)

fdrs_pCon_24h<-data.frame(p.adjust(res.ttest_pCon_24h$pval, method="BH"))

ready_pCon_24h<-cbind(raw_HCC44_24h_c[,1:3],res.ttest_pCon_24h,fdrs_pCon_24h,proteome_data[,c(2,95:96)])
colnames(ready_pCon_24h)[8] <- 'adj_pVal'
#####Proteome - Calculating Diff and pvalue - Control 48h####
Con_raw_48h <- raw_HCC44_48h_c[,2:4]
Con_raw_48h$mean_HCC44 <- rowMeans(log2(raw_HCC44_48h_c[4:7]), na.rm=TRUE)
Con_raw_48h$SD_HCC44 <- apply(log2(raw_HCC44_48h_c[4:7]), 1, sd, na.rm=TRUE)
Con_raw_48h$mean_H23 <- rowMeans(log2(raw_H23_48h_c[4:7]), na.rm=TRUE)
Con_raw_48h$SD_H23 <- apply(log2(raw_H23_48h_c[4:7]), 1, sd, na.rm=TRUE)

res.ttest_pCon_48h<- ttest_protti(mean1 = Con_raw_48h$mean_HCC44,
                                  mean2 = Con_raw_48h$mean_H23,
                                  sd1 = Con_raw_48h$SD_HCC44 ,
                                  sd2 = Con_raw_48h$SD_H23,
                                  n1 = 4,
                                  n2 = 4,
                                  log_values=TRUE)

fdrs_pCon_48h<-data.frame(p.adjust(res.ttest_pCon_48h$pval, method="BH"))

ready_pCon_48h<-cbind(raw_HCC44_48h_c[,1:3],res.ttest_pCon_48h,fdrs_pCon_48h,proteome_data[,c(2,95:96)])
colnames(ready_pCon_48h)[8] <- 'adj_pVal'
#####Saving results#####
setwd(save_wd)
write.csv(ready_pHCC44_0h,"stat_pHCC44_0h.csv")
write.csv(ready_pHCC44_2h,"stat_pHCC44_2h.csv")
write.csv(ready_pHCC44_6h,"stat_pHCC44_6h.csv")
write.csv(ready_pHCC44_24h,"stat_pHCC44_24h.csv")
write.csv(ready_pHCC44_48h,"stat_pHCC44_48h.csv")

write.csv(raw_HCC44_0h_c, "raw_pHCC44_0h_c.csv")
write.csv(raw_HCC44_0h_2Gy, "raw_pHCC44_0h_2Gy.csv")
write.csv(raw_HCC44_2h_c, "raw_pHCC44_2h_c.csv")
write.csv(raw_HCC44_2h_2Gy, "raw_pHCC44_2h_2Gy.csv")
write.csv(raw_HCC44_6h_c, "raw_pHCC44_6h_c.csv")
write.csv(raw_HCC44_6h_2Gy, "raw_pHCC44_6h_2Gy.csv")
write.csv(raw_HCC44_24h_c, "raw_pHCC44_24h_c.csv")
write.csv(raw_HCC44_24h_2Gy, "raw_pHCC44_24h_2Gy.csv")
write.csv(raw_HCC44_48h_c, "raw_pHCC44_48h_c.csv")
write.csv(raw_HCC44_48h_2Gy, "raw_pHCC44_48h_2Gy.csv")

write.csv(HCC44_0h, "mean_pHCC44_0h.csv")
write.csv(HCC44_2h, "mean_pHCC44_2h.csv")
write.csv(HCC44_6h, "mean_pHCC44_6h.csv")
write.csv(HCC44_24h, "mean_pHCC44_24h.csv")
write.csv(HCC44_48h, "mean_pHCC44_48h.csv")

write.csv(ready_pH23_0h,"stat_pH23_0h.csv")
write.csv(ready_pH23_2h,"stat_pH23_2h.csv")
write.csv(ready_pH23_6h,"stat_pH23_6h.csv")
write.csv(ready_pH23_24h,"stat_pH23_24h.csv")
write.csv(ready_pH23_48h,"stat_pH23_48h.csv")

write.csv(raw_H23_0h_c, "raw_pH23_0h_c.csv")
write.csv(raw_H23_0h_2Gy, "raw_pH23_0h_2Gy.csv")
write.csv(raw_H23_2h_c, "raw_pH23_2h_c.csv")
write.csv(raw_H23_2h_2Gy, "raw_pH23_2h_2Gy.csv")
write.csv(raw_H23_6h_c, "raw_pH23_6h_c.csv")
write.csv(raw_H23_6h_2Gy, "raw_pH23_6h_2Gy.csv")
write.csv(raw_H23_24h_c, "raw_pH23_24h_c.csv")
write.csv(raw_H23_24h_2Gy, "raw_pH23_24h_2Gy.csv")
write.csv(raw_H23_48h_c, "raw_pH23_48h_c.csv")
write.csv(raw_H23_48h_2Gy, "raw_pH23_48h_2Gy.csv")

write.csv(H23_0h, "mean_pH23_0h.csv")
write.csv(H23_2h, "mean_pH23_2h.csv")
write.csv(H23_6h, "mean_pH23_6h.csv")
write.csv(H23_24h, "mean_pH23_24h.csv")
write.csv(H23_48h, "mean_pH23_48h.csv")

#write.csv(ready_pRAD_0h,"stat_pRAD_0h.csv")
write.csv(ready_pRAD_2h,"stat_pRAD_2h.csv")
write.csv(ready_pRAD_6h,"stat_pRAD_6h.csv")
write.csv(ready_pRAD_24h,"stat_pRAD_24h.csv")
write.csv(ready_pRAD_48h,"stat_pRAD_48h.csv")

write.csv(ready_pbase_0h,"stat_pCon_0h.csv")
write.csv(ready_pCon_2h,"stat_pCon_2h.csv")
write.csv(ready_pCon_6h,"stat_pCon_6h.csv")
write.csv(ready_pCon_24h,"stat_pCon_24h.csv")
write.csv(ready_pCon_48h,"stat_pCon_48h.csv")
setwd(source_wd)