####Transcriptome - HCC44 24h - binding Data####
RHCC44_24h_2Gy <-cbind(Rdata[2],Rdata[24],Rdata[13:16],Rdata[,c(1,8)])
RHCC44_24h_c <- cbind(Rdata[2],Rdata[24],Rdata[9:12],Rdata[,c(1,8)])
#RHCC44_24h_c[,3:6][RHCC44_24h_c[,3:6]==0]<-NA

RHCC44_24h <- data.frame(cbind(RHCC44_24h_c[,1:2],log2(RHCC44_24h_c[,3:6]),log2(RHCC44_24h_2Gy[,3:6])),RHCC44_24h_c[,7:8])
colnames(RHCC44_24h)[11]<-"ENSEMBL"
colnames(RHCC44_24h)[12]<-"gene_biotype"
####Transcriptome - HCC44 24h - diff and Pval####
RHCC44_24h[["mean_c"]] <- rowMeans(RHCC44_24h[,3:6], na.rm=TRUE)
RHCC44_24h[["SD_c"]] <- apply(RHCC44_24h[,c(3:6)], 1, sd, na.rm=TRUE)
RHCC44_24h[["mean_2Gy"]] <- rowMeans(RHCC44_24h[,7:10], na.rm=TRUE)
RHCC44_24h[["SD_2Gy"]] <- apply(RHCC44_24h[,c(7:10)], 1, sd, na.rm=TRUE)

res.ttest_RHCC44_24h<- ttest_protti(mean1 = RHCC44_24h$mean_2Gy,
                                    mean2 = RHCC44_24h$mean_c,
                                    sd1 = RHCC44_24h$SD_2Gy,
                                    sd2 = RHCC44_24h$SD_c,
                                    n1 = 4,
                                    n2 = 4,
                                    log_values=TRUE)

fdrs_RHCC44_24h<-data.frame(p.adjust(res.ttest_RHCC44_24h$pval, method="BH"))

ready_RHCC44_24h<-cbind(RHCC44_24h[,1:2],res.ttest_RHCC44_24h,fdrs_RHCC44_24h,RHCC44_24h[,11:12])
colnames(ready_RHCC44_24h)[7] <- 'adj_pVal'

####Transcriptome - H23 24h - binding Data####
RH23_24h_2Gy <-cbind(Rdata[2],Rdata[24],Rdata[,c(20:23)],Rdata[,c(1,8)])
RH23_24h_c <- cbind(Rdata[2],Rdata[24],Rdata[,c(17:19)],Rdata[,c(1,8)])
RH23_24h <- data.frame(cbind(RH23_24h_c[,1:2],log2(RH23_24h_c[,3:5]),log2(RH23_24h_2Gy[,3:6]),RH23_24h_c[,6:7]))
colnames(RH23_24h)[10]<-"ENSEMBL"
colnames(RH23_24h)[11]<-"gene_biotype"
####Transcriptome - H23 24h - diff and Pval####
RH23_24h[["mean_c"]] <- rowMeans(RH23_24h[,3:5],na.rm = TRUE)
RH23_24h[["SD_c"]] <- apply(RH23_24h[,c(3:5)], 1, sd,na.rm = TRUE)
RH23_24h[["mean_2Gy"]] <- rowMeans(RH23_24h[,6:9],na.rm = TRUE)
RH23_24h[["SD_2Gy"]] <- apply(RH23_24h[,c(6:9)], 1, sd,na.rm = TRUE)

res.ttest_RH23_24h<- ttest_protti(mean1 = RH23_24h$mean_2Gy,
                                  mean2 = RH23_24h$mean_c,
                                  sd1 = RH23_24h$SD_2Gy,
                                  sd2 = RH23_24h$SD_c,
                                  n1 = 4,
                                  n2 = 3,
                                  log_values=TRUE)

fdrs_RH23_24h<-data.frame(p.adjust(res.ttest_RH23_24h$pval, method="BH"))

ready_RH23_24h<-cbind(RH23_24h[,1:2],res.ttest_RH23_24h,fdrs_RH23_24h,RH23_24h[,10:11])

colnames(ready_RH23_24h)[7] <- 'adj_pVal'
####Transcriptome - Radiation 24h - diff and Pval####
RNA_RAD_24h <- cbind(RH23_24h[,c(1,2,14,15)],RHCC44_24h[,c(15:16)])
colnames(RNA_RAD_24h)<- c("gene_name","ID","mean_H23","SD_H23","mean_HCC44","SD_HCC44")

res.ttest_RNA_RAD_24h<- ttest_protti(mean1 = RNA_RAD_24h$mean_HCC44,
                                     mean2 = RNA_RAD_24h$mean_H23,
                                     sd2 = RNA_RAD_24h$SD_H23,
                                     sd1 = RNA_RAD_24h$SD_HCC44,
                                     n1 = 4,
                                     n2 = 4,
                                     log_values=TRUE)

fdrs_RNA_RAD_24h<-data.frame(p.adjust(res.ttest_RNA_RAD_24h$pval, method="BH"))

ready_RRAD_24h<-cbind(RH23_24h[,c(1:2)],res.ttest_RNA_RAD_24h,fdrs_RNA_RAD_24h,RH23_24h[,10:11])
colnames(ready_RRAD_24h)[7] <- 'adj_pVal'
####Transcriptome - Control 24h - diff and Pval####
RNA_Con_24h <- cbind(RH23_24h[,c(1,2,12,13)],RHCC44_24h[,c(13:14)])
colnames(RNA_Con_24h)<- c("gene_name","ID","mean_H23","SD_H23","mean_HCC44","SD_HCC44")

res.ttest_RNA_Con_24h<- ttest_protti(mean1 = RNA_Con_24h$mean_HCC44,
                                     mean2 = RNA_Con_24h$mean_H23,
                                     sd2 = RNA_Con_24h$SD_H23,
                                     sd1 = RNA_Con_24h$SD_HCC44,
                                     n1 = 4,
                                     n2 = 4,
                                     log_values=TRUE)

fdrs_RNA_Con_24h<-data.frame(p.adjust(res.ttest_RNA_Con_24h$pval, method="BH"))

ready_RCon_24h<-cbind(RH23_24h[,c(1:2)],res.ttest_RNA_Con_24h,fdrs_RNA_Con_24h,RH23_24h[,10:11])
colnames(ready_RCon_24h)[7] <- 'adj_pVal'
####Transcriptome - Saving results####
setwd(save_wd)
ready_RHCC44_24h<- cbind(ready_RHCC44_24h[,c(1,2)], rowMeans(cbind(RHCC44_24h[,3:6],RHCC44_24h[,7:10]),na.rm=TRUE),ready_RHCC44_24h[,c(3:7)],RHCC44_24h[,c(13,15)],ready_RHCC44_24h[,c(8,9)])
colnames(ready_RHCC44_24h)<- c("GeneSymbol","GeneID","baseMean","log2FoldChange","std_error","t_statistic","pval","adj_pVal", "mean_control","mean_2Gy","ENSEMBL","gene_biotype")
write.table(ready_RHCC44_24h, "DE_HCC44_2Gy_vs_0Gy.csv",
            sep = ",", row.names = FALSE, quote = FALSE)

ready_RH23_24h<- cbind(ready_RH23_24h[,c(1,2)], rowMeans(cbind(RH23_24h[,3:5],RH23_24h[,6:9]),na.rm=TRUE),ready_RH23_24h[,c(3:7)],RH23_24h[,c(12,14)],ready_RH23_24h[,c(8,9)])
colnames(ready_RH23_24h)<- c("GeneSymbol","GeneID","baseMean","log2FoldChange","std_error","t_statistic","pval","adj_pVal", "mean_control","mean_2Gy","ENSEMBL","gene_biotype")
write.table(ready_RH23_24h, "DE_H23_2Gy_vs_0Gy.csv",
            sep = ",", row.names = FALSE, quote = FALSE)

ready_RRAD_24h<- cbind(ready_RRAD_24h[,c(1,2)], rowMeans(cbind(RH23_24h[,6:9],RHCC44_24h[,7:10]),na.rm=TRUE),ready_RRAD_24h[,c(3:7)],RHCC44_24h[,c(15)],RH23_24h[,c(14)],ready_RRAD_24h[,c(8,9)])
colnames(ready_RRAD_24h)<- c("GeneSymbol","GeneID","baseMean","log2FoldChange","std_error","t_statistic","pval","adj_pVal", "mean_HCC44","mean_H23","ENSEMBL","gene_biotype")
write.table(ready_RRAD_24h, "DE_HCC44_vs_H23_2Gy.csv",
            sep = ",", row.names = FALSE, quote = FALSE)

ready_RCon_24h<- cbind(ready_RCon_24h[,c(1,2)], rowMeans(cbind(RH23_24h[,3:5],RHCC44_24h[,3:6]),na.rm=TRUE),ready_RCon_24h[,c(3:7)],RHCC44_24h[,c(13)],RH23_24h[,c(12)],ready_RCon_24h[,c(8,9)])
colnames(ready_RCon_24h)<- c("GeneSymbol","GeneID","baseMean","log2FoldChange","std_error","t_statistic","pval","adj_pVal", "mean_HCC44","mean_H23","ENSEMBL","gene_biotype")
write.table(ready_RCon_24h, "DE_HCC44_vs_H23_0Gy.csv",
            sep = ",", row.names = FALSE, quote = FALSE)
setwd(source_wd)