#####HCC44 - binding Data####
setwd(save_wd)
HCC44_timeData_all<- cbind(raw_HCC44_0h_c[,1:3],log2(raw_HCC44_0h_c[,4:7]),log2(raw_HCC44_0h_2Gy[,4:7]),log2(raw_HCC44_2h_c[,4:7]),log2(raw_HCC44_2h_2Gy[,4:7]),log2(raw_HCC44_6h_c[,4:7]),log2(raw_HCC44_6h_2Gy[,4:7]),log2(raw_HCC44_24h_c[,4:7]),log2(raw_HCC44_24h_2Gy[,4:7]),log2(raw_HCC44_48h_c[,4:7]),log2(raw_HCC44_48h_2Gy[,4:7]))
#####HCC44 - Spline regression####
exprs_HCC44_time<- HCC44_timeData_all[,c(1:2,8:11,16:19, 24:27, 32:35,40:43)]
colnames(exprs_HCC44_time)<- c("Sugested.Symbol","ID",paste0("Sample", 1:16))
targets_HCC44_time <- data.frame(
  Sample = paste0("Sample", 1:20),
  Time = factor(rep(c(rep("h0",4), rep("h2",4), rep("h6",4), rep("h24",4), rep("h48",4)),1)), levels = c(rep("h0",4), rep("h2",4), rep("h6",4), rep("h24",4), rep("h48",4)))
design <- model.matrix(~ns(c(rep(0,4),rep(2,4),rep(6,4),rep(24,4),rep(48,4)), df = 3), data = targets_HCC44_time)
colnames(design)<- c("Intercept","Spline1","Spline2","Spline3")
fit_HCC44_time <- lmFit(exprs_HCC44_time[,3:22], design)
contrast.matrix <- makeContrasts(
  Spline1,
  Spline2,
  Spline3,
  levels = design
)

fit_HCC44_time <- contrasts.fit(fit_HCC44_time, contrast.matrix)
fit2_HCC44_time <- eBayes(fit_HCC44_time)
results_HCC44_time <- topTable(fit2_HCC44_time, coef =1:3, adjust = "BH", number = Inf)
res_sig_HCC44_time <- results_HCC44_time[which(results_HCC44_time$adj.P.Val<0.05),]
#####HCC44 - protein content selected with splines####
HCC44_Cluster <- data.frame(ID=c("1631","6775","7801","808","8329","8389","8927","1108","1142","1218","1801","2670","3400",
                                 "3576","3710","3910","4364","4639","4786","4862","5078","5123","5148","5283","5760","5813",
                                 "5922","6211","662","6822","6955","7082","7649","7726","7943","8581","9202","9345" ,                               
                                 "1","1736","2360","328","3473","3498","3511","3542","3699","3702","3885","3938","40",  
                                 "4250","4944","5180","5447","5808","5931","6149","7328","7699","8027","8364","8880","9876",                        
                                 "2812","4781","5336","1540","2078","2397","2739","2780","2914","3100","3285","334","3867",
                                 "4223","4722","5146","5287","5322","6341","7031","7170","7385","7459","8050","8871","8941",
                                 "9036","9162","944","9683","9750"))
res_sig_HCC44_all<- merge(x=HCC44_Cluster, y=HCC44_timeData_all, by="ID")
colnames(res_sig_HCC44_all)<- c("ID","suggested.Symbol","ENSEMBL",names_prot[1:40])
res_sig_HCC44_all<-res_sig_HCC44_all[complete.cases(res_sig_HCC44_all[4:43]),]
ready_Spline_time_HCC44_all_hm1<- t(apply(res_sig_HCC44_all[,4:43], 1, cal_z_score))
rownames(ready_Spline_time_HCC44_all_hm1)<- res_sig_HCC44_all[,2]
Graph_Spline_time_HCC44_all_hm1<-pheatmap(ready_Spline_time_HCC44_all_hm1, method="complete", names.arg=res_sig_HCC44_all[,2], label_row=res_sig_HCC44_all[,2],annotation_col = annotation_prot,  fontsize_row = 3,  fontsize_col = 15, show_rownames=TRUE, show_colnames=FALSE)
save_pheatmap_png(Graph_Spline_time_HCC44_all_hm1,"Graph_Spline_time_HCC44_all_hm1.png")

res_sig_HCC44_all<- merge(x=HCC44_Cluster, y=HCC44_timeData_all, by="ID")
colnames(res_sig_HCC44_all)<- c("ID","suggested.Symbol","ENSEMBL",names_prot[1:40])
res_sig_HCC44_all<-res_sig_HCC44_all[complete.cases(res_sig_HCC44_all[c(8:11,16:19,24:27,32:35,40:43)]),c(1:2,8:11,16:19,24:27,32:35,40:43)]
ready_Spline_time_HCC44_Rad_hm1<- t(apply(res_sig_HCC44_all[,3:22], 1, cal_z_score))
rownames(ready_Spline_time_HCC44_Rad_hm1)<- res_sig_HCC44_all[,2]
Graph_Spline_time_HCC44_Rad_hm1<-pheatmap(ready_Spline_time_HCC44_Rad_hm1,method="complete",col=my_colors, names.arg=res_sig_HCC44_all[,3], cluster_col=FALSE, annotation_col = annotation_prot,  fontsize_row =7,  fontsize_col = 15, show_rownames=TRUE, show_colnames=FALSE)
save_pheatmap_png(Graph_Spline_time_HCC44_Rad_hm1,"Graph_Spline_time_HCC44_Rad_hm1.png")
#####H23 - binding Data####
H23_timeData_all<- cbind(raw_H23_0h_c[,1:3],log2(raw_H23_0h_c[,4:7]),log2(raw_H23_0h_2Gy[,4:7]),log2(raw_H23_2h_c[,4:7]),log2(raw_H23_2h_2Gy[,4:7]),log2(raw_H23_6h_c[,4:7]),log2(raw_H23_6h_2Gy[,4:7]),log2(raw_H23_24h_c[,4:7]),log2(raw_H23_24h_2Gy[,4:7]),log2(raw_H23_48h_c[,4:7]),log2(raw_H23_48h_2Gy[,4:7]))
H23_timeData_all[,4:43] <- lapply(H23_timeData_all[,4:43], function(x) {
  x[x == "NaN"] <- NA           
  as.numeric(x)                
})
#####H23 - Spline regression####
exprs_H23_time <- H23_timeData_all[,c(1:2,8:11,16:19, 24:27, 32:35,40:43)]
colnames(exprs_H23_time)<- c("Sugested.Symbol","ID",paste0("Sample", 1:20))
targets_H23_time <- data.frame(
  Sample = paste0("Sample", 1:20),
  Time = factor(rep(c(rep("h0",4), rep("h2",4), rep("h6",4), rep("h24",4), rep("h48",4)),1)), levels = c(rep("h0",4), rep("h2",4), rep("h6",4), rep("h24",4), rep("h48",4)))
design <- model.matrix(~ ns(c(rep(0,4),rep(2,4),rep(6,4),rep(24,4),rep(48,4)), df = 3), data = targets_H23_time)
colnames(design)<- c("Intersept","Spline1" ,"Spline2","Spline3")
fit_H23_time <- lmFit(exprs_H23_time[,3:22], design)
contrast.matrix <- makeContrasts(
  Spline1,
  Spline2,
  Spline3,
  levels = design
)

fit_H23_time <- contrasts.fit(fit_H23_time, contrast.matrix)
fit2_H23_time <- eBayes(fit_H23_time)
results_H23_time <- topTable(fit2_H23_time, coef =1:3, adjust = "BH", number = Inf)
res_sig_H23_time <- results_H23_time[which(results_H23_time$P.Val<0.05),]
#####H23 - protein content selected with splines####
H23_Cluster <- data.frame(ID=c("1020","1125","1182","12","1662","1714","1905","2722","2853","3002","3030","3047","3190",
                              "3259","3395","3446","3461","3599","3657","3834","4033","4187","4202","4299","455","4638",
                               "4661","4719","5283","5412", "5509" ,"5660" ,"5697", "5713", "6051" ,"6531", "6532", "6549", "6621",
                              "7043", "717","7801", "7918", "806","81","8327","8389","9018","9313","9560","9606","9854",
                               "9889","1469","1530","1838","2397","2685","2731","3163","3417","3445","3673","3702","3827","4143",
                               "4259","4418","445","4654","4688","4769","5087","5099","5382","5401","5493","5534","5725",
                               "5780","5930","5931","5969","6306","6341","6748","7081","7103","7596","7718","7815","7998",
                               "8005","8007","8038","8046","8079","8181","8322","8347","8387","8404","8477","8488","8599",
                                "8703","8799","9007","9038","9185","944","9477","9701","9824"))
res_sig_H23_all<- merge(x=H23_Cluster, y=H23_timeData_all[,c(1:43)], by="ID")
res_sig_H23_all<-res_sig_H23_all[complete.cases(res_sig_H23_all[,4:43]),]
colnames(res_sig_H23_all)<- c("ID","suggested.Symbol","ENSEMBL",names_prot[41:80])
ready_Spline_time_H23_all_hm1<- t(apply(res_sig_H23_all[,4:43], 1, cal_z_score))
rownames(ready_Spline_time_H23_all_hm1)<- res_sig_H23_all[,2]
Graph_Spline_time_H23_all_hm1<-pheatmap(ready_Spline_time_H23_all_hm1, method="complete", names.arg=res_sig_H23_all[,2], label_row=res_sig_H23_all[,2],annotation_col = annotation_prot[41:80,],  fontsize_row = 3,  fontsize_col = 15, show_rownames=TRUE, show_colnames=FALSE)
save_pheatmap_png(Graph_Spline_time_H23_all_hm1,"Graph_Spline_time_H23_all_hm1.png")

res_sig_H23_all<- merge(x=H23_Cluster, y=H23_timeData_all[,c(1:43)], by="ID")
res_sig_H23_all<-res_sig_H23_all[complete.cases(res_sig_H23_all[,c(8:11,16:19,24:27,32:35,40:43)]),]
colnames(res_sig_H23_all)<- c("ID","suggested.Symbol","ENSEMBL",names_prot[41:80])
ready_Spline_time_H23_all_hm1<- t(apply(res_sig_H23_all[,c(8:11,16:19,24:27,32:35,40:43)], 1, cal_z_score))
rownames(ready_Spline_time_H23_all_hm1)<- res_sig_H23_all[,2]
Graph_Spline_time_H23_Rad_hm1<-pheatmap(ready_Spline_time_H23_Rad_hm1,method="complete",col=my_colors, names.arg=res_sig_H23_all[,3], cluster_col=FALSE, annotation_col = annotation_prot,  fontsize_row =7,  fontsize_col = 15, show_rownames=TRUE, show_colnames=FALSE)
save_pheatmap_png(Graph_Spline_time_H23_Rad_hm1,"Graph_Spline_time_H23_Rad_hm1.png")
dev.off()
setwd(source_wd)
