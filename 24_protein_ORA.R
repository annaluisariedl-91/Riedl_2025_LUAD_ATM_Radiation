setwd(save_wd)
#####HCC44 - binding Data####
HCC44_timeData<- cbind(ready_pHCC44_0h[,c(1,2,11,4,7)],ready_pHCC44_2h[,c(4,7)],ready_pHCC44_6h[,c(4,7)],ready_pHCC44_24h[,c(4,7)],ready_pHCC44_48h[,c(4,7)])
HCC44_timeData<-HCC44_timeData[complete.cases(HCC44_timeData[,3:4]),]
HCC44_timeData_all<- cbind(raw_HCC44_0h_c[,1:3],log2(raw_HCC44_0h_c[,4:7]),log2(raw_HCC44_0h_2Gy[,4:7]),log2(raw_HCC44_2h_c[,4:7]),log2(raw_HCC44_2h_2Gy[,4:7]),log2(raw_HCC44_6h_c[,4:7]),log2(raw_HCC44_6h_2Gy[,4:7]),log2(raw_HCC44_24h_c[,4:7]),log2(raw_HCC44_24h_2Gy[,4:7]),log2(raw_HCC44_48h_c[,4:7]),log2(raw_HCC44_48h_2Gy[,4:7]))
#####H23 - binding Data####
H23_timeData<- cbind(ready_pH23_0h[,c(1,2,11,4,7)],ready_pH23_2h[,c(4,7)],ready_pH23_6h[,c(4,7)],ready_pH23_24h[,c(4,7)],ready_pH23_48h[,c(4,7)])
H23_timeData_all<- cbind(raw_H23_0h_c[,1:3],log2(raw_H23_0h_c[,4:7]),log2(raw_H23_0h_2Gy[,4:7]),log2(raw_H23_2h_c[,4:7]),log2(raw_H23_2h_2Gy[,4:7]),log2(raw_H23_6h_c[,4:7]),log2(raw_H23_6h_2Gy[,4:7]),log2(raw_H23_24h_c[,4:7]),log2(raw_H23_24h_2Gy[,4:7]),log2(raw_H23_48h_c[,4:7]),log2(raw_H23_48h_2Gy[,4:7]))
#####Proteom - binding mean#####
all_mean <-cbind(H23_0h[,c(1:4,6)],H23_2h[,c(4,6)],H23_6h[,c(4,6)],H23_24h[,c(4,6)],H23_48h[,c(4,6)],
                 HCC44_0h[,c(4,6)],HCC44_2h[,c(4,6)],HCC44_6h[,c(4,6)],HCC44_24h[,c(4,6)],HCC44_48h[,c(4,6)])
Data_all <-merge(x=H23_timeData_all, y=HCC44_timeData_all[,c(2,4:43)], by="ID")
#####HCC44 - GO ORA#####
GO_HCC44 <- merge(x=HCC44_Cluster,y=HCC44_timeData[,c(1:2)], by="ID")

entrez_ids_pHCC44_Cluster <- bitr(GO_HCC44[,2], 
                                  fromType = "SYMBOL", 
                                  toType = "ENTREZID", 
                                  OrgDb = org.Hs.eg.db)
res_pHCC44_Cluster<-enrichGO(gene = entrez_ids_pHCC44_Cluster[,2], OrgDb = org.Hs.eg.db,ont = "BP",    
                             pAdjustMethod = "BH",
                             readable = TRUE)
go_pHCC44_Cluster<- as.data.frame(res_pHCC44_Cluster)
go_pHCC44_Cluster<-go_pHCC44_Cluster[!duplicated(go_pHCC44_Cluster$geneID),]
go_pHCC44_Cluster<-go_pHCC44_Cluster[order(go_pHCC44_Cluster$Count, decreasing = TRUE),]
go_pHCC44_Cluster<-go_pHCC44_Cluster[1:10,]
go_pHCC44_Cluster[["p.adjust_log"]]<- -log10(go_pHCC44_Cluster$p.adjust)
go_pHCC44_Cluster[["gene_ratio_numeric"]] <- sapply(go_pHCC44_Cluster$GeneRatio, function(s) eval(parse(text = s)))
jpeg("GO_ORA_HCC44.jpeg",  width=41, height=21, unit="cm",quality=100,res=70)
go_pHCC44_Cluster %>%
  ggplot(aes(x=gene_ratio_numeric, y=10:1,  group =  p.adjust_log, size = Count,fill =p.adjust_log )) +
  geom_point(alpha=0.75, shape=21,stroke=1.5)+
  theme_bw()+
  guides(fill = guide_legend(override.aes = list(size=6)))+
  scale_fill_distiller(palette = "Greens", breaks = c(3,5,7,9), limits=c(2,10),direction=FALSE,name = bquote(-log[10]~"pValue"))+
  scale_size(range = c(5, 15), name="# Proteins")+
  xlab("Gene ratio") +
  ylab("") +
  labs(title="HCC44")+
  scale_x_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, by = 0.05))+
  scale_y_continuous(labels=go_pHCC44_Cluster$Description,breaks = seq(1, 10, by = 1))+
  theme(axis.text.x = element_text(color = "black", size = 28, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 24, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 34, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 24, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.text = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        title= element_text(color = "black", size = 36, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(color = "black", size = 1.25),
        axis.ticks.length = unit(0.25, "cm")
  )
dev.off()
#####H23 - GO ORA####
GO_H23 <- merge(x=H23_Cluster,y=H23_timeData[,1:2], by="ID")

entrez_ids_pH23_cluster_24h <- bitr(GO_H23[,2], 
                                    fromType = "SYMBOL", 
                                    toType = "ENTREZID", 
                                    OrgDb = org.Hs.eg.db)
res_pH23_cluster<-enrichGO(gene = entrez_ids_pH23_cluster_24h[,2], OrgDb = org.Hs.eg.db,ont = "BP",    
                           pAdjustMethod = "BH",
                           readable = TRUE)
go_pH23_cluster<- as.data.frame(res_pH23_cluster)
go_pH23_cluster<-go_pH23_cluster[1:8,]
go_pH23_cluster[["p.adjust_log"]]<- -log10(go_pH23_cluster$p.adjust)
go_pH23_cluster[["gene_ratio_numeric"]] <- sapply(go_pH23_cluster$GeneRatio, function(s) eval(parse(text = s)))
jpeg("GO_ORA_H23.jpeg",  width=41, height=21, unit="cm",quality=100,res=70)
go_pH23_cluster %>%
  ggplot(aes(x=gene_ratio_numeric, y=8:1,  group =  p.adjust_log, size = Count,fill =p.adjust_log )) +
  geom_point(alpha=0.75, shape=21,stroke=1.5) +
  theme_bw()+
  guides(fill = guide_legend(override.aes = list(size=6)))+
  scale_fill_distiller(palette = "Greens", breaks = c(1, 2,4,6), limits=c(1,7),direction=FALSE,name = bquote(-log[10]~"pValue"))+
  scale_size(range = c(5, 15), name="# Proteins")+
  xlab("Gene ratio") +
  ylab("") +
  labs(title="NCI-H23")+
  scale_x_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, by = 0.05))+
  scale_y_continuous(labels=c(go_pH23_cluster[8:1,2]),breaks = seq(1, 8, by = 1))+
  theme(axis.text.x = element_text(color = "black", size = 28, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 24, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 34, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 24, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.text = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        title= element_text(color = "black", size = 36, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(color = "black", size = 1.25),
        axis.ticks.length = unit(0.25, "cm")
  )
dev.off()
#####Plottinge cell cylce cluster#####
HCC44_Cluster3<-data.frame(ID=c("2812","4781","5336","1540","2078","2397","2739","2780","2914","3100","3285","334","3867","4223","4722",
                                "5146","5287","5322","6341","7031","7170","7385","7459","8050","8871","8941","9036",
                                "9162","944","9683","9750"))

HCC44_Cluster3_sym<-data.frame(Suggested.Symbol=c("ANLN","AURKA","BUB1","CCNB1","CCNB2","CDC20","CENPE","CKAP2","GTSE1","HMMR","PLK1","TROAP"))
HCC44_Cluster3<-data.frame(ID=c("8871" , "2739", "2914" ,"3867" ,"3285","5287", "5146", "7385" , "9162" ,"3100" ,"4722", "2780", "2078","5370","7530"))
HCC44_Cluster3_all_mean<-merge(x=HCC44_Cluster3,y=all_mean, by="ID",all.x=TRUE) 
HCC44_Cluster3_all_mean <- HCC44_Cluster3_all_mean[!grepl("8330", HCC44_Cluster3_all_mean$ID),] 
HCC44_Cluster3_all_mean<- HCC44_Cluster3_all_mean[!duplicated(HCC44_Cluster3_all_mean[,1:2]),]
HCC44_Cluster3_all_mean<- HCC44_Cluster3_all_mean[order(HCC44_Cluster3_all_mean[,2]),]
#####HCC44 and H23  - time course cluster 3####
jpeg(filename = "Cluster3_ATM.jpeg",  width=41, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1:2),ncol=2,nrow =1,byrow = TRUE))
leg.txt<- c("HCC44 2Gy","NCI-H23 2Gy")
par(mar = c(10, 13, 5, 6),xpd=TRUE)
for(i in c(2)){
  plot(c(0,2,6,24,48),(as.numeric(unlist(HCC44_Cluster3_all_mean[i,c(9:13)]))/as.numeric(unlist(HCC44_Cluster3_all_mean[i,c(9)]))),ann=FALSE,axes=FALSE,col= "#F46D43",type="l",lwd=3,ylim=c(0.95,1.15))
  lines(c(0,2,6,24,48),(as.numeric(unlist(HCC44_Cluster3_all_mean[i,c(19:23)]))/as.numeric(unlist(HCC44_Cluster3_all_mean[i,c(19)]))), lwd=3,col= "#3288BD")
  box(which = "plot",lwd=3)
  axis(1,lwd=3,at=c(0,24,48),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,mgp=c(12,3,0))
  axis(2,lwd=3,at=c(1,1.1),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,las=2,mgp=c(12,3,0))
  lines(c(0,48),c(1,1),lty=3,col="grey60",lwd=2)
  title(ylab= bquote(atop(log[2]~italic("(Intensity)"), bold("normalized"))), cex.lab=3,cex.main=2,line=7,adj=c(0.5))
  title(xlab="Time [h]", cex.lab=3,cex.main=3,line=6,adj=c(0.5))
  title( main=HCC44_Cluster3_all_mean[i,2] , cex.lab=3,cex.main=3,line=1.5,adj=c(0.5))
  #minor.tick(nx =1,ny=2/1,tick.ratio =1, x.args = list( lwd = 4,tcl = -1), y.args = list( lwd = 4,tcl = -1)) 
  #legend(x=c(0,48),y=c(1.15,1.05),ncol=1 ,bty = "n", bg = NA  ,leg.txt,pch=21,cex =2.25,pt.cex = 2.5,pt.bg = c("#3288BD","#F46D43"),  col="black",box.lwd =3,pt.lwd = 3,y.intersp=1,x.intersp=0.75)
}
par(mar = c(10, 2, 5, 6),xpd=TRUE)
plot(1:10,1:10,ann=FALSE, axes=FALSE, col="white")
legend(x=c(0,4),y=c(3,10),ncol=1 ,bty = "n", bg = NA  ,leg.txt,pch=21,cex =2.25,pt.cex = 2.5,pt.bg = c("#3288BD","#F46D43"), col="black",box.lwd =3,pt.lwd = 3,y.intersp=1,x.intersp=0.75)
dev.off()

exprss_H23<-NULL
for(i in c(1,3:nrow(HCC44_Cluster3_all_mean))){
  cur <- data.frame(HCC44_Cluster3_all_mean[i,9:13])
  cur <- cur/cur[,1]
  cur<- data.frame(t(cur))
  colnames(cur) <- "mean"
  cur$Time <- c(0,2,6,24,48)
  cur$Protein <- HCC44_Cluster3_all_mean[i,2]
  exprss_H23<- rbind(exprss_H23, cur)
}
exprss_H23$cellLine <- "NCI-H23"

exprss_HCC44<-NULL
for(i in c(1,3:nrow(HCC44_Cluster3_all_mean))){
  cur <- data.frame(HCC44_Cluster3_all_mean[i,19:23])
  cur <- cur/cur[,1]
  cur<- data.frame(t(cur))
  colnames(cur) <- "mean"
  cur$Time <- c(0,2,6,24,48)
  cur$Protein <- HCC44_Cluster3_all_mean[i,2]
  exprss_HCC44<- rbind(exprss_HCC44, cur)
}
exprss_HCC44$cellLine <- "HCC44"

exprss <- rbind(exprss_H23,exprss_HCC44)

jpeg(filename = "Cluster3_all_ggplot.jpeg",  width=41, height=28, unit="cm",quality=100,res=70)
ggplot(exprss, aes(x = Time, y = mean, color = cellLine, group = cellLine)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Protein) +
  theme_bw()+
  scale_color_manual(values = c("NCI-H23" = "#3288BD","HCC44" =  "#F46D43"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size=0.75)+
  labs(
    title = "Expression per Protein",
    x = "Time point",
    y = "Expression",
    color = "Cell line")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 28, angle = 0, hjust = .5, vjust = 0, face = "plain",margin = margin(r = 20)),
        axis.title.y = element_text(color = "black", size = 28, angle = 90, hjust = .5, vjust = .5, face = "plain",margin = margin(r = 20)),
        legend.title = element_text(color = "black", size = 24, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.text = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        title= element_text(color = "black", size = 36, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(color = "black", size = 1.25),
        axis.ticks.length = unit(0.25, "cm"),
        strip.text = element_text(size = 20, face = "bold"),
        panel.spacing = unit(1.5, "lines")  
  )
dev.off()
#####

exprss_H23<-NULL
for(i in c(2)){
  cur <- data.frame(HCC44_Cluster3_all_mean[i,9:13])
  #  cur <- cur/cur[,1]
  cur<- data.frame(t(cur))
  colnames(cur) <- "mean"
  cur$Time <- c(0,2,6,24,48)
  cur$Protein <- HCC44_Cluster3_all_mean[i,2]
  exprss_H23<- rbind(exprss_H23, cur)
}
exprss_H23$cellLine <- "NCI-H23"

exprss_HCC44<-NULL
for(i in c(2)){
  cur <- data.frame(HCC44_Cluster3_all_mean[i,19:23])
  #  cur <- cur/cur[,1]
  cur<- data.frame(t(cur))
  colnames(cur) <- "mean"
  cur$Time <- c(0,2,6,24,48)
  cur$Protein <- HCC44_Cluster3_all_mean[i,2]
  exprss_HCC44<- rbind(exprss_HCC44, cur)
}
exprss_HCC44$cellLine <- "HCC44"

exprss <- rbind(exprss_H23,exprss_HCC44)

jpeg(filename = "Cluster3_ATM_ggplot.jpeg",  width=15.15, height=7, unit="cm",quality=100,res=70)
ggplot(exprss, aes(x = Time, y = mean, color = cellLine, group = cellLine)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_bw()+
  ylim(13, 16) +
  scale_color_manual(values = c("NCI-H23" = "#3288BD","HCC44" =  "#F46D43"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size=0.75)+
  labs(
    title = "ATM",
    x = "Time point",
    y = "Expression",
    color = "Cell line")+
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain",margin = margin(r = 20)),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain",margin = margin(r = 20)),
        legend.title = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.text = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        title= element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(color = "black", size = 1.25),
        axis.ticks.length = unit(0.25, "cm"),
        strip.text = element_text(size = 20, face = "bold")
  )
dev.off()
#####
exprss_H23<-NULL
for(i in c(3,4,12)){
  cur <- data.frame(HCC44_Cluster3_all_mean[i,9:13])
  cur <- cur/cur[,1]
  cur<- data.frame(t(cur))
  colnames(cur) <- "mean"
  cur$Time <- c(0,2,6,24,48)
  cur$Protein <- HCC44_Cluster3_all_mean[i,2]
  exprss_H23<- rbind(exprss_H23, cur)
}
exprss_H23$cellLine <- "NCI-H23"

exprss_HCC44<-NULL
for(i in c(3,4,12)){
  cur <- data.frame(HCC44_Cluster3_all_mean[i,19:23])
  cur <- cur/cur[,1]
  cur<- data.frame(t(cur))
  colnames(cur) <- "mean"
  cur$Time <- c(0,2,6,24,48)
  cur$Protein <- HCC44_Cluster3_all_mean[i,2]
  exprss_HCC44<- rbind(exprss_HCC44, cur)
}
exprss_HCC44$cellLine <- "HCC44"

exprss <- rbind(exprss_H23,exprss_HCC44)

jpeg(filename = "Cluster3_Trio_ggplot.jpeg",  width=30.75, height=10, unit="cm",quality=100,res=70)
ggplot(exprss, aes(x = Time, y = mean, color = cellLine, group = cellLine)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Protein) +
  theme_bw()+
  scale_color_manual(values = c("NCI-H23" = "#3288BD","HCC44" =  "#F46D43"))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size=0.75)+
  labs(
    title = "Expression per protein",
    x = "Time point",
    y = "Expression",
    color = "Cell line")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 28, angle = 0, hjust = .5, vjust = 0, face = "plain",margin = margin(r = 20)),
        axis.title.y = element_text(color = "black", size = 28, angle = 90, hjust = .5, vjust = .5, face = "plain",margin = margin(r = 20)),
        legend.title = element_text(color = "black", size = 28, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.text = element_text(color = "black", size = 24, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        title= element_text(color = "black", size = 36, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.ticks = element_line(color = "black", size = 1.25),
        axis.ticks.length = unit(0.25, "cm"),
        strip.text = element_text(size = 20, face = "bold"),
        panel.spacing = unit(1.5, "lines")  
  )
dev.off()
setwd(source_wd)