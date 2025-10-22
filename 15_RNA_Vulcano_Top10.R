setwd(save_wd)
#####Transcriptome - Volcano Plot - Control#####
ready_RCon_24h[ready_RCon_24h == "NaN"] <- NA
down_RCon_24h<-ready_RCon_24h[which(ready_RCon_24h[["pval"]] < 0.05&ready_RCon_24h[["log2FoldChange"]]<= 	-4),]
down_RCon_24h<-down_RCon_24h[order((down_RCon_24h[,4]*-1),decreasing = TRUE),]
down_RCon_24h$Nr <- 1:nrow(down_RCon_24h)

up_RCon_24h<-ready_RCon_24h[which(ready_RCon_24h[["pval"]] < 0.05&ready_RCon_24h[["log2FoldChange"]]  >= 	4),]
up_RCon_24h<-up_RCon_24h[order(up_RCon_24h[,4],decreasing = TRUE),]
up_RCon_24h$Nr <- 1:nrow(up_RCon_24h)
#####Transcriptome - Plotting volcano - Control####
jpeg(filename = "Transcriptome_Control_volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(8, 10, 6, 4),xpd=TRUE)
plot((ready_RCon_24h[,4]),-log10(ready_RCon_24h[,7]),pch=16,ylim=c(0,15), xlim = c(-20,20), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="Transcriptome - Control 24h",mgp=c(7,2,0))
points(up_RCon_24h[,4],-log10(up_RCon_24h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_RCon_24h [,4],-log10(down_RCon_24h [,7]), col="#D53E4F",pch=16, cex=2)
lines(c(-4,-4),c(0,15.5), lty=3,lwd=2)
lines(c(4,4),c(0,15.5), lty=3,lwd=2)
lines(c(-20,20),c(1.3,1.3),lwd=2,lty=3)
lines(c(0,0),c(0,16),lty=2,col="grey60",lwd=2)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
axis(1,lwd=3,at = c(-20,-10,0,10,20),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10,15),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=3)
title(ylab=bquote(log[10]~"(pValue)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(HCC44/H23)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_RCon_24h[c(41),4])+8,(-log10(up_RCon_24h[c(41),7]))+1.5,labels=up_RCon_24h[c(41),1],cex=2.5)
text((up_RCon_24h[c(3),4])+1.5,(-log10(up_RCon_24h[c(3),7]))+1.0,labels=up_RCon_24h[c(3),1],cex=2.5)
text((down_RCon_24h[c(1),4])-2.5,(-log10(down_RCon_24h[c(1),7]))+1.0,labels=down_RCon_24h[c(1),1],cex=2.5)
text((down_RCon_24h[c(260),4])-6.5,(-log10(down_RCon_24h[c(260),7]))+1,labels=down_RCon_24h[c(260),1],cex=2.5)
dev.off()
#####Transcriptome - Top 10 - Control####
sig.pVal_RCon_24h_down<-ready_RCon_24h[which(ready_RCon_24h[["pval"]] < 0.05&ready_RCon_24h[["log2FoldChange"]] <(-1.25)),]
sig.pVal_RCon_24h_down<-sig.pVal_RCon_24h_down[order(sig.pVal_RCon_24h_down[,4], decreasing=FALSE),]
sig.pVal_RCon_24h_down<- sig.pVal_RCon_24h_down[1:10,]

sig.pVal_RCon_24h_up<-ready_RCon_24h[which(ready_RCon_24h[["pval"]] < 0.05&ready_RCon_24h[["log2FoldChange"]] >(1.25)),]
sig.pVal_RCon_24h_up<-sig.pVal_RCon_24h_up[order(sig.pVal_RCon_24h_up[,4], decreasing=TRUE),]
sig.pVal_RCon_24h_up <- sig.pVal_RCon_24h_up[1:10,]
sig.pVal_RCon_24h_up<- sig.pVal_RCon_24h_up[order(sig.pVal_RCon_24h_up[,4],decreasing=TRUE),]
#####
sig.con_pbase_0h_up<-ready_pbase_0h[which(ready_pbase_0h[["pval"]] < 0.05),]
sig.con_pbase_0h_up<-sig.con_pbase_0h_up[order(sig.con_pbase_0h_up[,5], decreasing=TRUE),]
sig.con_pbase_0h_up$rank <- 1:nrow(sig.con_pbase_0h_up)

sig.pVal_RCon_24h_up<-merge(x= sig.pVal_RCon_24h_up,y=sig.con_pbase_0h_up[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RCon_24h_up<-sig.pVal_RCon_24h_up[order(sig.pVal_RCon_24h_up[,5], decreasing=TRUE),]
#####
sig.con_pbase_0h_down<-ready_pbase_0h[which(ready_pbase_0h[["pval"]] < 0.05),]
sig.con_pbase_0h_down<-sig.con_pbase_0h_down[order(sig.con_pbase_0h_down[,5], decreasing=TRUE),]
sig.con_pbase_0h_down$rank <- 1:nrow(sig.con_pbase_0h_down)

sig.pVal_RCon_24h_down<-merge(x= sig.pVal_RCon_24h_down,y=sig.con_pbase_0h_down[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RCon_24h_down<-sig.pVal_RCon_24h_down[order(sig.pVal_RCon_24h_down[,5], decreasing=FALSE),]
#####
sig.pVal_RCon_24h_all <- rbind(sig.pVal_RCon_24h_up[10:1,c(2,5,8,13:16)],sig.pVal_RCon_24h_down[10:1,c(2,5,8,13:16)])
pval_RCon_24h<- data.frame(-log10(sig.pVal_RCon_24h_all[,3]))
pval_RCon_24h$order = findInterval(pval_RCon_24h[,1], sort(pval_RCon_24h[,1]))
v_pval_RCon_24h<- c(0,(max(pval_RCon_24h[,1])*(2/10)),(max(pval_RCon_24h[,1])*(4/10)),(max(pval_RCon_24h[,1])*(6/10)),(max(pval_RCon_24h[,1])*(8/10)))
pval_RCon_24h$order = findInterval(pval_RCon_24h[,1], v_pval_RCon_24h)
lgd_pval_RCon_24h = rep(NA, 5)
lgd_pval_RCon_24h[c(1,5)] = c(0,substring(max(pval_RCon_24h[,1]),  1,4))

pval_RCon_24h_2<- data.frame(-log10(sig.pVal_RCon_24h_all[,6]))
pval_RCon_24h_2$order = findInterval(pval_RCon_24h_2[,1], sort(pval_RCon_24h_2[,1]))
v_pval_RCon_24h_2<- c(0,(max(pval_RCon_24h_2[,1],na.rm=TRUE)*(2/10)),(max(pval_RCon_24h_2[,1],na.rm=TRUE)*(4/10)),(max(pval_RCon_24h_2[,1],na.rm=TRUE)*(6/10)),(max(pval_RCon_24h_2[,1],na.rm=TRUE)*(8/10)))
v_pval_RCon_24h_2<- v_pval_RCon_24h_2[complete.cases(v_pval_RCon_24h_2)]
pval_RCon_24h_2$order = findInterval(pval_RCon_24h_2[,2], v_pval_RCon_24h_2)
lgd_pval_RCon_24h_2 = rep(NA, 5)
lgd_pval_RCon_24h_2[c(1,5)] = c(0,substring(max(pval_RCon_24h_2[,1],na.rm=TRUE),  1,3))
#####Transcriptome - Plotting Top10 Control####
jpeg(filename = "Transcriptome_Control_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(8, 10, 6, 4),xpd=TRUE)
sig.pVal_RCon_24h_all<-sig.pVal_RCon_24h_all[order(sig.pVal_RCon_24h_all[,2], decreasing=TRUE),]
plot(sig.pVal_RCon_24h_all[c(20:11,10:1),2],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_Con(5)[pval_RCon_24h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=4,line=2, cex.main=3,main="Transcriptome - Control 24h",outer=FALSE,adj=c(0.5))
title(xlab=bquote(log[2]~"(HCC44/NCI-H23)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
sig.pVal_RCon_24h_all<-sig.pVal_RCon_24h_all[order(sig.pVal_RCon_24h_all[,2], decreasing=FALSE),]
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.25),labels = sig.pVal_RCon_24h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_RCon_24h_all[i,2],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  lines(c(sig.pVal_RCon_24h_all[i,2],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  text(c(16.5),c(i,i), srt=0, cex=2.5,labels=ifelse(!is.na(sig.pVal_RCon_24h_all[i,7]), sig.pVal_RCon_24h_all[i,7],"*"))
}
points(sig.pVal_RCon_24h_all[20:1,5],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_RCon_24h_2[20:1,2]],cex=3,lwd=2)
sig.pVal_RCon_24h_all<-sig.pVal_RCon_24h_all[order(sig.pVal_RCon_24h_all[,2], decreasing=TRUE),]
points(sig.pVal_RCon_24h_all[20:1,5],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_RCon_24h_2[20:1,2]],cex=3,lwd=2)
for (i in 11:20) {
  text(x = 1,y = c(21-i),adj=c(0,0.4),labels = sig.pVal_RCon_24h_all[i,1],srt = 0,cex = 2.5)
}
legend(x = -28, y =6, title = "\n Transcriptome \n -log10(p-Value)\n",
       legend = lgd_pval_RCon_24h,
       fill = colorRampPalette(colors =cols_RAD)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
legend(x = -28, y =12, title = "\n Proteom\n -log10(p-Value)\n",
       legend = lgd_pval_RCon_24h_2,
       fill = colorRampPalette(colors =cols_Con)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
text(c(-27),c(0), srt=0, cex=2,labels="* not significant",adj=c(0,0))
text(c(15.5),c(20.5), srt=0, cex=2.5,labels="Position",adj=c(0.5,0))
#####
dev.off()
#####Transcriptome - Volcano - Radiation#####
ready_RRAD_24h[ready_RRAD_24h == "NaN"] <- NA
down_RRAD_24h<-ready_RRAD_24h[which(ready_RRAD_24h[["pval"]] < 0.05&ready_RRAD_24h[["log2FoldChange"]]<= 	-4),]
down_RRAD_24h<-down_RRAD_24h[order((down_RRAD_24h[,4]*-1),decreasing = TRUE),]
down_RRAD_24h$Nr <- 1:nrow(down_RRAD_24h)

up_RRAD_24h<-ready_RRAD_24h[which(ready_RRAD_24h[["pval"]] < 0.05&ready_RRAD_24h[["log2FoldChange"]]  >= 	4),]
up_RRAD_24h<-up_RRAD_24h[order(up_RRAD_24h[,4],decreasing = TRUE),]
up_RRAD_24h$Nr <- 1:nrow(up_RRAD_24h)
#####Transcriptome - Plotting volcano - Radiation#####
jpeg(filename = "Transcriptome_Radiation_volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2, 1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(10, 10, 6, 4),xpd=TRUE)
plot((ready_RRAD_24h[,4]),-log10(ready_RRAD_24h[,7]),pch=16,ylim=c(0,15), xlim = c(-20,20), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="",mgp=c(7,2,0))
points(up_RRAD_24h[,4],-log10(up_RRAD_24h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_RRAD_24h [,4],-log10(down_RRAD_24h [,7]), col="#D53E4F",pch=16, cex=2)
lines(c(-4,-4),c(0,15.5), lty=3,lwd=2)
lines(c(4,4),c(0,15.5), lty=3,lwd=2)
lines(c(-21,21),c(1.3,1.3),lwd=2, lty=3)
lines(c(0,0),c(0,15),lty=2,col="grey",lwd=2)
minor.tick(nx =6/3,ny=3,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
axis(1,lwd=3,at = c(-20,-10,0,10,20),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10,15),lwd.ticks = 3,cex=3, cex.axis=2.5,tcl = -1.5,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=3)
title(main="Transcriptome - Radiation 24h",cex.lab=3,cex.main=3,line=2,adj=c(1))
title(ylab=bquote(-log[10]~"(p-Value)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(HCC44/NCI-H23)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_RRAD_24h[c(1),4])+2.5,(-log10(up_RRAD_24h[c(1),7]))-1,labels=up_RRAD_24h[c(1),1],cex=2.5)
text((up_RRAD_24h[c(6),4])+2.5,(-log10(up_RRAD_24h[c(6),7]))+1,labels=up_RRAD_24h[c(6),1],cex=2.5)
text((down_RRAD_24h[c(13),4])-5,(-log10(down_RRAD_24h[c(13),7]))+1,labels=down_RRAD_24h[c(13),1],cex=2.5)
text((down_RRAD_24h[c(75),4])-7.75,(-log10(down_RRAD_24h[c(75),7]))+1,labels=down_RRAD_24h[c(75),1],cex=2.5)
dev.off()
#####Transcriptome - Top 10 - Radiation####
sig.pVal_RRAD_24h_down<-ready_RRAD_24h[which(ready_RRAD_24h[["pval"]] < 0.05&ready_RRAD_24h[["log2FoldChange"]] <(-1.25)),]
sig.pVal_RRAD_24h_down<-sig.pVal_RRAD_24h_down[order(sig.pVal_RRAD_24h_down[,4], decreasing=FALSE),]
sig.pVal_RRAD_24h_down<- sig.pVal_RRAD_24h_down[1:10,]

sig.pVal_RRAD_24h_up<-ready_RRAD_24h[which(ready_RRAD_24h[["pval"]] < 0.05&ready_RRAD_24h[["log2FoldChange"]] >(1.25)),]
sig.pVal_RRAD_24h_up<-sig.pVal_RRAD_24h_up[order(sig.pVal_RRAD_24h_up[,4], decreasing=TRUE),]
sig.pVal_RRAD_24h_up <- sig.pVal_RRAD_24h_up[1:10,]
sig.pVal_RRAD_24h_up<- sig.pVal_RRAD_24h_up[order(sig.pVal_RRAD_24h_up[,4],decreasing=TRUE),]
#####
sig.con_pRAD_24h_up<-ready_pRAD_24h[which(ready_pRAD_24h[["adj_pVal"]] < 0.05),]
sig.con_pRAD_24h_up<-sig.con_pRAD_24h_up[order(sig.con_pRAD_24h_up[,5], decreasing=TRUE),]
sig.con_pRAD_24h_up$rank <- 1:nrow(sig.con_pRAD_24h_up)

sig.pVal_RRAD_24h_up<-merge(x= sig.pVal_RRAD_24h_up,y=sig.con_pRAD_24h_up[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RRAD_24h_up<-sig.pVal_RRAD_24h_up[order(sig.pVal_RRAD_24h_up[,5], decreasing=TRUE),]
#####
sig.con_pRAD_24h_down<-ready_pRAD_24h[which(ready_pRAD_24h[["adj_pVal"]] < 0.05),]
sig.con_pRAD_24h_down<-sig.con_pRAD_24h_down[order(sig.con_pRAD_24h_down[,5], decreasing=TRUE),]
sig.con_pRAD_24h_down$rank <- 1:nrow(sig.con_pRAD_24h_down)

sig.pVal_RRAD_24h_down<-merge(x= sig.pVal_RRAD_24h_down,y=sig.con_pRAD_24h_down[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RRAD_24h_down<-sig.pVal_RRAD_24h_down[order(sig.pVal_RRAD_24h_down[,5], decreasing=FALSE),]
#####
sig.pVal_RRAD_24h_all <- rbind(sig.pVal_RRAD_24h_up[10:1,c(2,5,8,14:16)],sig.pVal_RRAD_24h_down[10:1,c(2,5,8,14:16)])
pval_RRAD_24h<- data.frame(-log10(sig.pVal_RRAD_24h_all[,3]))
pval_RRAD_24h$order = findInterval(pval_RRAD_24h[,1], sort(pval_RRAD_24h[,1]))
v_pval_RRAD_24h<- c(0,(max(pval_RRAD_24h[,1])*(2/10)),(max(pval_RRAD_24h[,1])*(4/10)),(max(pval_RRAD_24h[,1])*(6/10)),(max(pval_RRAD_24h[,1])*(8/10)))
pval_RRAD_24h$order = findInterval(pval_RRAD_24h[,1], v_pval_RRAD_24h)
lgd_pval_RRAD_24h = rep(NA, 5)
lgd_pval_RRAD_24h[c(1,5)] = c(0,substring(max(pval_RRAD_24h[,1]),  1,4))

pval_RRAD_24h_2<- data.frame(-log10(as.numeric(sig.pVal_RRAD_24h_all[,6])))
pval_RRAD_24h_2$order = findInterval(pval_RRAD_24h_2[,1], sort(pval_RRAD_24h_2[,1]))
#v_pval_RRAD_24h_2<- c(0,(max(pval_RRAD_24h_2[,1],na.rm=TRUE)*(2/10)),(max(pval_RRAD_24h_2[,1],na.rm=TRUE)*(4/10)),(max(pval_RRAD_24h_2[,1],na.rm=TRUE)*(6/10)),(max(pval_RRAD_24h_2[,1],na.rm=TRUE)*(8/10)))
#v_pval_RRAD_24h_2<- v_pval_RRAD_24h_2[complete.cases(v_pval_RRAD_24h_2)]
#pval_RRAD_24h_2$order = findInterval(pval_RRAD_24h_2[,2], v_pval_RRAD_24h_2)
lgd_pval_RRAD_24h_2 = rep(NA, 5)
#lgd_pval_RRAD_24h_2[c(1,5)] = c(0,substring(max(pval_RRAD_24h_2[,1],na.rm=TRUE),  1,4))
#####Transcriptome - Plotting Top10 - Radiation####
jpeg(filename = "Transcriptome_Radiation_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(8, 10, 6, 4),xpd=TRUE)
sig.pVal_RRAD_24h_all<-sig.pVal_RRAD_24h_all[order(sig.pVal_RRAD_24h_all[,2], decreasing=TRUE),]
plot(sig.pVal_RRAD_24h_all[c(20:11,10:1),2],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_Con(5)[pval_RRAD_24h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=3,line=2, cex.main=3,main="Transcriptome - Radiation 24h",outer=FALSE,adj=c(0.6))
title(xlab=bquote(log[2]~"(HCC44/NCI-H23)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
sig.pVal_RRAD_24h_all<-sig.pVal_RRAD_24h_all[order(sig.pVal_RRAD_24h_all[,2], decreasing=TRUE),]
for (i in 1:10) {
  text(x = -1,y = c(10+i),adj=c(1,0.2),labels = sig.pVal_RRAD_24h_all[11-i,1],srt = 0,cex = 2.5)
}
for (i in 1:10) {
  lines(c(sig.pVal_RRAD_24h_all[11-i,2],0),c(10+i,10+i), lty=1,lwd=2)
}

sig.pVal_RRAD_24h_all<-sig.pVal_RRAD_24h_all[order(sig.pVal_RRAD_24h_all[,2], decreasing=FALSE),]
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.65),labels = sig.pVal_RRAD_24h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 1:10) {
  lines(c(sig.pVal_RRAD_24h_all[i,2],0),c(i,i), lty=1,lwd=2)
}

points(sig.pVal_RRAD_24h_all[1:20,5],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_RRAD_24h_2[20:1,2]],cex=3,lwd=2)
for (i in 1:20) {
  lines(c(sig.pVal_RRAD_24h_all[i,5],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  text(c(16.5),c(i,i), srt=0, cex=2.5,labels=ifelse(!is.na(sig.pVal_RRAD_24h_all[i,6]), sig.pVal_RRAD_24h_all[i,6],"*"))
}
legend(x = -28, y =12, title = "\n Proteome\n -log10(p-Value)\n",
       legend = lgd_pval_RRAD_24h,
       fill = colorRampPalette(colors =cols_RAD)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
#legend(x = -28, y =6, title = "\n Transcriptome \n -log10(p-Value)\n",
#       legend = lgd_pval_RRAD_24h_2,
#       fill = colorRampPalette(colors =cols_Con)(5),
#       border = "black",
#       y.intersp = 0.5,title.cex =2,
#       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
text(c(-27),c(0), srt=0, cex=2,labels="* not significant",adj=c(0,0))
text(c(16.5),c(20.5), srt=0, cex=2.5,labels="Position",adj=c(0.5,0))
#####
dev.off()
#####Transcriptome - Volcano - HCC44#####
ready_RHCC44_24h[ready_RHCC44_24h == "NaN"] <- NA
down_RHCC44_24h<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05&ready_RHCC44_24h[["log2FoldChange"]]<= 	-1),]
down_RHCC44_24h<-down_RHCC44_24h[order((down_RHCC44_24h[,4]*-1),decreasing = TRUE),]
down_RHCC44_24h$Nr <- 1:nrow(down_RHCC44_24h)

up_RHCC44_24h<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05&ready_RHCC44_24h[["log2FoldChange"]]  >= 	1),]
up_RHCC44_24h<-up_RHCC44_24h[order(up_RHCC44_24h[,4],decreasing = TRUE),]
up_RHCC44_24h$Nr <- 1:nrow(up_RHCC44_24h)
#####Transcriptome - Plotting  volcano HCC44#####
jpeg(filename = "Transcriptome_HCC44_volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2, 1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(10, 10, 6, 4),xpd=TRUE)
plot((ready_RHCC44_24h[,4]),-log10(ready_RHCC44_24h[,7]),pch=16,ylim=c(0,10), xlim = c(-10,10), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="RNA - HCC44 24h",mgp=c(7,2,0))
points(up_RHCC44_24h[,4],-log10(up_RHCC44_24h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_RHCC44_24h [,4],-log10(down_RHCC44_24h [,7]), col= "#D53E4F",pch=16, cex=2)
lines(c(-1,-1),c(0,10.5), lty=3,lwd=2)
lines(c(1,1),c(0,10.5), lty=3,lwd=2)
lines(c(-11,11),c(1.3,1.3),lwd=2, lty=3)
lines(c(0,0),c(0,10),lty=2,col="grey60",lwd=2)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 3,tcl =  -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=3)
title(ylab=bquote(-log[10]~"(pValue)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_RHCC44_24h[c(1),4])+3.25,(-log10(up_RHCC44_24h[c(1),7]))+1,labels=up_RHCC44_24h[c(1),1],cex=2.5)
text((up_RHCC44_24h[c(11),4])+2.25,(-log10(up_RHCC44_24h[c(11),7]))+1,labels=up_RHCC44_24h[c(11),1],cex=2.5)
text((down_RHCC44_24h[c(2),4])-2.5,(-log10(down_RHCC44_24h[c(2),7]))-1.0,labels=down_RHCC44_24h[c(2),1],cex=2.5)
text((down_RHCC44_24h[c(1),4])-2.5,(-log10(down_RHCC44_24h[c(1),7]))+1.0,labels=down_RHCC44_24h[c(1),1],cex=2.5)
dev.off()
#####Transcriptome - Top 10 HCC44####
sig.pVal_RHCC44_24h_down<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05&ready_RHCC44_24h[["log2FoldChange"]] <(-1)),]
sig.pVal_RHCC44_24h_down<-sig.pVal_RHCC44_24h_down[order(sig.pVal_RHCC44_24h_down[,4], decreasing=FALSE),]
sig.pVal_RHCC44_24h_down<- sig.pVal_RHCC44_24h_down[1:10,]

sig.pVal_RHCC44_24h_up<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05&ready_RHCC44_24h[["log2FoldChange"]] >(1)),]
sig.pVal_RHCC44_24h_up<-sig.pVal_RHCC44_24h_up[order(sig.pVal_RHCC44_24h_up[,5], decreasing=TRUE),]
sig.pVal_RHCC44_24h_up <- sig.pVal_RHCC44_24h_up[1:10,]
sig.pVal_RHCC44_24h_up<- sig.pVal_RHCC44_24h_up[order(sig.pVal_RHCC44_24h_up[,5],decreasing=TRUE),]
#####
sig.con_pHCC44_24h_up<-ready_pHCC44_24h[which(ready_pHCC44_24h[["pval"]] < 0.05),]
sig.con_pHCC44_24h_up<-sig.con_pHCC44_24h_up[order(sig.con_pHCC44_24h_up[,5], decreasing=TRUE),]
sig.con_pHCC44_24h_up$rank <- 1:nrow(sig.con_pHCC44_24h_up)

sig.pVal_RHCC44_24h_up<-merge(x= sig.pVal_RHCC44_24h_up,y=sig.con_pHCC44_24h_up[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RHCC44_24h_up<-sig.pVal_RHCC44_24h_up[order(sig.pVal_RHCC44_24h_up[,5], decreasing=TRUE),]
#####
sig.con_pHCC44_24h_down<-ready_pHCC44_24h[which(ready_pHCC44_24h[["pval"]] < 0.05),]
sig.con_pHCC44_24h_down<-sig.con_pHCC44_24h_down[order(sig.con_pHCC44_24h_down[,5], decreasing=TRUE),]
sig.con_pHCC44_24h_down$rank <- 1:nrow(sig.con_pHCC44_24h_down)

sig.pVal_RHCC44_24h_down<-merge(x= sig.pVal_RHCC44_24h_down,y=sig.con_pHCC44_24h_down[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RHCC44_24h_down<-sig.pVal_RHCC44_24h_down[order(sig.pVal_RHCC44_24h_down[,5], decreasing=FALSE),]
#####
sig.pVal_RHCC44_24h_all <- rbind(sig.pVal_RHCC44_24h_up[10:1,c(2,5,8,13,14,16)],sig.pVal_RHCC44_24h_down[10:1,c(2,5,8,13,14,16)])
pval_RHCC44_24h<- data.frame(-log10(sig.pVal_RHCC44_24h_all[,3]))
pval_RHCC44_24h$order = findInterval(pval_RHCC44_24h[,1], sort(pval_RHCC44_24h[,1]))
v_pval_RHCC44_24h<- c(0,(max(pval_RHCC44_24h[,1],na.rm=TRUE)*(2/10)),(max(pval_RHCC44_24h[,1],na.rm=TRUE)*(4/10)),(max(pval_RHCC44_24h[,1],na.rm=TRUE)*(6/10)),(max(pval_RHCC44_24h[,1],na.rm=TRUE)*(8/10)))
pval_RHCC44_24h$order = findInterval(pval_RHCC44_24h[,1], v_pval_RHCC44_24h)
lgd_pval_RHCC44_24h = rep(NA, 5)
lgd_pval_RHCC44_24h[c(1,5)] = c(0,substring(max(pval_RHCC44_24h[,1]),  1,4))
#####Transcriptome - Plotting Top10 - HCC44#####
jpeg(filename = "Transcriptome_HCC44_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2, 1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(8, 10, 6, 4),xpd=TRUE)
sig.pVal_RHCC44_24h_all<-sig.pVal_RHCC44_24h_all[order(sig.pVal_RHCC44_24h_all[,2], decreasing = FALSE),]
plot(sig.pVal_RHCC44_24h_all[1:20,2],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_Con(5)[pval_RHCC44_24h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=3,line=2, cex.main=3,main="Transcriptome - NCI-HCC44 24h",outer=FALSE,adj=c(0.7))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
sig.pVal_RHCC44_24h_all<-sig.pVal_RHCC44_24h_all[order(sig.pVal_RHCC44_24h_all[,2], decreasing = FALSE),]
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.85),labels = sig.pVal_RHCC44_24h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.1),labels = sig.pVal_RHCC44_24h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_RHCC44_24h_all[i,2],0),c(i,i), lty=1,lwd=2)
}
#points(sig.pVal_RHCC44_24h_all[1:20,4],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_RHCC44_24h_2$order],cex=3,lwd=2)
for (i in 1:20) {
  lines(c(sig.pVal_RHCC44_24h_all[i,4],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  text(c(15.5),c(i,i), srt=0, cex=2.5,labels=ifelse(!is.na(sig.pVal_RHCC44_24h_all[i,6]), sig.pVal_RHCC44_24h_all[i,6],"*"))
}
legend(x = -28, y =12, title = "\n Transcriptome\n -log10(p-Value)\n",
       legend = lgd_pval_RHCC44_24h,
       fill = colorRampPalette(colors =cols_Con)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
text(c(-27),c(0), srt=0, cex=2,labels="* not significant",adj=c(0,0))
text(c(15.5),c(20.5), srt=0, cex=2.5,labels="Position",adj=c(0.5,0))
#####
dev.off()
#####Transcriptome - Volcano - NCI-HCC44#####
ready_RH23_24h[ready_RH23_24h == "NaN"] <- NA
down_RH23_24h<-ready_RH23_24h[which(ready_RH23_24h[["pval"]] < 0.05&ready_RH23_24h[["log2FoldChange"]]<= 	-1),]
down_RH23_24h<-down_RH23_24h[order((down_RH23_24h[,4]*-1),decreasing = TRUE),]
down_RH23_24h$Nr <- 1:nrow(down_RH23_24h)

up_RH23_24h<-ready_RH23_24h[which(ready_RH23_24h[["pval"]] < 0.05&ready_RH23_24h[["log2FoldChange"]]  >= 	1),]
up_RH23_24h<-up_RH23_24h[order(up_RH23_24h[,4],decreasing = TRUE),]
up_RH23_24h$Nr <- 1:nrow(up_RH23_24h)
#####Transcriptome - Plotting Top10 - NCI-H23####
jpeg(filename = "Transcriptome_NCI-H23_Volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(8, 10, 6, 4),xpd=TRUE)
plot((ready_RH23_24h[,4]),-log10(ready_RH23_24h[,7]),pch=16,ylim=c(0,10), xlim = c(-10,10), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="RNA - NCI-H23 24h",mgp=c(7,2,0))
points(up_RH23_24h[,4],-log10(up_RH23_24h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_RH23_24h [,4],-log10(down_RH23_24h [,7]), col="#D53E4F",pch=16, cex=2)
lines(c(-1,-1),c(0,10.5), lty=3,lwd=3)
lines(c(1,1),c(0,10.5), lty=3,lwd=3)
lines(c(-11,11),c(1.3,1.3),lwd=3, lty=3)
lines(c(0,0),c(0,10),lty=2,col="grey",lwd=3)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -0.75,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -0.75,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=3)
title(ylab=bquote(-log[10]~"(pValue)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_RH23_24h[c(1),4])+3.8,(-log10(up_RH23_24h[c(1),7]))+0.5,labels=up_RH23_24h[c(1),1],cex=2.5)
text((up_RH23_24h[c(33),4])+4.8,(-log10(up_RH23_24h[c(33),7]))+0.5,labels=up_RH23_24h[c(33),1],cex=2.5)
text((down_RH23_24h[c(22),4])-3.35,(-log10(down_RH23_24h[c(22),7]))+0.50,labels=down_RH23_24h[c(22),1],cex=2.5)
text((down_RH23_24h[c(1),4])-3.75,(-log10(down_RH23_24h[c(1),7]))+0.50,labels=down_RH23_24h[c(1),1],cex=2.5)
dev.off()

#####Transcriptome - Top 10 - NCI-H23####
sig.pVal_RH23_24h_down<-ready_RH23_24h[which(ready_RH23_24h[["pval"]] < 0.05&ready_RH23_24h[["log2FoldChange"]] <(-1)),]
sig.pVal_RH23_24h_down<-sig.pVal_RH23_24h_down[order(sig.pVal_RH23_24h_down[,4], decreasing=FALSE),]
sig.pVal_RH23_24h_down<- sig.pVal_RH23_24h_down[1:10,]

sig.pVal_RH23_24h_up<-ready_RH23_24h[which(ready_RH23_24h[["pval"]] < 0.05&ready_RH23_24h[["log2FoldChange"]] >(1)),]
sig.pVal_RH23_24h_up<-sig.pVal_RH23_24h_up[order(sig.pVal_RH23_24h_up[,4], decreasing=TRUE),]
sig.pVal_RH23_24h_up <- sig.pVal_RH23_24h_up[1:10,]
sig.pVal_RH23_24h_up<- sig.pVal_RH23_24h_up[order(sig.pVal_RH23_24h_up[,4],decreasing=TRUE),]
#####
sig.con_pH23_24h_up<-ready_pH23_24h[which(ready_pH23_24h[["pval"]] < 0.05),]
sig.con_pH23_24h_up<-sig.con_pH23_24h_up[order(sig.con_pH23_24h_up[,5], decreasing=TRUE),]
sig.con_pH23_24h_up$rank <- 1:nrow(sig.con_pH23_24h_up)

sig.pVal_RH23_24h_up<-merge(x= sig.pVal_RH23_24h_up,y=sig.con_pH23_24h_up[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RH23_24h_up<-sig.pVal_RH23_24h_up[order(sig.pVal_RH23_24h_up[,5], decreasing=TRUE),]
#####
sig.con_pH23_24h_down<-ready_pH23_24h[which(ready_pH23_24h[["pval"]] < 0.05),]
sig.con_pH23_24h_down<-sig.con_pH23_24h_down[order(sig.con_pH23_24h_down[,5], decreasing=TRUE),]
sig.con_pH23_24h_down$rank <- 1:nrow(sig.con_pH23_24h_down)

sig.pVal_RH23_24h_down<-merge(x= sig.pVal_RH23_24h_down,y=sig.con_pH23_24h_down[,c(3,4,5,8,12)],by="ENSEMBL",all.x = TRUE)
sig.pVal_RH23_24h_down<-sig.pVal_RH23_24h_down[order(sig.pVal_RH23_24h_down[,5], decreasing=FALSE),]
#####
sig.pVal_RH23_24h_all <- rbind(sig.pVal_RH23_24h_up[10:1,c(2,5,8,11:14)],sig.pVal_RH23_24h_down[10:1,c(2,5,8,11:14)])
pval_RH23_24h<- data.frame(-log10(sig.pVal_RH23_24h_all[,3]))
pval_RH23_24h$order = findInterval(pval_RH23_24h[,1], sort(pval_RH23_24h[,1]))
v_pval_RH23_24h<- c(0,(max(pval_RH23_24h[,1])*(2/10)),(max(pval_RH23_24h[,1])*(4/10)),(max(pval_RH23_24h[,1])*(6/10)),(max(pval_RH23_24h[,1])*(8/10)))
pval_RH23_24h$order = findInterval(pval_RH23_24h[,1], v_pval_RH23_24h)
lgd_pval_RH23_24h = rep(NA, 5)
lgd_pval_RH23_24h[c(1,5)] = c(0,substring(max(pval_RH23_24h[,1]),  1,4))
#####Transcriptome - Plotting Top10 - NCI-H23#####
jpeg(filename = "Transcriptome_NCI-H23_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2, 1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(8, 10, 6, 4),xpd=TRUE)
sig.pVal_RH23_24h_all<-sig.pVal_RH23_24h_all[order(sig.pVal_RH23_24h_all[,2], decreasing = FALSE),]
plot(sig.pVal_RH23_24h_all[1:20,2],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_Con(5)[pval_RH23_24h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=3,line=2, cex.main=3,main="Transcriptome - NCI-H23 24h",outer=FALSE,adj=c(0.6))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
sig.pVal_RH23_24h_all<-sig.pVal_RH23_24h_all[order(sig.pVal_RH23_24h_all[,2], decreasing = FALSE),]
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.85),labels = sig.pVal_RH23_24h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.1),labels = sig.pVal_RH23_24h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_RH23_24h_all[i,2],0),c(i,i), lty=1,lwd=2)
}
#points(sig.pVal_RH23_24h_all[1:20,4],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_RH23_24h_2$order],cex=3,lwd=2)
for (i in 1:20) {
  lines(c(sig.pVal_RH23_24h_all[i,6],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  text(c(15.5),c(i,i), srt=0, cex=2.5,labels=ifelse(!is.na(sig.pVal_RH23_24h_all[i,6]), sig.pVal_RH23_24h_all[i,6],"*"))
}
legend(x = -28, y =12, title = "\n Transcriptome\n -log10(p-Value)\n",
       legend = lgd_pval_RH23_24h,
       fill = colorRampPalette(colors =cols_Con)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
text(c(-27),c(0), srt=0, cex=2,labels="* not significant",adj=c(0,0))
text(c(15.5),c(20.5), srt=0, cex=2.5,labels="Position",adj=c(0.5,0))
#####
dev.off()
setwd(source_wd)