setwd(save_wd)
#####Proteome - Volcano HCC44 0h####
ready_pHCC44_0h[ready_pHCC44_0h == "NaN"] <- NA
down_pHCC44_0h<-ready_pHCC44_0h[which(ready_pHCC44_0h[["pval"]] < 0.05&ready_pHCC44_0h[["diff"]]<= 	-1),]
down_pHCC44_0h<-down_pHCC44_0h[order((down_pHCC44_0h[,4]*-1),decreasing = TRUE),]
down_pHCC44_0h$Nr <- 1:nrow(down_pHCC44_0h)

up_pHCC44_0h<-ready_pHCC44_0h[which(ready_pHCC44_0h[["pval"]] < 0.05&ready_pHCC44_0h[["diff"]]  >= 	1),]
up_pHCC44_0h<-up_pHCC44_0h[order(up_pHCC44_0h[,4],decreasing = TRUE),]
up_pHCC44_0h$Nr <- 1:nrow(up_pHCC44_0h)
#####Proteome - Plot Volcano HCC44 0h#####
jpeg(filename = "Proteome_HCC44_0h_volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(10, 10, 6, 4),xpd=TRUE)
plot((ready_pHCC44_0h[,4]),-log10(ready_pHCC44_0h[,7]),pch=16,ylim=c(0,20), xlim = c(-20,20), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="Proteome - HCC44 0h",mgp=c(7,2,0))
points(up_pHCC44_0h[,4],-log10(up_pHCC44_0h[,7]), col= "#D53E4F",pch=16, cex=2)
points(down_pHCC44_0h [,4],-log10(down_pHCC44_0h [,7]), col= "#D53E4F",pch=16, cex=2)
lines(c(-1,-1),c(0,21), lty=3,lwd=3)
lines(c(1,1),c(0,21), lty=3,lwd=3)
lines(c(-21,21),c(1.3,1.3),lwd=3, lty=3)
lines(c(0,0),c(0,21),lty=3,col="grey60",lwd=2)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
axis(1,lwd=3,at = c(-20,-10,0,10,20),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10,15),lwd.ticks = 3,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=3)
title(ylab=bquote(-log[10]~"(p-Value)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_pHCC44_0h[c(1),4])+3,(-log10(up_pHCC44_0h[c(1),7]))+1,labels=up_pHCC44_0h[c(1),1],cex=2)
text((up_pHCC44_0h[c(2),4])+5.75,(-log10(up_pHCC44_0h[c(2),7]))+1.5,labels=up_pHCC44_0h[c(2),1],cex=2)
text((down_pHCC44_0h[c(1),4])-1,(-log10(down_pHCC44_0h[c(1),7]))+1.50,labels=down_pHCC44_0h[c(1),1],cex=2)
text((down_pHCC44_0h[c(2),4])-3.75,(-log10(down_pHCC44_0h[c(2),7]))+1.5,labels=down_pHCC44_0h[c(2),1],cex=2)
dev.off()

#####Proteom - Top 10 HCC44 0h#####
sig.pVal_pHCC44_0h_down<-ready_pHCC44_0h[which(ready_pHCC44_0h[["pval"]] < 0.05&ready_pHCC44_0h[["diff"]] <(-1)),]
sig.pVal_pHCC44_0h_down<-sig.pVal_pHCC44_0h_down[order(sig.pVal_pHCC44_0h_down[,4], decreasing=FALSE),]
sig.pVal_pHCC44_0h_down <- sig.pVal_pHCC44_0h_down[1:10,]

sig.pVal_pHCC44_0h_up<-ready_pHCC44_0h[which(ready_pHCC44_0h[["pval"]] < 0.05&ready_pHCC44_0h[["diff"]] >(1)),]
sig.pVal_pHCC44_0h_up<-sig.pVal_pHCC44_0h_up[order(sig.pVal_pHCC44_0h_up[,4], decreasing=TRUE),]
sig.pVal_pHCC44_0h_up <- sig.pVal_pHCC44_0h_up[1:10,]
sig.pVal_pHCC44_0h_up<- sig.pVal_pHCC44_0h_up[order(sig.pVal_pHCC44_0h_up[,4],decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_up<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_up<-sig.HCC44_RHCC44_24h_up[order(sig.HCC44_RHCC44_24h_up[,4], decreasing=TRUE),]
sig.HCC44_RHCC44_24h_up$rank <- 1:nrow(sig.HCC44_RHCC44_24h_up)

sig.pVal_pHCC44_0h_up<-merge(x= sig.pVal_pHCC44_0h_up,y=sig.HCC44_RHCC44_24h_up[,c(4,7,11,13)],by="ENSEMBL",all.x = TRUE)
#sig.pVal_pHCC44_0h_up[10,10:12]<-sig.HCC44_RHCC44_24h_up[which(sig.HCC44_RHCC44_24h_up$gene_name=="S100A13"),c(5,8,11)]
sig.pVal_pHCC44_0h_up<-sig.pVal_pHCC44_0h_up[order(sig.pVal_pHCC44_0h_up[,4], decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_down<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_down<-sig.HCC44_RHCC44_24h_down[order(sig.HCC44_RHCC44_24h_down[,4], decreasing=FALSE),]
sig.HCC44_RHCC44_24h_down$rank <- 1:nrow(sig.HCC44_RHCC44_24h_down)

sig.pVal_pHCC44_0h_down<-merge(x= sig.pVal_pHCC44_0h_down,y=sig.HCC44_RHCC44_24h_down[,c(4,7,11,13)],by="ENSEMBL",all.x = TRUE)
sig.pVal_pHCC44_0h_down<-sig.pVal_pHCC44_0h_down[order(sig.pVal_pHCC44_0h_down[,4], decreasing=FALSE),]
#####
sig.pVal_pHCC44_0h_all <- rbind(sig.pVal_pHCC44_0h_down[1:10, c(2,4,7,12:14)],sig.pVal_pHCC44_0h_up[, c(2,4,7,12:14)])
pval_pHCC44_0h<- data.frame(-log10(sig.pVal_pHCC44_0h_all[,3]))
pval_pHCC44_0h$order = findInterval(pval_pHCC44_0h[,1], sort(pval_pHCC44_0h[,1]))
v_pval_pHCC44_0h<- c(0,(max(pval_pHCC44_0h[,1],na.rm=TRUE)*(2/10)),(max(pval_pHCC44_0h[,1],na.rm=TRUE)*(4/10)),(max(pval_pHCC44_0h[,1],na.rm=TRUE)*(6/10)),(max(pval_pHCC44_0h[,1],na.rm=TRUE)*(8/10)))
pval_pHCC44_0h$order = findInterval(pval_pHCC44_0h[,1], v_pval_pHCC44_0h)
lgd_pval_pHCC44_0h = rep(NA, 5)
lgd_pval_pHCC44_0h[c(1,5)] = c(0,substring(max(pval_pHCC44_0h[,1]),  1,4))

pval_pHCC44_0h_2<- data.frame(-log10(sig.pVal_pHCC44_0h_all[,5]))
pval_pHCC44_0h_2$order = findInterval(pval_pHCC44_0h_2[,1], sort(pval_pHCC44_0h_2[,1]))
v_pval_pHCC44_0h_2<- c(0,(max(pval_pHCC44_0h_2[,1],na.rm=TRUE)*(2/10)),(max(pval_pHCC44_0h_2[,1],na.rm=TRUE)*(4/10)),(max(pval_pHCC44_0h_2[,1],na.rm=TRUE)*(6/10)),(max(pval_pHCC44_0h_2[,1],na.rm=TRUE)*(8/10)))
v_pval_pHCC44_0h_2<- v_pval_pHCC44_0h_2[complete.cases(v_pval_pHCC44_0h_2)]
pval_pHCC44_0h_2$order = findInterval(pval_pHCC44_0h_2[,2], v_pval_pHCC44_0h_2)
lgd_pval_pHCC44_0h_2 = rep(NA, 5)
lgd_pval_pHCC44_0h_2[c(1,5)] = c(0,substring(max(pval_pHCC44_0h_2[,1],na.rm=TRUE),  1,3))

#####Proteome- Plot Top10 HCC44 0h####
jpeg(filename = "Proteome_HCC44_0h_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(10, 14, 6, 4),xpd=TRUE)
plot(sig.pVal_pHCC44_0h_all[1:20,2],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_HCC44(5)[pval_pHCC44_0h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
title( cex.lab=3.5,cex.main=4,line=2, cex.main=3,main="Proteom - HCC44 0h",outer=FALSE,adj=c(0.9))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.65),labels = sig.pVal_pHCC44_0h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.2),labels = sig.pVal_pHCC44_0h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_0h_all[i,2],0),c(i,i), lty=1,lwd=2)
}
points(sig.pVal_pHCC44_0h_all[1:20,4],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_pHCC44_0h_2$order],cex=3,lwd=2)
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_0h_all[i,4],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  text(c(17),c(i,i), srt=0, cex=2,labels=ifelse(!is.na(sig.pVal_pHCC44_0h_all[i,6]), sig.pVal_pHCC44_0h_all[i,6],"*"))
}
legend(x = -32, y =12, title = "\n Proteome\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_0h,
       fill = colorRampPalette(colors =cols_HCC44)(5),
       border = "black",bg = NA ,
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
legend(x = -32, y =6.5, title = "\n Transcriptome\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_0h_2,
       fill = colorRampPalette(colors =cols_RAD)(5),
       border = "black",bg = NA ,
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
text(c(-32),c(0), srt=0, cex=1.5,labels="* not significant",adj=c(0,0))
text(c(17),c(20.5), srt=0, cex=2,labels="Position",adj=c(0.5,0))
#####
dev.off()
#####Proteome - Volcano HCC44 2h####
ready_pHCC44_2h[ready_pHCC44_2h == "NaN"] <- NA
down_pHCC44_2h<-ready_pHCC44_2h[which(ready_pHCC44_2h[["pval"]] < 0.05&ready_pHCC44_2h[["diff"]]<= 	-1),]
down_pHCC44_2h<-down_pHCC44_2h[order((down_pHCC44_2h[,4]*-1),decreasing = TRUE),]
down_pHCC44_2h$Nr <- 1:nrow(down_pHCC44_2h)

up_pHCC44_2h<-ready_pHCC44_2h[which(ready_pHCC44_2h[["pval"]] < 0.05&ready_pHCC44_2h[["diff"]]  >= 	1),]
up_pHCC44_2h<-up_pHCC44_2h[order(up_pHCC44_2h[,4],decreasing = TRUE),]
up_pHCC44_2h$Nr <- 1:nrow(up_pHCC44_2h)
#####Proteome - Plot Volcano HCC44 2h#####
jpeg(filename = "Proteome_HCC44_2h_Volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(10, 10, 6, 4),xpd=TRUE)
plot((ready_pHCC44_2h[,4]),-log10(ready_pHCC44_2h[,7]),pch=16,ylim=c(0,20), xlim = c(-20,20), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="Proteome - HCC44 2h",mgp=c(7,2,0))
points(up_pHCC44_2h[,4],-log10(up_pHCC44_2h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_pHCC44_2h [,4],-log10(down_pHCC44_2h [,7]), col="#D53E4F",pch=16, cex=2)
lines(c(-1,-1),c(0,21), lty=3,lwd=3)
lines(c(1,1),c(0,21), lty=3,lwd=3)
lines(c(-21,21),c(1.3,1.3),lwd=3, lty=3)
lines(c(0,0),c(0,21),lty=2,col="grey60",lwd=3)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 4,tcl = -1), y.args = list( lwd = 4,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10,15),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=4)
title(ylab=bquote(-log[10]~"(pValue)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_pHCC44_2h[c(1),4])+3,(-log10(up_pHCC44_2h[c(1),7]))+1,labels=up_pHCC44_2h[c(1),1],cex=2)
text((up_pHCC44_2h[c(2),4])+3.75,(-log10(up_pHCC44_2h[c(2),7]))+1.5,labels=up_pHCC44_2h[c(2),1],cex=2)
text((down_pHCC44_2h[c(2),4])-2,(-log10(down_pHCC44_2h[c(2),7]))+1.50,labels=down_pHCC44_2h[c(2),1],cex=2)
text((down_pHCC44_2h[c(1),4])-2.25,(-log10(down_pHCC44_2h[c(1),7]))+1.5,labels=down_pHCC44_2h[c(1),1],cex=2)
dev.off()

#####Proteome - Top 10 HCC44 2h#####
sig.pVal_pHCC44_2h_down<-ready_pHCC44_2h[which(ready_pHCC44_2h[["pval"]] < 0.05&ready_pHCC44_2h[["diff"]] <(-1)),]
sig.pVal_pHCC44_2h_down<-sig.pVal_pHCC44_2h_down[order(sig.pVal_pHCC44_2h_down[,4], decreasing=FALSE),]
sig.pVal_pHCC44_2h_down <- sig.pVal_pHCC44_2h_down[1:10,]

sig.pVal_pHCC44_2h_up<-ready_pHCC44_2h[which(ready_pHCC44_2h[["pval"]] < 0.05&ready_pHCC44_2h[["diff"]] >(1)),]
sig.pVal_pHCC44_2h_up<-sig.pVal_pHCC44_2h_up[order(sig.pVal_pHCC44_2h_up[,4], decreasing=TRUE),]
sig.pVal_pHCC44_2h_up <- sig.pVal_pHCC44_2h_up[1:10,]
sig.pVal_pHCC44_2h_up<- sig.pVal_pHCC44_2h_up[order(sig.pVal_pHCC44_2h_up[,4],decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_up<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_up<-sig.HCC44_RHCC44_24h_up[order(sig.HCC44_RHCC44_24h_up[,5], decreasing=TRUE),]
sig.HCC44_RHCC44_24h_up$rank <- 1:nrow(sig.HCC44_RHCC44_24h_up)

sig.pVal_pHCC44_2h_up<-merge(x= sig.pVal_pHCC44_2h_up,y=sig.HCC44_RHCC44_24h_up[,c(4,7,11,13)],by="ENSEMBL",all.x = TRUE)
#sig.pVal_pHCC44_2h_up[10,10:12]<-sig.HCC44_RHCC44_24h_up[which(sig.HCC44_RHCC44_24h_up$gene_name=="S100A13"),c(5,8,11)]
sig.pVal_pHCC44_2h_up<-sig.pVal_pHCC44_2h_up[order(sig.pVal_pHCC44_2h_up[,5], decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_down<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_down<-sig.HCC44_RHCC44_24h_down[order(sig.HCC44_RHCC44_24h_down[,4], decreasing=FALSE),]
sig.HCC44_RHCC44_24h_down$rank <- 1:nrow(sig.HCC44_RHCC44_24h_down)

sig.pVal_pHCC44_2h_down<-merge(x= sig.pVal_pHCC44_2h_down,y=sig.HCC44_RHCC44_24h_down[,c(4,7,11,13)],by="ENSEMBL",all.x = TRUE)
sig.pVal_pHCC44_2h_down<-sig.pVal_pHCC44_2h_down[order(sig.pVal_pHCC44_2h_down[,4], decreasing=FALSE),]
#####
sig.pVal_pHCC44_2h_all <- rbind(sig.pVal_pHCC44_2h_down[1:10, c(2,4,7,12:14)],sig.pVal_pHCC44_2h_up[, c(2,4,7,12:14)])
pval_pHCC44_2h<- data.frame(-log10(sig.pVal_pHCC44_2h_all[,3]))
pval_pHCC44_2h$order = findInterval(pval_pHCC44_2h[,1], sort(pval_pHCC44_2h[,1]))
v_pval_pHCC44_2h<- c(0,(max(pval_pHCC44_2h[,1])*(2/10)),(max(pval_pHCC44_2h[,1])*(4/10)),(max(pval_pHCC44_2h[,1])*(6/10)),(max(pval_pHCC44_2h[,1])*(8/10)))
pval_pHCC44_2h$order = findInterval(pval_pHCC44_2h[,1], v_pval_pHCC44_2h)
lgd_pval_pHCC44_2h = rep(NA, 5)
lgd_pval_pHCC44_2h[c(1,5)] = c(0,substring(max(pval_pHCC44_2h[,1]),  1,4))

pval_pHCC44_2h_2<- data.frame(-log10(sig.pVal_pHCC44_2h_all[,5]))
pval_pHCC44_2h_2$order = findInterval(pval_pHCC44_2h_2[,1], sort(pval_pHCC44_2h_2[,1]))
v_pval_pHCC44_2h_2<- c(0,(max(pval_pHCC44_2h_2[,1],na.rm=TRUE)*(2/10)),(max(pval_pHCC44_2h_2[,1],na.rm=TRUE)*(4/10)),(max(pval_pHCC44_2h_2[,1],na.rm=TRUE)*(6/10)),(max(pval_pHCC44_2h_2[,1],na.rm=TRUE)*(8/10)))
v_pval_pHCC44_2h_2<- v_pval_pHCC44_2h_2[complete.cases(v_pval_pHCC44_2h_2)]
pval_pHCC44_2h_2$order = findInterval(pval_pHCC44_2h_2[,2], v_pval_pHCC44_2h_2)
lgd_pval_pHCC44_2h_2 = rep(NA, 5)
lgd_pval_pHCC44_2h_2[c(1,5)] = c(0,substring(max(pval_pHCC44_2h_2[,1],na.rm=TRUE),  1,3))

#####Proteome - Plot Top10 HCC44 2h#####
jpeg(filename = "Proteome_HCC44_2h_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(10, 14, 6, 4),xpd=TRUE)
sig.pVal_pHCC44_2h_all<-sig.pVal_pHCC44_2h_all[order(sig.pVal_pHCC44_2h_all[,2], decreasing = FALSE),]
plot(sig.pVal_pHCC44_2h_all[1:20,2],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_HCC44(5)[pval_pHCC44_2h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=3,line=2, cex.main=3,main="Proteome - HCC44 2h",outer=FALSE,adj=c(0.7))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
sig.pVal_pHCC44_2h_all<-sig.pVal_pHCC44_2h_all[order(sig.pVal_pHCC44_2h_all[,2], decreasing = FALSE),]
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.65),labels = sig.pVal_pHCC44_2h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.2),labels = sig.pVal_pHCC44_2h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_2h_all[i,2],0),c(i,i), lty=1,lwd=2)
}
points(sig.pVal_pHCC44_2h_all[1:20,4],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_pHCC44_2h_2$order],cex=3,lwd=2)
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_2h_all[i,4],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  text(c(17),c(i,i), srt=0, cex=2,labels=ifelse(!is.na(sig.pVal_pHCC44_2h_all[i,6]), sig.pVal_pHCC44_2h_all[i,6],"*"))
}
legend(x = -32, y =12, title = "\n Proteome\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_2h,
       fill = colorRampPalette(colors =cols_HCC44)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
legend(x = -32, y =6.5, title = "\n Transcriptome\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_2h_2,
       fill = colorRampPalette(colors =cols_RAD)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
text(c(-32),c(0), srt=0, cex=1.5,labels="* not significant",adj=c(0,0))
text(c(17),c(20.5), srt=0, cex=2,labels="Position",adj=c(0.5,0))
#####
dev.off()
#####Proteome - Volcano HCC44 6h####
ready_pHCC44_6h[ready_pHCC44_6h == "NaN"] <- NA
down_pHCC44_6h<-ready_pHCC44_6h[which(ready_pHCC44_6h[["pval"]] < 0.05&ready_pHCC44_6h[["diff"]]<= 	-1),]
down_pHCC44_6h<-down_pHCC44_6h[order((down_pHCC44_6h[,4]*-1),decreasing = TRUE),]
down_pHCC44_6h$Nr <- 1:nrow(down_pHCC44_6h)

up_pHCC44_6h<-ready_pHCC44_6h[which(ready_pHCC44_6h[["pval"]] < 0.05&ready_pHCC44_6h[["diff"]]  >= 	1),]
up_pHCC44_6h<-up_pHCC44_6h[order(up_pHCC44_6h[,4],decreasing = TRUE),]
up_pHCC44_6h$Nr <- 1:nrow(up_pHCC44_6h)
#####Proteome - Plot Volcano HCC44 6h#####
jpeg(filename = "Proteome_HCC44_6h_volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(10, 10, 6, 4),xpd=TRUE)
plot((ready_pHCC44_6h[,4]),-log10(ready_pHCC44_6h[,7]),pch=16,ylim=c(0,20), xlim = c(-20,20), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="Proteome - HCC44 6h",mgp=c(7,2,0))
points(up_pHCC44_6h[,4],-log10(up_pHCC44_6h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_pHCC44_6h [,4],-log10(down_pHCC44_6h [,7]), col="#D53E4F",pch=16, cex=2)
lines(c(-1,-1),c(0,21), lty=3,lwd=3)
lines(c(1,1),c(0,21), lty=3,lwd=3)
lines(c(-21,21),c(1.3,1.3),lwd=3, lty=3)
lines(c(0,0),c(0,21),lty=2,col="grey60",lwd=3)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 4,tcl = -1), y.args = list( lwd = 4,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10,15),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=4)
title(ylab=bquote(-log[10]~"(pValue)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_pHCC44_6h[c(1),4])+4,(-log10(up_pHCC44_6h[c(1),7]))+1,labels=up_pHCC44_6h[c(1),1],cex=2)
text((up_pHCC44_6h[c(2),4])+5.75,(-log10(up_pHCC44_6h[c(2),7]))+1.5,labels=up_pHCC44_6h[c(2),1],cex=2)
text((down_pHCC44_6h[c(2),4])-3,(-log10(down_pHCC44_6h[c(2),7]))+1.50,labels=down_pHCC44_6h[c(2),1],cex=2)
text((down_pHCC44_6h[c(1),4])-3.25,(-log10(down_pHCC44_6h[c(1),7]))+1.5,labels=down_pHCC44_6h[c(1),1],cex=2)
dev.off()
#####Proteome - Top 10 HCC44 6h#####
sig.pVal_pHCC44_6h_down<-ready_pHCC44_6h[which(ready_pHCC44_6h[["pval"]] < 0.05&ready_pHCC44_6h[["diff"]] <(-1)),]
sig.pVal_pHCC44_6h_down<-sig.pVal_pHCC44_6h_down[order(sig.pVal_pHCC44_6h_down[,4], decreasing=FALSE),]
sig.pVal_pHCC44_6h_down <- sig.pVal_pHCC44_6h_down[1:10,]

sig.pVal_pHCC44_6h_up<-ready_pHCC44_6h[which(ready_pHCC44_6h[["pval"]] < 0.05&ready_pHCC44_6h[["diff"]] >(1)),]
sig.pVal_pHCC44_6h_up<-sig.pVal_pHCC44_6h_up[order(sig.pVal_pHCC44_6h_up[,4], decreasing=TRUE),]
sig.pVal_pHCC44_6h_up <- sig.pVal_pHCC44_6h_up[1:10,]
sig.pVal_pHCC44_6h_up<- sig.pVal_pHCC44_6h_up[order(sig.pVal_pHCC44_6h_up[,4],decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_up<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_up<-sig.HCC44_RHCC44_24h_up[order(sig.HCC44_RHCC44_24h_up[,4], decreasing=TRUE),]
sig.HCC44_RHCC44_24h_up$rank <- 1:nrow(sig.HCC44_RHCC44_24h_up)

sig.pVal_pHCC44_6h_up<-merge(x= sig.pVal_pHCC44_6h_up,y=sig.HCC44_RHCC44_24h_up[,c(4,7,11,13)],by="ENSEMBL",all.x = TRUE)
#sig.pVal_pHCC44_6h_up[10,10:12]<-sig.HCC44_RHCC44_24h_up[which(sig.HCC44_RHCC44_24h_up$gene_name=="S100A13"),c(5,8,11)]
sig.pVal_pHCC44_6h_up<-sig.pVal_pHCC44_6h_up[order(sig.pVal_pHCC44_6h_up[,4], decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_down<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_down<-sig.HCC44_RHCC44_24h_down[order(sig.HCC44_RHCC44_24h_down[,4], decreasing=FALSE),]
sig.HCC44_RHCC44_24h_down$rank <- 1:nrow(sig.HCC44_RHCC44_24h_down)

sig.pVal_pHCC44_6h_down<-merge(x= sig.pVal_pHCC44_6h_down,y=sig.HCC44_RHCC44_24h_down[,c(4,7,11,13)],by="ENSEMBL",all.x = TRUE)
sig.pVal_pHCC44_6h_down<-sig.pVal_pHCC44_6h_down[order(sig.pVal_pHCC44_6h_down[,4], decreasing=FALSE),]
#####
sig.pVal_pHCC44_6h_all <- rbind(sig.pVal_pHCC44_6h_down[1:10, c(2,4,7,12:14)],sig.pVal_pHCC44_6h_up[, c(2,4,7,12:14)])
pval_pHCC44_6h<- data.frame(-log10(sig.pVal_pHCC44_6h_all[,3]))
pval_pHCC44_6h$order = findInterval(pval_pHCC44_6h[,1], sort(pval_pHCC44_6h[,1]))
v_pval_pHCC44_6h<- c(0,(max(pval_pHCC44_6h[,1],na.rm=TRUE)*(2/10)),(max(pval_pHCC44_6h[,1],na.rm=TRUE)*(4/10)),(max(pval_pHCC44_6h[,1],na.rm=TRUE)*(6/10)),(max(pval_pHCC44_6h[,1],na.rm=TRUE)*(8/10)))
pval_pHCC44_6h$order = findInterval(pval_pHCC44_6h[,1], v_pval_pHCC44_6h)
lgd_pval_pHCC44_6h = rep(NA, 5)
lgd_pval_pHCC44_6h[c(1,5)] = c(0,substring(max(pval_pHCC44_6h[,1]),  1,4))

#pval_pHCC44_6h_2<- data.frame(-log10(sig.pVal_pHCC44_6h_all[,5]))
#pval_pHCC44_6h_2$order = findInterval(pval_pHCC44_6h_2[,1], sort(pval_pHCC44_6h_2[,1]))
#v_pval_pHCC44_6h_2<- c(0,(max(pval_pHCC44_6h_2[,1],na.rm=TRUE)*(2/10)),(max(pval_pHCC44_6h_2[,1],na.rm=TRUE)*(4/10)),(max(pval_pHCC44_6h_2[,1],na.rm=TRUE)*(6/10)),(max(pval_pHCC44_6h_2[,1],na.rm=TRUE)*(8/10)))
#v_pval_pHCC44_6h_2<- v_pval_pHCC44_6h_2[complete.cases(v_pval_pHCC44_6h_2)]
#pval_pHCC44_6h_2$order = findInterval(pval_pHCC44_6h_2[,2], v_pval_pHCC44_6h_2)
#lgd_pval_pHCC44_6h_2 = rep(NA, 5)
#lgd_pval_pHCC44_6h_2[c(1,5)] = c(0,substring(max(pval_pHCC44_6h_2[,1],na.rm=TRUE),  1,3))

#####Proteome - Plot Top 10 HCC44 6h####
jpeg(filename = "Proteome_HCC44_6h_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(10, 14, 6, 4),xpd=TRUE)
sig.pVal_pHCC44_6h_all<-sig.pVal_pHCC44_6h_all[order(sig.pVal_pHCC44_6h_all[,2], decreasing = FALSE),]
plot(sig.pVal_pHCC44_6h_all[1:20,2],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_HCC44(5)[pval_pHCC44_6h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=3,line=2, cex.main=3,main="Proteome - HCC44 6h",outer=FALSE,adj=c(0.6))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
sig.pVal_pHCC44_6h_all<-sig.pVal_pHCC44_6h_all[order(sig.pVal_pHCC44_6h_all[,2], decreasing = FALSE),]
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.65),labels = sig.pVal_pHCC44_6h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.2),labels = sig.pVal_pHCC44_6h_all[i,1],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_6h_all[i,2],0),c(i,i), lty=1,lwd=2)
}
#points(sig.pVal_pHCC44_6h_all[1:20,4],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_pHCC44_6h_2$order],cex=3,lwd=2)
#for (i in 1:20) {
#  lines(c(sig.pVal_pHCC44_6h_all[i,4],0),c(i,i), lty=1,lwd=2)
#}
#for (i in 1:20) {
#  text(c(16),c(i,i), srt=0, cex=2,labels=ifelse(!is.na(sig.pVal_pHCC44_6h_all[i,6]), sig.pVal_pHCC44_6h_all[i,6],"*"))
#}
legend(x = -32, y =12, title = "\n Proteome\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_6h,
       fill = colorRampPalette(colors =cols_HCC44)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
#legend(x = -32, y =6.5, title = "\n Transcriptome\n -log10(pValue)\n",
#       legend = lgd_pval_pHCC44_6h_2,
#       fill = colorRampPalette(colors =cols_RAD)(5),
#       border = "black",
#       y.intersp = 0.5,title.cex =2,
#       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
#text(c(-32),c(0), srt=0, cex=1.5,labels="* not significant",adj=c(0,0))
#text(c(17),c(20.5), srt=0, cex=2,labels="Position",adj=c(0.5,0))
#####
dev.off()
#####Proteome - Volcano HCC44 24h####
ready_pHCC44_24h[ready_pHCC44_24h == "NaN"] <- NA
down_pHCC44_24h<-ready_pHCC44_24h[which(ready_pHCC44_24h[["pval"]] < 0.05&ready_pHCC44_24h[["diff"]]<= 	-1),]
down_pHCC44_24h<-down_pHCC44_24h[order((down_pHCC44_24h[,4]*-1),decreasing = TRUE),]
down_pHCC44_24h$Nr <- 1:nrow(down_pHCC44_24h)

up_pHCC44_24h<-ready_pHCC44_24h[which(ready_pHCC44_24h[["pval"]] < 0.05&ready_pHCC44_24h[["diff"]]  >= 	1),]
up_pHCC44_24h<-up_pHCC44_24h[order(up_pHCC44_24h[,4],decreasing = TRUE),]
up_pHCC44_24h$Nr <- 1:nrow(up_pHCC44_24h)
#####Proteome - Plot Volcano HCC44 24h#####
jpeg(filename = "Proteome_HCC44_24h_Volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(10, 10, 6, 4),xpd=TRUE)
plot((ready_pHCC44_24h[,4]),-log10(ready_pHCC44_24h[,7]),pch=16,ylim=c(0,20), xlim = c(-20,20), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="Proteom - HCC44 24h",mgp=c(7,2,0))
points(up_pHCC44_24h[,4],-log10(up_pHCC44_24h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_pHCC44_24h [,4],-log10(down_pHCC44_24h [,7]), col="#D53E4F",pch=16, cex=2)
lines(c(-1,-1),c(0,21), lty=3,lwd=3)
lines(c(1,1),c(0,21), lty=3,lwd=3)
lines(c(-21,21),c(1.3,1.3),lwd=3, lty=3)
lines(c(0,0),c(0,21),lty=2,col="grey60",lwd=3)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 4,tcl = -1), y.args = list( lwd = 4,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10,15),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=4)
title(ylab=bquote(-log[10]~"(pValue)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_pHCC44_24h[c(1),4])+3,(-log10(up_pHCC44_24h[c(1),7]))+1,labels=up_pHCC44_24h[c(1),1],cex=2)
text((up_pHCC44_24h[c(2),4])+4.75,(-log10(up_pHCC44_24h[c(2),7]))+1.5,labels=up_pHCC44_24h[c(2),1],cex=2)
text((down_pHCC44_24h[c(3),4])-3.25,(-log10(down_pHCC44_24h[c(3),7]))+1.50,labels=down_pHCC44_24h[c(3),1],cex=2)
text((down_pHCC44_24h[c(1),4])-4.75,(-log10(down_pHCC44_24h[c(1),7]))+1,labels=down_pHCC44_24h[c(1),1],cex=2)
dev.off()
#####Proteome - Top 10 HCC44 24h#####
sig.pVal_pHCC44_24h_down<-ready_pHCC44_24h[which(ready_pHCC44_24h[["pval"]] < 0.05&ready_pHCC44_24h[["diff"]] <(-1)),]
sig.pVal_pHCC44_24h_down<-sig.pVal_pHCC44_24h_down[order(sig.pVal_pHCC44_24h_down[,4], decreasing=FALSE),]
sig.pVal_pHCC44_24h_down <- sig.pVal_pHCC44_24h_down[1:10,]

sig.pVal_pHCC44_24h_up<-ready_pHCC44_24h[which(ready_pHCC44_24h[["pval"]] < 0.05&ready_pHCC44_24h[["diff"]] >(1)),]
sig.pVal_pHCC44_24h_up<-sig.pVal_pHCC44_24h_up[order(sig.pVal_pHCC44_24h_up[,4], decreasing=TRUE),]
sig.pVal_pHCC44_24h_up <- sig.pVal_pHCC44_24h_up[1:10,]
sig.pVal_pHCC44_24h_up<- sig.pVal_pHCC44_24h_up[order(sig.pVal_pHCC44_24h_up[,4],decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_up<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_up<-sig.HCC44_RHCC44_24h_up[order(sig.HCC44_RHCC44_24h_up[,4], decreasing=TRUE),]
sig.HCC44_RHCC44_24h_up$rank <- 1:nrow(sig.HCC44_RHCC44_24h_up)

sig.pVal_pHCC44_24h_up<-merge(x= sig.pVal_pHCC44_24h_up,y=sig.HCC44_RHCC44_24h_up[,c(1,4,7,11,13)],by="ENSEMBL",all.x = TRUE)
sig.pVal_pHCC44_24h_up<-sig.pVal_pHCC44_24h_up[order(sig.pVal_pHCC44_24h_up[,4], decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_down<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_down<-sig.HCC44_RHCC44_24h_down[order(sig.HCC44_RHCC44_24h_down[,4], decreasing=FALSE),]
sig.HCC44_RHCC44_24h_down$rank <- 1:nrow(sig.HCC44_RHCC44_24h_down)

sig.pVal_pHCC44_24h_down<-merge(x= sig.pVal_pHCC44_24h_down,y=sig.HCC44_RHCC44_24h_down[,c(1,4,7,11,13)],by="ENSEMBL",all.x = TRUE)
sig.pVal_pHCC44_24h_down<-sig.pVal_pHCC44_24h_down[order(sig.pVal_pHCC44_24h_down[,5], decreasing=FALSE),]
#####
sig.pVal_pHCC44_24h_all <- rbind(sig.pVal_pHCC44_24h_down[1:10, c(1,2,4,7,11,13:15)],sig.pVal_pHCC44_24h_up[, c(1,2,4,7,11,13:15)])
pval_pHCC44_24h<- data.frame(-log10(sig.pVal_pHCC44_24h_all[,4]))
pval_pHCC44_24h$order = findInterval(pval_pHCC44_24h[,1], sort(pval_pHCC44_24h[,1]))
v_pval_pHCC44_24h<- c(0,(max(pval_pHCC44_24h[,1])*(2/10)),(max(pval_pHCC44_24h[,1])*(4/10)),(max(pval_pHCC44_24h[,1])*(6/10)),(max(pval_pHCC44_24h[,1])*(8/10)))
pval_pHCC44_24h$order = findInterval(pval_pHCC44_24h[,1], v_pval_pHCC44_24h)
lgd_pval_pHCC44_24h = rep(NA, 5)
lgd_pval_pHCC44_24h[c(1,5)] = c(0,substring(max(pval_pHCC44_24h[,1]),  1,4))

#pval_pHCC44_24h_2<- data.frame(-log10(sig.pVal_pHCC44_24h_all[,7]))
#pval_pHCC44_24h_2$order = findInterval(pval_pHCC44_24h_2[,1], sort(pval_pHCC44_24h_2[,1]))
#v_pval_pHCC44_24h_2<- c(0,(max(pval_pHCC44_24h_2[,1],na.rm=TRUE)*(2/10)),(max(pval_pHCC44_24h_2[,1],na.rm=TRUE)*(4/10)),(max(pval_pHCC44_24h_2[,1],na.rm=TRUE)*(6/10)),(max(pval_pHCC44_24h_2[,1],na.rm=TRUE)*(8/10)))
#v_pval_pHCC44_24h_2<- v_pval_pHCC44_24h_2[complete.cases(v_pval_pHCC44_24h_2)]
#pval_pHCC44_24h_2$order = findInterval(pval_pHCC44_24h_2[,2], v_pval_pHCC44_24h_2)
#lgd_pval_pHCC44_24h_2 = rep(NA, 5)
#lgd_pval_pHCC44_24h_2[c(1,5)] = c(0,substring(max(pval_pHCC44_24h_2[,1],na.rm=TRUE),  1,3))
#####Proteome - Plot Top 10 HCC44 24h#####
jpeg(filename = "Proteome_HCC44_24h_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 1,byrow = TRUE))
par(mar =c(10, 14, 6, 4),xpd=TRUE)
sig.pVal_pHCC44_24h_all<-sig.pVal_pHCC44_24h_all[order(sig.pVal_pHCC44_24h_all[,3], decreasing = FALSE),]
plot(sig.pVal_pHCC44_24h_all[1:20,3],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_HCC44(5)[pval_pHCC44_24h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=3,line=2, cex.main=3,main="Proteom - HCC44 24h",outer=FALSE,adj=c(0.6))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
#####
sig.pVal_pHCC44_24h_all<-sig.pVal_pHCC44_24h_all[order(sig.pVal_pHCC44_24h_all[,3], decreasing = FALSE),]
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.65),labels = sig.pVal_pHCC44_24h_all[i,2],srt = 0,cex = 2.5)
}
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.2),labels = sig.pVal_pHCC44_24h_all[i,2],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_24h_all[i,3],0),c(i,i), lty=1,lwd=2)
}
#points(sig.pVal_pHCC44_24h_all[1:20,6],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_pHCC44_24h_2$order],cex=3,lwd=2)
#for (i in 1:20) {
#  lines(c(sig.pVal_pHCC44_24h_all[i,4],0),c(i,i), lty=1,lwd=2)
#}
#for (i in 1:20) {
#  text(c(17),c(i,i), srt=0, cex=2,labels=ifelse(!is.na(sig.pVal_pHCC44_24h_all[i,8]), sig.pVal_pHCC44_24h_all[i,8],"*"))
#}
legend(x = -32, y =12, title = "\n Proteom\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_24h,
       fill = colorRampPalette(colors =cols_HCC44)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
#legend(x = -32, y =6.5, title = "\n Transcriptome\n -log10(pValue)\n",
#       legend = lgd_pval_pHCC44_24h_2,
#       fill = colorRampPalette(colors =cols_RAD)(5),
#       border = "black",
#       y.intersp = 0.5,title.cex =2,
#       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
#text(c(-32),c(0), srt=0, cex=1.5,labels="* not significant",adj=c(0,0))
#text(c(17),c(20.5), srt=0, cex=2,labels="Position",adj=c(0.5,0))
#####
dev.off()
#####Proteome - Volcano HCC44 48h####
ready_pHCC44_48h[ready_pHCC44_48h == "NaN"] <- NA
down_pHCC44_48h<-ready_pHCC44_48h[which(ready_pHCC44_48h[["pval"]] < 0.05&ready_pHCC44_48h[["diff"]]<= 	-1),]
down_pHCC44_48h<-down_pHCC44_48h[order((down_pHCC44_48h[,4]*-1),decreasing = TRUE),]
down_pHCC44_48h$Nr <- 1:nrow(down_pHCC44_48h)

up_pHCC44_48h<-ready_pHCC44_48h[which(ready_pHCC44_48h[["pval"]] < 0.05&ready_pHCC44_48h[["diff"]]  >= 	1),]
up_pHCC44_48h<-up_pHCC44_48h[order(up_pHCC44_48h[,4],decreasing = TRUE),]
up_pHCC44_48h$Nr <- 1:nrow(up_pHCC44_48h)
#####Proteome - Plot Volcano HCC44 48h#####
jpeg(filename = "Proteome_HCC44_48h_Volcano.jpeg",  width=20.5, height=14, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(12, 12, 6, 4),xpd=TRUE)
plot((ready_pHCC44_48h[,4]),-log10(ready_pHCC44_48h[,7]),pch=16,ylim=c(0,20), xlim = c(-20,20), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main="Proteome - HCC44 48h",mgp=c(7,2,0))
points(up_pHCC44_48h[,4],-log10(up_pHCC44_48h[,7]), col="#D53E4F",pch=16, cex=2)
points(down_pHCC44_48h [,4],-log10(down_pHCC44_48h [,7]), col="#D53E4F",pch=16, cex=2)
lines(c(-1,-1),c(0,21), lty=3,lwd=3)
lines(c(1,1),c(0,21), lty=3,lwd=3)
lines(c(-21,21),c(1.3,1.3),lwd=3, lty=3)
lines(c(0,0),c(0,21),lty=2,col="grey60",lwd=3)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 4,tcl = -1), y.args = list( lwd = 4,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,3,0))
axis(2,lwd=3,at = c(0,5,10,15),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -0.75,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=4)
title(ylab=bquote(-log[10]~"(pValue)"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text((up_pHCC44_48h[c(1),4])+3.5,(-log10(up_pHCC44_48h[c(1),7]))-1,labels=up_pHCC44_48h[c(1),1],cex=2)
text((up_pHCC44_48h[c(2),4])+3.75,(-log10(up_pHCC44_48h[c(2),7]))+1.5,labels=up_pHCC44_48h[c(2),1],cex=2)
text((down_pHCC44_48h[c(2),4])-2,(-log10(down_pHCC44_48h[c(2),7]))+1.50,labels=down_pHCC44_48h[c(2),1],cex=2)
text((down_pHCC44_48h[c(1),4])-3.25,(-log10(down_pHCC44_48h[c(1),7]))+1.5,labels=down_pHCC44_48h[c(1),1],cex=2)
dev.off()
#####Proteome - Top 10 HCC44 48h#####
sig.pVal_pHCC44_48h_down<-ready_pHCC44_48h[which(ready_pHCC44_48h[["pval"]] < 0.05&ready_pHCC44_48h[["diff"]] <(-1)),]
sig.pVal_pHCC44_48h_down<-sig.pVal_pHCC44_48h_down[order(sig.pVal_pHCC44_48h_down[,4], decreasing=FALSE),]
sig.pVal_pHCC44_48h_down <- sig.pVal_pHCC44_48h_down[1:10,]

sig.pVal_pHCC44_48h_up<-ready_pHCC44_48h[which(ready_pHCC44_48h[["pval"]] < 0.05&ready_pHCC44_48h[["diff"]] >(1)),]
sig.pVal_pHCC44_48h_up<-sig.pVal_pHCC44_48h_up[order(sig.pVal_pHCC44_48h_up[,4], decreasing=TRUE),]
sig.pVal_pHCC44_48h_up <- sig.pVal_pHCC44_48h_up[1:10,]
sig.pVal_pHCC44_48h_up<- sig.pVal_pHCC44_48h_up[order(sig.pVal_pHCC44_48h_up[,4],decreasing=FALSE),]
#####
sig.HCC44_RHCC44_24h_up<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_24h_up<-sig.HCC44_RHCC44_24h_up[order(sig.HCC44_RHCC44_24h_up[,4], decreasing=TRUE),]
sig.HCC44_RHCC44_24h_up$rank <- 1:nrow(sig.HCC44_RHCC44_24h_up)

sig.pVal_pHCC44_48h_up<-merge(x= sig.pVal_pHCC44_48h_up,y=sig.HCC44_RHCC44_24h_up[,c(1,4,7,11,13)],by="ENSEMBL",all.x = TRUE)
sig.pVal_pHCC44_48h_up<-sig.pVal_pHCC44_48h_up[order(sig.pVal_pHCC44_24h_up[,4], decreasing=FALSE),]
#####
sig.HCC44_RHCC44_48h_down<-ready_RHCC44_24h[which(ready_RHCC44_24h[["pval"]] < 0.05),]
sig.HCC44_RHCC44_48h_down<-sig.HCC44_RHCC44_48h_down[order(sig.HCC44_RHCC44_48h_down[,4], decreasing=FALSE),]
sig.HCC44_RHCC44_48h_down$rank <- 1:nrow(sig.HCC44_RHCC44_48h_down)

sig.pVal_pHCC44_48h_down<-merge(x= sig.pVal_pHCC44_48h_down,y=sig.HCC44_RHCC44_48h_down[,c(1,4,7,11,13)],by="ENSEMBL",all.x = TRUE)
sig.pVal_pHCC44_48h_down<-sig.pVal_pHCC44_48h_down[order(sig.pVal_pHCC44_48h_down[,5], decreasing=FALSE),]
#####
sig.pVal_pHCC44_48h_all <- rbind(sig.pVal_pHCC44_48h_down[1:10, c(1,2,4,7,11,13:15)],sig.pVal_pHCC44_48h_up[, c(1,2,4,7,11,13:15)])
pval_pHCC44_48h<- data.frame(-log10(sig.pVal_pHCC44_48h_all[,4]))
pval_pHCC44_48h$order = findInterval(pval_pHCC44_48h[,1], sort(pval_pHCC44_48h[,1]))
v_pval_pHCC44_48h<- c(0,(max(pval_pHCC44_48h[,1])*(2/10)),(max(pval_pHCC44_48h[,1])*(4/10)),(max(pval_pHCC44_48h[,1])*(6/10)),(max(pval_pHCC44_48h[,1])*(8/10)))
pval_pHCC44_48h$order = findInterval(pval_pHCC44_48h[,1], v_pval_pHCC44_48h)
lgd_pval_pHCC44_48h = rep(NA, 5)
lgd_pval_pHCC44_48h[c(1,5)] = c(0,substring(max(pval_pHCC44_48h[,1]),  1,4))

pval_pHCC44_48h_2<- data.frame(-log10(sig.pVal_pHCC44_48h_all[,7]))
pval_pHCC44_48h_2$order = findInterval(pval_pHCC44_48h_2[,1], sort(pval_pHCC44_48h_2[,1]))
v_pval_pHCC44_48h_2<- c(0,(max(pval_pHCC44_48h_2[,1],na.rm=TRUE)*(2/10)),(max(pval_pHCC44_48h_2[,1],na.rm=TRUE)*(4/10)),(max(pval_pHCC44_48h_2[,1],na.rm=TRUE)*(6/10)),(max(pval_pHCC44_48h_2[,1],na.rm=TRUE)*(8/10)))
v_pval_pHCC44_48h_2<- v_pval_pHCC44_48h_2[complete.cases(v_pval_pHCC44_48h_2)]
pval_pHCC44_48h_2$order = findInterval(pval_pHCC44_48h_2[,2], v_pval_pHCC44_48h_2)
lgd_pval_pHCC44_48h_2 = rep(NA, 5)
lgd_pval_pHCC44_48h_2[c(1,5)] = c(0,substring(max(pval_pHCC44_48h_2[,1],na.rm=TRUE),  1,3))
#####Proteome - Plot Top10 HCC44 48h####
jpeg(filename = "Proteome_HCC44_48h_Top10.jpeg",  width=20.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,1,2),ncol=3,nrow = 2,byrow = TRUE))
par(mar =c(10, 14, 6, 4),xpd=TRUE)
sig.pVal_pHCC44_48h_all<-sig.pVal_pHCC44_48h_all[order(sig.pVal_pHCC44_48h_all[,3], decreasing = FALSE),]
plot(sig.pVal_pHCC44_48h_all[1:20,3],1:20,xlim=c(-15,15),ylim=c(0,20),pch=21,ann=FALSE, axes=FALSE,col="black",bg=pal_HCC44(5)[pval_pHCC44_48h$order],cex=3,lwd=2)
lines(c(0,0),c(0,20), lty=1,lwd=3)
minor.tick(nx =6/3,ny=1,tick.ratio =1, x.args = list( lwd = 3,tcl = -1), y.args = list( lwd = 3,tcl = -1)) 
axis(1,lwd=3,at = c(-10,0,10),lwd.ticks = 4,cex=4, cex.axis=2.5,tcl = -1.5,mgp=c(3,3,0))
#box(which = "plot",lwd=4)
title( cex.lab=3.5,cex.main=3,line=2, cex.main=3,main="Proteome - HCC44 48h",outer=FALSE,adj=c(0.6))
title(xlab=bquote(log[2]~"(2Gy/Control)"),cex.lab=3,cex.main=3,line=7,adj=c(0.5))
sig.pVal_pHCC44_48h_all<-sig.pVal_pHCC44_48h_all[order(sig.pVal_pHCC44_48h_all[,3], decreasing = FALSE),]
for (i in 1:10) {
  text(x = 1,y = c(i),adj=c(0,0.65),labels = sig.pVal_pHCC44_48h_all[i,2],srt = 0,cex = 2.5)
}
for (i in 11:20) {
  text(x = -1,y = c(i),adj=c(1,0.2),labels = sig.pVal_pHCC44_48h_all[i,2],srt = 0,cex = 2.5)
}
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_48h_all[i,3],0),c(i,i), lty=1,lwd=2)
}
points(sig.pVal_pHCC44_48h_all[1:20,6],(1:20)+0.25,pch=21,col="black",bg=pal_RAD(5)[pval_pHCC44_48h_2$order],cex=3,lwd=2)
for (i in 1:20) {
  lines(c(sig.pVal_pHCC44_48h_all[i,6],0),c(i,i), lty=1,lwd=2)
}
for (i in 1:20) {
  text(c(17),c(i,i), srt=0, cex=1.5,labels=ifelse(!is.na(sig.pVal_pHCC44_48h_all[i,8]), sig.pVal_pHCC44_48h_all[i,8],"*"))
}
legend(x = -32, y =12, title = "\n Proteome\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_48h,
       fill = colorRampPalette(colors =cols_HCC44)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
legend(x = -32, y =6.5, title = "\n Transcriptome\n -log10(pValue)\n",
       legend = lgd_pval_pHCC44_48h_2,
       fill = colorRampPalette(colors =cols_RAD)(5),
       border = "black",
       y.intersp = 0.5,title.cex =2,
       cex = 2,pt.cex=2,box.lwd=1,box.col = "white")
text(c(-32),c(0), srt=0, cex=1.5,labels="* not significant",adj=c(0,0))
text(c(17),c(20.5), srt=0, cex=1.75,labels="Position",adj=c(0.5,0))
#####
dev.off()
setwd(source_wd)
