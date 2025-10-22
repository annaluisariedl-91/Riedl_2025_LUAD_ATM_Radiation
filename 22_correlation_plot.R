#####correlation - preparing data####
compare_H23 <- merge(x=ready_pH23_24h[,c(1,4,7)],y=ready_RH23_24h[,c(1,4,7)], by.x="Suggested.Symbol",by.y="GeneSymbol")
compare_H23 <- compare_H23[complete.cases(compare_H23),]
compare_H23 <- compare_H23[compare_H23$Suggested.Symbol != "", ]
compare_H23 <- compare_H23 %>%  mutate(SignifikanzStatus = case_when(
  pval.x < 0.05 & pval.y < 0.05 ~ 3,
  pval.x < 0.05 ~ 1,
  pval.y < 0.05 ~ 2,
  TRUE ~ 4))

col_h23<- c("orange2","#5E4FA2", "#66C2A5", "black")[compare_H23$SignifikanzStatus]
corr_H23 <- cor(compare_H23[,2],compare_H23[,4], method = "pearson")
compare_HCC44 <- merge(x=ready_pHCC44_24h[,c(1,4,7)],y=ready_RHCC44_24h[,c(1,4,7)], by.x="Suggested.Symbol",by.y="GeneSymbol")
compare_HCC44 <- compare_HCC44[complete.cases(compare_HCC44),]
compare_HCC44 <- compare_HCC44[compare_HCC44$Suggested.Symbol != "", ]
compare_HCC44 <- compare_HCC44 %>%  mutate(SignifikanzStatus = case_when(
  pval.x < 0.05 & pval.y < 0.05 ~ 3,
  pval.x < 0.05 ~ 1,
  pval.y < 0.05 ~ 2,
  TRUE ~ 4))

col_HCC44<- c("orange2","#5E4FA2", "#66C2A5", "black")[compare_HCC44$SignifikanzStatus]
corr_HCC44 <- cor(compare_HCC44[,2],compare_HCC44[,4], method = "pearson")
#####Plotting correlation####
setwd(save_wd)
jpeg(filename = "Correlation_Proteom_Trancriptome.jpeg",  width=61.5, height=21, unit="cm",quality=100,res=70)
layout(matrix(c(1,2),ncol=2,nrow = 1,byrow = TRUE))
par(mar = c(10, 14, 8, 4),xpd=TRUE)
leg.txt <- c("significant in proteome","significant in transcriptome","significant in both")
plot((compare_HCC44[,2]),compare_HCC44[,4],pch=16, col=col_HCC44,ylim=c(-5,5), xlim = c(-10,10), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main=,mgp=c(7,2,0))
lines(c(-10,10),c(0,0),lty=3,col="grey40",lwd=2)
lines(c(0,0),c(-5,5),lty=3,col="grey40",lwd=2)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd =3,tcl = -0.75)) 
axis(1,lwd=3,at = c(-10,-5,0,5,10),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -1.5,mgp=c(3,3,0))
axis(2,lwd=3,at = c(-5,0,5),lwd.ticks = 4,cex=4, cex.axis=3,tcl = -1.5,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=3)
title(main="Comparison Proteome vs Transcriptome \n HCC44",cex.lab=3,cex.main=3,line=2,adj=c(0.5))
title(ylab=bquote("Proteome["~log[2]~"(2Gy/Control)]"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote("Proteome["~log[2]~"(2Gy/Control)]"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text(c(2),c(4), srt=0, cex=2.5,labels= bquote(R^2 == .(round(summary(lm(compare_HCC44[,2] ~ compare_HCC44[,4]))$r.squared, 5))),adj=c(0,0))
legend(x=c(-10.5,0),y=c(-5,-2),leg.txt,ncol=1,cex=2,pch=19,bty="n",  col=c("orange2","#5E4FA2", "#66C2A5", "black"),inset=0.02,box.lwd =2,x.intersp = 0.7)
model <- lm(compare_HCC44[,2] ~ compare_HCC44[,4])
clip(x1 = -10, x2 = 10, y1 = -5, y2 = 5)
abline(model, col ="#D53E4F", lwd = 4)

plot((compare_H23[,2]),compare_H23[,4],pch=16, col=col_h23,ylim=c(-5,5), xlim = c(-10,10), cex.axis=4,cex.lab=4,cex=2,
     xlab="", cex.main=3,ylab = "",axes=FALSE, main=,mgp=c(7,2,0))
lines(c(-10,10),c(0,0),lty=3,col="grey40",lwd=2)
lines(c(0,0),c(-5,5),lty=3,col="grey40",lwd=2)
minor.tick(nx =6/3,ny=2,tick.ratio =1, x.args = list( lwd = 3,tcl = -0.75), y.args = list( lwd = 3,tcl = -0.75)) 
axis(1,lwd=3,at = c(-10,-5,0,5,10),lwd.ticks = 3,cex=2.5, cex.axis=3,tcl = -1.5,mgp=c(3,3,0))
axis(2,lwd=3,at = c(-5,0,5),lwd.ticks = 3,cex=2.5, cex.axis=3,tcl = -1.5,mgp=c(3,2,0),las=2)
box(which = "plot",lwd=3)
title(main="Comparison Proteome vs Transcriptome \n NCI-H23",cex.lab=3,cex.main=3,line=2,adj=c(0.5))
title(ylab=bquote("Proteome["~log[2]~"(2Gy/Control)]"),cex.lab=3,cex.main=3,line=6,adj=c(0.5))
title(xlab=bquote("Proteome["~log[2]~"(2Gy/Control)]"),cex.lab=3,cex.main=3,line=7,adj=c(0.6))
text(c(2),c(4), srt=0, cex=2.5,labels= bquote(R^2 == .(round(summary(lm(compare_H23[,2] ~ compare_H23[,4]))$r.squared, 5))),adj=c(0,0))
legend(x=c(-10.5,0),y=c(-5,-2),leg.txt,ncol=1,cex=2,pch=19,bty="n",  col=c("orange2","#5E4FA2", "#66C2A5", "black"),inset=0.02,box.lwd =3,x.intersp = 0.7)
model <- lm(compare_H23[,2] ~ compare_H23[,4])
clip(x1 = -10, x2 = 10, y1 = -5, y2 = 5)
abline(model, col ="#D53E4F", lwd = 4)
dev.off()
setwd(source_wd)