#####loading data####
setwd("C:/Users/annar/Desktop/Codes+Data/GEO")

Rdata <- data.frame(read.table("NormData_Bohnenberger.csv", sep="\t",  header=TRUE))
Rdata$ID <- 1:nrow(Rdata) 
identifier<-checkGeneSymbols(Rdata[,2])
Rdata<- cbind(Rdata, identifier)

Rdata[,9:23][Rdata[,9:23]==0]<-NA
#####loading raw Data####
setwd("C:/Users/annar/Desktop/Codes+Data/Pride/results")
proteome_data <- read.csv("Intensities_all.csv", header=TRUE,sep="\t")
proteome_data[,7:90] <- lapply(proteome_data[,7:90], function(x) {
  x[x == "NaN"] <- NA           
  as.numeric(x)                
})
proteome_data$ID <- 1:nrow(proteome_data) 
identifier<-checkGeneSymbols(proteome_data[,2],unmapped.as.na=TRUE)
proteome_data<- cbind(proteome_data, identifier)
proteome_data<-merge(x=proteome_data,y=res,by.x="Suggested.Symbol" ,by.y="SYMBOL",all.x=TRUE)
proteome_data<-proteome_data[!duplicated(proteome_data[,92]),]
proteome_data<-proteome_data[,c(2:92,1,94:97)]
#####
ready_pbase_0h<-read.csv("stat_pCon_0h.csv")
ready_pRAD_24h<-read.csv("stat_pRAD_24h.csv")
ready_pH23_24h<-read.csv("stat_pH23_24h.csv")
ready_pHCC44_24h<-read.csv("stat_pHCC44_24h.csv")
setwd(source_wd)