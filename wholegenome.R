 #This file is to calculate the ibd shared by any two contigs of the whole genome
#read all the ibd.sorted.file
setwd("/home/xinzhou/ibdsorted")
chr<-list()
id<-integer()
name<-character()
contig_table<-data.frame(id,name,stringsAsFactors=FALSE)
f<-list.files(pattern="ibd.sorted")
for (name in f) {
  if (file.exists(name)) {
    if (file.info(name)$size!=0) {
      chr[[length(chr)+1]]<-read.table(name,head=F)
      contig_table[nrow(contig_table)+1,]<-c(nrow(contig_table)+1,name)
    }
  }
}

#calculate the ibd shared by any two contigs
setwd("/home/xinzhou")
write.table(contig_table,"contig_table.txt")
#read the file that records the country of each individual
country<-read.csv("sample-info1.csv")
inds<-read.table("tmp.ids",head=F)
setwd("/home/xinzhou/ibdshared")
country.table <- country[,c(2,11)]
country.table <- subset(country.table,SUBJID%in%inds[,1]&(!is.na(COUNTRY_SELF)))
n <- length(chr)
for (j in 1:(n-1)) {
  for (i in (j+1):n) {
    B <- intersect(chr[[i]][,9],chr[[j]][,9])
    B1 <- as.integer(B%%(10^8))
    B2 <- as.integer(B%/%(10^8))
    #change id to countries and sum up the ibd shared between same pair of countries
    country1 <- country.table[match(B1,country.table[,1]),2]
    country2 <- country.table[match(B2,country.table[,1]),2]
    B.new <- data.frame(country1,country2,contig1=rep(j,length(B1)),contig2=rep(i,length(B1)),N=rep(1,length(B1)))
    if (nrow(B.new)>0) {
      N.agg <- aggregate(B.new$N,by=list(B.new$country1,B.new$country2,B.new$contig1,B.new$contig2),FUN="sum")
      write.table(N.agg,paste("contig_",j,"_contig_",i,"_ibdshared.csv",sep=""),row.names=FALSE,sep=',')  
    }
  }
}