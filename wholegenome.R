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
setwd("/home/xinzhou")
write.table(contig_table,"contig_table.txt")
id1<-character()
id2<-character()
contig1<-integer()
contig2<-integer()
N<-integer()
A<-data.frame(id1,id2,contig1,contig2,N,stringsAsFactors=FALSE)
n<-length(chr)
for (j in 1:(n-1)) {
  for (i in (j+1):n) {
    B<-intersect(chr[[i]][,9],chr[[j]][,9])
    B1<-B%%(10^8)
    B2<-B%/%(10^8)
    A<-rbind(A,cbind(B1,B2,rep(j,length(B1)),rep(i,length(B2)),rep(1,length(B1))))
    
  }
}
write.table(A,"wholegenome.txt")