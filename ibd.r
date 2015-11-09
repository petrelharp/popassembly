#The script calculating the ibd shared on chrom 22 and how adjacent contigs share more ibds 
#read the data and the country of each individual, we have 69 contigs this time
country<-read.csv("sample-info1.csv")
inds<-read.table("tmp.ids",head=F)
c<-list()
l<-matrix(0,1,69)
for (i in 1:69) {
  name<-paste(i,"result.ibd",sep="")
  if (file.exists(name)) {
    if (file.info(name)$size!=0) {
      c[[i]]<-read.table(name,head=F)
      l[i]<-nrow(c[[i]])
    }
    else l[i]<-0
  }
  else l[i]<-0
}
country.table<-country[,c(2,11)]
country.table<-subset(country.table,SUBJID%in%inds[,1])
country.table<-subset(country.table,!is.na(COUNTRY_SELF))
idtocountry<-function(pair,table) {
  c1<-table[match(pair[1],table[,1]),2]
  c2<-table[match(pair[2],table[,1]),2]
  c(c1,c2,pair[3],pair[4],pair[5])
}
#Here we find that many contigs have none all few ibds, so we merge the contigs into 16 large contigs according to the number of ibds
m<-l[l>0]
p<-match(m,l)
d<-list()
d[[1]]<-rbind(c[[p[1]]],c[[p[2]]])
d[[2]]<-rbind(c[[p[3]]],c[[p[4]]],c[[p[5]]])
d[[3]]<-c[[p[6]]]
d[[4]]<-c[[p[7]]]
d[[5]]<-c[[p[8]]]
d[[6]]<-rbind(c[[p[9]]],c[[p[10]]],c[[p[11]]],c[[p[12]]])
d[[7]]<-c[[p[13]]]
d[[8]]<-rbind(c[[p[14]]],c[[p[15]]])
d[[9]]<-c[[p[16]]]
d[[10]]<-c[[p[17]]]
d[[11]]<-rbind(c[[p[18]]],c[[p[19]]],c[[p[20]]],c[[p[21]]])
d[[12]]<-rbind(c[[p[22]]],c[[p[23]]])
d[[13]]<-c[[p[24]]]
d[[14]]<-c[[p[25]]]
d[[15]]<-c[[p[26]]]
d[[16]]<-rbind(c[[p[27]]],c[[p[28]]],c[[p[29]]])
id1<-character()
id2<-character()
contig1<-integer()
contig2<-integer()
N<-integer()
A<-data.frame(id1,id2,contig1,contig2,N,stringsAsFactors=FALSE)


#calculate the number of ibd shared by each pair of contigs and individuals
N<-matrix(0,69,69)
for (j in 1:28) {
  for (i in (j+1):29) {
    left<-c[[p[j]]][,c(1,3)]
    right<-c[[p[i]]][,c(1,3)]
    for (k in 1:nrow(left)) {
      pair<-c(left[k,1],left[k,2])
      compare.1<-right[,1]==pair[1]
      compare.2<-right[,2]==pair[2]
      compare.3<-right[,1]==pair[2]
      compare.4<-right[,2]==pair[1]
      if (sum(compare.1*compare.2)+sum(compare.3*compare.4)>0) {
        A[nrow(A)+1,]<-c(pair[1],pair[2],j,i,1)
      }
    }    
  }
}

country1<-country.table[match(A[,1],country.table[,1]),2]
country2<-country.table[match(A[,2],country.table[,1]),2]

B.new<-cbind(country1,country2,A[,c(3:5)])
colnames(N)<-c("id1","id2","contig1","contig2","N")
B<-as.data.frame(B)

N<-aggregate(B.new$N,by=list(B.new$country1,B.new$country2,B.new$contig1,B.new$contig2),FUN="sum")
M<-matrix(0,29,29)
for (j in 1:28) {
  for (i in (j+1):29) {
    M[i,j]<-N[p[i],p[j]]
    
  }
}
Z <- (M-mean(M))*29/sqrt(sum((M-mean(M))^2))
matplot(Z,pch=1,col="black")
for (i in 1:28) {
  points(i,Z[i+1,i],col="red")
}
dx<-row(Z)-col(Z)
dx[abs(dx)>10]<-NA
table(dx,apply(Z,1,rank))
u<-upper.tri(M,diag=FALSE)
M[u]<-t(M)[u]
for (i in 1:29) {
  for (j in 1:29) {
    if (sum(M[i,])*sum(M[,j])!=0) {
      Z[i,j]<-(M[i,j]-(sum(M[i,])*sum(M[,j]))/sum(M))/(sum(M[i,])*sum(M[,j])/sum(M))
    } else {
      Z[i,j]<-0
    }
  }
}
P<-matrix(0,29,29)
for (i in 1:29) {
  for (j in 1:29) {
    
      P[i,j]<-dbinom(M[i,j],size=sum(M),prob=(sum(M[i,])*sum(M[,j]))/(sum(M)*sum(M)))
    
  }
}

table <- BuildIBDTable(10, 29, N)

the.glm <- glm(table$N~table$contig1+table$contig2+table$id1*table$id2,family=poisson("log"))
resid(the.glm)->resid
cbind(table,resid)->N.resid
N.adja<-subset(N.resid,(as.numeric(contig2)-as.numeric(contig1)==0))
N.resid$resid2 <- N.resid$N-fitted(the.glm)
tmp <- ( tapply( N.resid$resid, N.resid[c("contig1","contig2")], mean, na.rm=TRUE))
tmp2 <- ( tapply( N.resid$resid2, N.resid[c("contig1","contig2")], mean, na.rm=TRUE))
plot( row(tmp)-col(tmp), as.vector(tmp2), xlim=c(-30,0), col=row(tmp), pch=col(tmp)%%23+1)