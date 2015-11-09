#This part simulates the process that if two contigs share one ibd, the ibd is devided to two parts. 
#We assume that the length of each part follows poisson distribution exp(T), where T follows a Gamma distribution Gamma(2,5*exp(-4)). 
#How many of the ibds have both parts of length above 0.02.
#start to simulate data, we focus on chrom 1, so we assume there are 262 contigs(1mb on average)
P <- c(1000,5000,10000,20000)
Nreps <- 5
Simulate <- function(ncontig,nibd) {
  T <- matrix(rgamma(nibd*ncontig,2,0.0005),nibd)
  L <- matrix(rexp(nibd*ncontig,T),nibd)
  R <- matrix(rexp(nibd*ncontig,T),nibd)
  LL <- (sign(L-0.02)+1)/2
  RR <- (sign(R-0.02)+1)/2
  N <- t(RR)%*%LL
  nn <- sum(N)
}
W.rep <- replicate(Nreps, sapply(P,Simulate,ncontig=262) )

setwd("/Users/xzge/百度云同步盘/dy")
#start to deal with human data
#read the file about ibd which contains the start and end position, id and country of each pair of person and the geographic distance between them
a <- read.csv("blocks-for-xinzhou.csv",head=TRUE)
b <- a[,-6]
b <- b[,-6]

c <- list() # divide the data according to chromesome number
for (j in 1:22){
  c[[j]] <- b[b$chrom==j,]
}

l <- matrix(0,1,10) #the length of the chromesome
for (i in 1:10) {
  l[i] <- max(c[[i]][,5])
}
# the function calculating Nij with given length restriction r and contigs
# The return value has 
#   N[i,j] = the number of pairs 
#            that have an IBD block overlapping the *left* end of contig *i* by at least distance select.r
#            *and* an IBD block overlapping the *right* end of contig *j* by at least distance select.r
# So: IBD blocks that span the break between contig i and i+1 will be counted in N[i+1,i] 

OverLapsCountAll <- function(data.mapstart,data.mapend,contig.div,select.r){
  OverLaps.left <- sapply(contig.div, function (cdiv) {
    (data.mapend > cdiv) & (data.mapstart < cdiv) & (data.mapend-cdiv > select.r)
  } )
  OverLaps.right <- sapply(contig.div, function (cdiv) {
    (data.mapend > cdiv) & (data.mapstart < cdiv) & (cdiv-data.mapstart>select.r)
  } )
  return( t(OverLaps.left)%*%OverLaps.right )
}
#repete the process 20 times at 11 different r 
d <- c[[1]]
indivs <- sort(unique(c(a$id1,a$id2)))
select.r <- seq(1,3,0.2)
Nreps <- 20
Nindivs <- 400
Ncontig <- 262

#This function calculates W with given ibd information and r by sampling 400 individuals from all, and divide the chromosome into contigs. 
#W implies how adjacent contigs share more ibd blocks
WAccordingToR <- function(indivs,select.r) {
    contig.div <- cumsum( rexp(1*Ncontig,1) )#sample from the dateset
    # USE THESE to restrict to contigs at least of length select.r
    sub.indivs <- sample( indivs, Nindivs )
    d.sub <- subset(d, ( id1 %in% sub.indivs ) & ( id2 %in% sub.indivs ) )
    pair.factor <- with( d.sub, paste( ifelse( id1<id2, id1, id2 ),
                                       ifelse( id1<id2, id2, id1 ), sep="-" ) )
    data.mapstart <- d.sub$mapstart
    data.mapend <- d.sub$mapend
    N.pairs <- OverLapsCountAll(data.mapstart,data.mapend,contig.div,select.r)
    # calculate w
    Z <- (N.pairs-mean(N.pairs))*Ncontig/sqrt(sum((N.pairs-mean(N.pairs))^2))
    W <- sum( Z[ (row(Z)+1)==col(Z)] )
    return( W/Ncontig/2 )
 }
W.rep <- replicate(Nreps, sapply(select.r,WAccordingToR,indivs=indivs) )
plot(select.r,rowMeans(W.rep),'l')
for (i in 1:20) {
  points(select.r,W.rep[,i])
}

#With different geographic distance restriction
select.dis<-seq(100,3100,200)
WAccordingToDis <- function(indivs,select.r,select.dis) {
  contig.div <- cumsum( rexp(1*Ncontig,1) )#sample from the dateset
  d.sub <- subset(d, gdist <= select.dis )
  pair.factor <- with( d.sub, paste( ifelse( id1<id2, id1, id2 ),
                                     ifelse( id1<id2, id2, id1 ), sep="-" ) )
  data.mapstart <- d.sub$mapstart
  data.mapend <- d.sub$mapend
  N.pairs <- OverLapsCountAll(data.mapstart,data.mapend,contig.div,select.r)
  # calculate w
  Z <- (N.pairs-mean(N.pairs))*Ncontig/sqrt(sum((N.pairs-mean(N.pairs))^2))
  W <- sum( Z[ (row(Z)+1)==col(Z) ] )
  return( W/Ncontig/2 )
}

W.rep <- replicate(Nreps, sapply(select.dis,WAccordingToDis,indivs=indivs,select.r=2) )
plot(select.dis,W.selectdis)
write.csv(W.rep,"WaccR.csv")
write.csv(W.selectdis,"WaccDis.csv")