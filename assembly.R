#this script servrs assembly the contgigs by assemblying pair of contigs which has the highest value of residual
#read the table that stores the number of ibd shared by each pair of country on each pair
# of contigs
#setwd("/Users/xzge/Desktop/codes_ibd")

# read in the data
A <- read.csv("contig_all_ibdshared.csv",header=FALSE)  # peter added a header

# collapse to fewer populations
divide<-function (country) {  
  ifelse (country %in% easterns, 1,
  ifelse (country %in% northerns, 2,
  ifelse (country %in% mideasterns, 3,
  ifelse (country %in% italy, 4,
  ifelse (country %in% iberia, 5,
  ifelse (country %in% france, 6,
  ifelse (country %in% UK, 7,
  ifelse (country %in% bene, 8,
  ifelse (country %in% swiss, 9,
  ifelse (country %in% german, 10,
  ifelse (country %in% swiss.french, 11,
  ifelse (country %in% swiss.german, 12, NA
          ) ) ) ) ) ) ) ) ) ) ) )
}

easterns <- c( 'Albania', 'Kosovo', 'Croatia', 'Bosnia', 'Montenegro',
               'Serbia', 'Yugoslavia', 'Slovenia', 'Macedonia', 'Greece', 'Bulgaria',
               'Romania', 'Poland', 'Austria', 'Hungary', 'Slovakia', 'Czech Republic', 'Russia', 'Ukraine' )
northerns <- c('Latvia', 'Finland', 'Sweden', 'Norway', 'Denmark')
mideasterns <- c('Cyprus','Turkey')
italy <- c('Italy')
iberia <- c('Portugal','Spain')
france <- c('France','Europe')
UK <- c('United Kingdom','England', 'Scotland', 'Ireland')
bene <- c('Belgium','Netherlands')
swiss <- c('Switzerland')
german <- c('Germany')
swiss.french<-c('Swiss French')
swiss.german<-c('Swiss German')

# recode countries  (note this is vectorized now)
country.divide1 <- divide(A[,1])
country.divide2 <- divide(A[,2])

#A stores the positive number of ibd shared by each pair of contig and each pair of country
#As we want to do the generalized linear regression, we have to make a larger table by adding all the 'zeros' into the table
#below ncontigs will be 442
ncontigs <- max(c(A[,3],A[,4]))
ibdtable<-data.frame(contig1=integer(),contig2=integer(),N=integer(),stringsAsFactors=FALSE)
div <- c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german")
id1 <- div[rep(1:12,each=12*choose(ncontigs,2))]
id2 <- div[rep(rep(1:12,each=choose(ncontigs,2)),12)]
ij <- t(combn(1:ncontigs,2))
contig1 <- rep(ij[,1],12*12)
contig2 <- rep(ij[,2],12*12)
N <- rep(0,choose(ncontigs,2)*12*12)
ibdtable <- data.frame(id1,id2,contig1,contig2,N)
N.n <- data.frame(country.divide1,country.divide2,contig1=A[,3],contig2=A[,4],N=A[,5])
M <- aggregate(N.n$N,by=list(N.n$country.divide1,N.n$country.divide2,N.n$contig1,N.n$contig2),FUN = 'sum')
col.names <- as.numeric((M[,1]-1)*(12*choose(ncontigs,2))+(M[,2]-1)*choose(ncontigs,2)+(2*ncontigs-M[,3])*(M[,3]-1)/2+M[,4]-M[,3])
ibdtable[col.names,5] <- M[,5]
colnames(ibdtable) <- c("id1","id2","contig1","contig2","N")
ibdtable$id1 <- factor(ibdtable$id1,levels=c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german"))
ibdtable$id2 <- factor(ibdtable$id2,levels=c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german"))
contig.names <- sort(unique(as.numeric(c(ibdtable$contig1,ibdtable$contig2))))
ibdtable$contig1 <- factor(ibdtable$contig1,levels=contig.names)
ibdtable$contig2 <- factor(ibdtable$contig2,levels=contig.names)
ibdtable$N<-as.integer(ibdtable$N)
write.csv(ibdtable,"ibdshared_table.csv",row.names=FALSE)

#calculate there different scores of residual with the generalized linear model. The three scores are:
#residual of GLM, the exact value minus fitted value of GLM, and normalized value of residual of GLM

# TOO BIG: model matrix is of size 10^10
# the.glm <- glm(N~contig1+contig2+id1*id2,family=poisson("log"),data=ibdtable)
# this works:
install.packages("MatrixModels")
library(MatrixModels)
the.glm <- glm4(N~contig1+contig2+id1*id2,family=poisson("log"),data=ibdtable,sparse=TRUE)

ibdtable$resid <- resid(the.glm)
ibdtable$resid2 <- N.resid$N-fitted(the.glm)
ibdtable$resid3 <- N.resid$resid2/sqrt(fitted(the.glm))
write.csv(ibdtable,"ibdshared_table_with_resids.csv",row.names=FALSE)

#here we use the second score alone which performs best in former test
tmp <- aggregate( ibdtable$resid2,by=list(ibdtable$contig1,ibdtable$contig2),FUN="mean")
colnames(tmp) <- c("contig1","contig2","resid")
tmp <- tmp[order(tmp$resid,decreasing = TRUE),]
write.table(tmp,"A.sorted.txt",row.names = FALSE)
#the function that assembly the contigs with a given number of chromosomes 'n' 
#and a given minimun score of the residual 'min.score'. 'tmp' is the data frame that stores
#the scores of residual of each pair of contigs. The output is a data frame which contains five columns
#'contig.id' is the id of each contig, 'left' is the id of the left neighbour, 'right' the right neighbour, 'chrom.id' is the id of the chromosome each contig belongs to
assembly <- function(n, min.score, tmp) {
   n.assembled <- 0
   n.all <- length(levels(tmp$contig1))
   contig.id <- 1:n.all
   left <- rep(0,n.all)
   right <- rep(0,n.all)
   chrom.id <- rep(0,n.all)
   chrom <- 0
   while (length(left[left==0])>n&tmp$resid[1]>min.score) {
     if (chrom.id[tmp$contig1[1]]==chrom.id[tmp$contig2[1]]&chrom.id[tmp$contig1[1]]!=0)  {
       tmp <- tmp[-1,]
       } else {
       if (chrom.id[tmp$contig1[1]]==0&chrom.id[tmp$contig2[1]]==0) {
         chrom <- chrom+1
         chrom.id[tmp$contig1[1]] <- chrom
         chrom.id[tmp$contig2[1]] <- chrom
         right[tmp$contig1[1]] <- tmp$contig2[1]
         left[tmp$contig2[1]] <- tmp$contig1[1]
       
         } else if (chrom.id[tmp$contig1[1]]!=0&chrom.id[tmp$contig2[1]]==0) {
         chrom.id[tmp$contig2[1]] <- chrom.id[tmp$contig1[1]]
         n.assembled <- n.assembled+1
         if (right[tmp$contig1[1]]==0) {
           left[tmp$contig2[1]] <- tmp$contig1[1]
           right[tmp$contig1[1]] <- tmp$contig2[1]
         } else {
           right[tmp$contig2[1]] <- tmp$contig1[1]
           left[tmp$contig1[1]] <- tmp$contig2[1]
         }
         tmp <- subset(tmp,contig1!=tmp$contig1[1]&contig2!=tmp$contig1[1])
       
         } else if (chrom.id[tmp$contig1[1]]==0&chrom.id[tmp$contig2[1]]!=0) {
         chrom.id[tmp$contig1[1]] <- chrom.id[tmp$contig2[1]]
         n.assembled <- n.assembled+1
         if (right[tmp$contig2[1]]==0) {
           left[tmp$contig1[1]] <- tmp$contig2[1]
           right[tmp$contig2[1]] <- tmp$contig1[1]
         } else {
           right[tmp$contig1[1]] <- tmp$contig2[1]
           left[tmp$contig2[1]] <- tmp$contig1[1]
         }
         tmp <- subset(tmp,contig1!=tmp$contig2[1]&contig2!=tmp$contig2[1])
       
         } else if (chrom.id[tmp$contig1[1]]!=0&chrom.id[tmp$contig2[1]]!=0) {
         n.assembled <- n.assembled+2
         if (right[tmp$contig2[1]]==0&left[tmp$contig1[1]]==0) {
           left[tmp$contig1[1]] <- tmp$contig2[1]
           right[tmp$contig2[1]] <- tmp$contig1[1]
         } else if (right[tmp$contig1[1]]==0&left[tmp$contig2[1]]==0){
           right[tmp$contig1[1]] <- tmp$contig2[1]
           left[tmp$contig2[1]] <- tmp$contig1[1]
         
           } else if (right[tmp$contig2[1]]==0&right[tmp$contig1[1]]==0) {
           
           for (i in 1:n.all) {
             if (chrom.id[i]==chrom.id[tmp$contig2[1]]) {
               a <- left[i]
               left[i] <- right[i]
               right[i] <- a
             }
           }
             right[tmp$contig1[1]] <- tmp$contig2[1]
             left[tmp$contig2[1]] <- tmp$contig1[1]
         
           } else if (left[tmp$contig1[1]]==0&left[tmp$contig2[1]]==0){
           
           for (i in 1:n.all) {
             if (chrom.id[i]==chrom.id[tmp$contig2[1]]) {
               a <- left[i]
               left[i] <- right[i]
               right[i] <- a
             }
           } 
            right[tmp$contig2[1]] <- tmp$contig1[1]
            left[tmp$contig1[1]] <- tmp$contig2[1]
         
           }
         b <- chrom.id[tmp$contig1[1]]
         for (i in 1:n.all) {
           if (chrom.id[i]==b) {
             chrom.id[i] <- chrom.id[tmp$contig2[1]]
           }
         }
         tmp <- subset(tmp,contig1!=tmp$contig2[1]&contig2!=tmp$contig2[1]&contig1!=tmp$contig1[1]&contig2!=tmp$contig1[1])
         }
       }
   }
   sort.table<-data.frame(contig.id,left,right,chrom.id)
} 
result <- assembly(50,0.05,tmp)
write.table(result,"Result.txt",row.names = FALSE)
result.to.list <- function( assembly.result ) {
  chrom.names<-unique(assembly.result$chrom.id)
  chrom.list<-list()
  for (i in 1:length(chrom.names)) {
    chrom.indi<-vector()
    contig.now <- subset(assembly.result,left==0&chrom.id==chrom.names[i])$contig.id[1]
    contig.next <- subset(assembly.result,contig.id==contig.now)$right[1]
    chrom.indi[1]<-contig.now
    while (contig.next!=0) {
      contig.now <- contig.next
      contig.next <- subset(assembly.result,contig.id==contig.now)$right[1]
      chrom.indi[length(chrom.indi)+1]<-contig.now
    }
    if (chrom.names[i]!=0) chrom.list[[length(chrom.list)+1]] <- chrom.indi
  }
  for (i in 1:nrow(assembly.result)) {
    if (assembly.result$left[i]==0&assembly.result$right[i]==0) chrom.list[[length(chrom.list)+1]] <- result$contig.id[i]
  }
  chrom.list
}
