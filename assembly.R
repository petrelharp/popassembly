#this script servrs assembly the contgigs by assemblying pair of contigs which has the highest value of residue
#read the table that stores the number of ibd shared by each pair of country on each pair
# of contigs
#setwd("/Users/xzge/Desktop/codes_ibd")
A <- read.csv("contig_all_ibdshared.csv",header=FALSE)
divide<-function (country) {  
  if (country %in% easterns) 1
  else if (country %in% northerns) 2
  else if (country %in% mideasterns) 3
  else if (country %in% italy) 4
  else if (country %in% iberia) 5
  else if (country %in% france) 6
  else if (country %in% UK) 7
  else if (country %in% bene) 8
  else if (country %in% swiss) 9
  else if (country %in% german) 10
  else if (country %in% swiss.french) 11
  else if (country %in% swiss.german) 12
  else NA
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
country.divide1<-sapply(A[,1],divide)
country.divide2<-sapply(A[,2],divide)
#A stores the positive number of ibd shared by each pair of contig and each pair of country
#As we want to do the generalized linear regression, we have to make a larger table by adding all the 'zeros' into the table
#the integer 432 below means the number of the contigs
table<-data.frame(contig1=integer(),contig2=integer(),N=integer(),stringsAsFactors=FALSE)
div <- c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german")
id1 <- div[rep(1:12,each=12*choose(442,2))]
id2 <- div[rep(rep(1:12,each=choose(442,2)),12)]
ij <- t(combn(1:442,2))
contig1 <- rep(ij[,1],12*12)
contig2 <- rep(ij[,2],12*12)
N <- rep(0,choose(442,2)*12*12)
table <- data.frame(id1,id2,contig1,contig2,N)
N.n <- data.frame(country.divide1,country.divide2,contig1=A[,3],contig2=A[,4],N=A[,5])
M <- aggregate(N.n$N,by=list(N.n$country.divide1,N.n$country.divide2,N.n$contig1,N.n$contig2),FUN = 'sum')
col.names <- as.numeric((M[,1]-1)*(12*choose(442,2))+(M[,2]-1)*choose(442,2)+(441+443-M[,3])*(M[,3]-1)/2+M[,4]-M[,3])
table[col.names,5] <- M[,5]
colnames(table) <- c("id1","id2","contig1","contig2","N")
table$id1 <- factor(table$id1,levels=c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german"))
table$id2 <- factor(table$id2,levels=c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german"))
contig.names <- sort(unique(as.numeric(c(table$contig1,table$contig2))))
table$contig1 <- factor(table$contig1,levels=contig.names)
table$contig2 <- factor(table$contig2,levels=contig.names)
table$N<-as.integer(table$N)
write.csv(table,"ibdshared_table.csv")
#calculate there different scores of residue with the generalized linear model. The three scores are:
#residue of GLM, the exact value minus fitted value of GLM, and normalized value of residue of GLM
the.glm <- glm(table$N~table$contig1+table$contig2+table$id1*table$id2,family=poisson("log"))
resid(the.glm) -> resid
cbind(table,resid) -> N.resid
N.resid$resid2 <- N.resid$N-fitted(the.glm)
N.resid$resid3 <- N.resid$resid2/sqrt(fitted(the.glm))
#here we use the second score alone which performs best in former test
tmp <- aggregate( N.resid$resid2,by=list(N.resid$contig1,N.resid$contig2),FUN="mean")
colnames(tmp) <- c("contig1","contig2","resid")
tmp <- tmp[order(tmp$resid,decreasing = TRUE),]
write.table(tmp,"A.sorted.txt",row.names = FALSE)
#the function that assembly the contigs with a given number of chromosomes 'n' 
#and a given minimun score of the residue 'min.score'. 'tmp' is the data frame that stores
#the scores of residue of each pair of contigs. The output is a data frame which contains five columns
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
result <- assembly(800,0.02,tmp)
write.table(result,"Result.txt",row.names = FALSE)