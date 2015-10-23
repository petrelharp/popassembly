#this script servrs assembly the contgigs by assemblying pair of contigs which has the highest value of residue
#read the table that stores the number of ibd shared by each pair of country on each pair
# of contigs
setwd("/Users/xzge/百度云同步盘/dy/chr1ibd")
A <- read.table("chr1table.txt")
A$id1 <- as.factor(A$id1,levels=c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german"))
A$id2 <- as.factor(A$id2,levels=c("easterns","northerns","mideasterns","italy","iberia","france","UK","bene","swiss","german","swiss.french","swiss.german"))
contig.names <- sort(unique(as.numeric(c(A$contig1,A$contig2))))
A$contig1 <- factor(A$contig1,levels=contig.names)
A$contig2 <- factor(A$contig2,levels=contig.names)
A$N<-as.integer(A$N)
the.glm <- glm(A$N~A$contig1+A$contig2+A$id1*A$id2,family=poisson("log"))
resid(the.glm) -> resid
cbind(A,resid) -> N.resid
N.resid$resid2 <- N.resid$N-fitted(the.glm)
N.resid$resid3 <- N.resid$resid2/sqrt(fitted(the.glm))
tmp <- aggregate( N.resid$resid2,by=list(N.resid$contig1,N.resid$contig2),FUN="mean")
colnames(tmp) <- c("contig1","contig2","resid")
tmp <- tmp[order(tmp$resid,decreasing = TRUE),]

assembly <- function(n, tmp) {
   n.assembled <- 0
   n.all <- length(levels(tmp$contig1))
   contig.id <- 1:n.all
   left <- rep(0,n.all)
   right <- rep(0,n.all)
   chrom.id <- rep(0,n.all)
   chrom <- 0
   while (n.assembled<n.all-2*n) {
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
         for (i in 1:n.all) {
           b <- chrom.id[tmp$contig1[1]]
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
