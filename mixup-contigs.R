# zcat ../thedata/POPRES_Genotypes_QC2_v2_TXT.tped.gz | awk '$1 == 1' > chr1.tped

pos <- scan( pipe("awk '{print $4}' <chr1.tped") )
breakpoints <- c(0,cumsum(rexp(60,40/max(pos))))
npos <- cut(pos,breaks=breakpoints)
breakpoints <- breakpoints[c(TRUE,table(npos)>200)]
npos <- cut(pos,breaks=breakpoints)

minpoints <- tapply( pos, npos, min )
maxpoints <- tapply( pos, npos, max )
minpos <- pos - minpoints[npos]
maxpos <- maxpoints[npos] - pos

contig.order <- sample.int(length(breakpoints)-1)

skip.lines <- as.numeric(c(0,cumsum(table(npos)))[contig.order])
n.lines <- as.numeric(table(npos)[contig.order])


cmds <- data.frame( contig=seq_along(contig.order), skip=skip.lines, nl=n.lines, rev=rbinom(length(contig.order),1,1/2) )
# with( cmds, paste( "cat chr1.tped | tail -n +", skip.lines+1, " | head -n ", n.lines, ifelse(rev==1,"| tac",""), " | gzip -c >chr1_contig_", contig, ".tped.gz", sep= '' ) )
cmdlines <- with( cmds, paste( "cat chr1.tped | tail -n +", skip.lines+1, " | head -n ", n.lines, 
                "| awk '{$4 = ", ifelse(rev==1, paste(maxpoints[contig.order],"- $4"), paste("$4 -",minpoints[contig.order])), "; print}'",
                ifelse(rev==1,"| tac",""), " | gzip -c >chr1_contig_", contig, ".tped.gz", sep= '' ) )

cat(paste(cmdlines,collapse="\n"))

write.table( data.frame(
                endbreak=breakpoints[-1],
                min=minpoints,
                max=maxpoints,
                skip=skip.lines,
                nlines=n.lines,
                order=contig.order
            ), file="chr1_contig_order.tsv", sep="\t" )
