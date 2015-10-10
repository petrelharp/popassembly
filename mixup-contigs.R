# zcat ../thedata/POPRES_Genotypes_QC2_v2_TXT.tped.gz | awk '$1 == 1' > chr1.tped
# note that the column IDs is in ../thedata/POPRES_Genotypes_QC2_v2_TXT.tfam.gz

all.mixup.data <- data.frame()

mean.length <- 5e6  # in bp

# 23 is X, 25 is Y
for (chrom in c(1:23,25)) {
    new.seed <- round(1e6*runif(1))
    set.seed(new.seed)
    outdir <- paste("chrom",chrom,new.seed,sep="_")
    dir.create(outdir)

    # pos <- scan( pipe("awk '{print $4}' <chr1.tped") )
    this.chrom <- paste("zcat ../thedata/POPRES_Genotypes_QC2_v2_TXT.tped.gz | awk '$1 == ",chrom,"'")
    pos <- scan( pipe(paste(this.chrom, " | awk '{print $4}'")) )
    breakpoints <- floor(c(0,cumsum(rexp(2*max(pos)/mean.length,1/mean.length))))
    npos <- cut(pos,breaks=breakpoints)
    breakpoints <- breakpoints[c(TRUE,table(npos)>0)]
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
    cmdlines <- with( cmds, paste( this.chrom, "| tail -n +", skip.lines+1, " | head -n ", n.lines, 
                    "| awk '{$4 = ", ifelse(rev==1, paste(maxpoints[contig.order],"- $4"), paste("$4 -",minpoints[contig.order])), 
                    "; $1 = ", new.seed,
                    "; print}'",
                    ifelse(rev==1,"| tac",""), " | gzip -c >", file.path(outdir,paste("chrom_",new.seed,"_contig_", contig, ".tped.gz", sep= '')), sep='') )
    # cat(paste(cmdlines,collapse="\n"))

    for (x in cmdlines) { system(x) }

    mixup.data <- data.frame(
                chrom=chrom,
                chrom.code=new.seed,
                    endbreak=breakpoints[-1],
                    min=minpoints,
                    max=maxpoints,
                    skip=skip.lines,
                    nlines=n.lines,
                    orig_order=contig.order,
                    contig=seq_along(contig.order)
                    )
    rownames(mixup.data) <- NULL
    write.table( mixup.data, file=file.path(outdir,paste("chrom_",new.seed,"_contig_order.tsv",sep='')), row.names=FALSE, sep="\t" )
    all.mixup.data <- rbind( all.mixup.data, mixup.data )
}

write.table( all.mixup.data, file=paste("all-mixup-data-",format(Sys.Date(),"%d-%m-%Y"),".tsv",sep=''), sep='\t' )

