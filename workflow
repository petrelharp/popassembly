1. Made `mixedup_1Mb/` with `mixup-contigs.R` and mean contig size set at 1e6
2. ran `parallel ../tped-to-vcf.sh ::: */*.tped.gz` to convert everything to vcf
3. ran `qsub run-beagle-on-everything.pbs` to run `beagle.sh` on everything

