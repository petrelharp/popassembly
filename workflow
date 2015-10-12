1. Made `mixedup_1Mb/` with `mixup-contigs.R` and mean contig size set at 1e6
2. ran `parallel ../tped-to-vcf.sh ::: */*.tped.gz` to convert everything to vcf
3. ran `qsub run-beagle-on-everything.pbs` to run `beagle.sh` on everything


Notes:

* Some contigs (three?) just sit there with java using no CPU or memory, e.g.  `/panfs/cmb-panasas2/pralph/popassembly/mixedup_1Mb/chrom_8_604778/chrom_604778_contig_6.beagle.log`.
