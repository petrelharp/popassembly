1. Made `mixedup_1Mb/` with `mixup-contigs.R` and mean contig size set at 1e6
2. ran `parallel ../tped-to-vcf.sh ::: */*.tped.gz` to convert everything to vcf
3. ran `qsub run-beagle-on-everything.pbs` to run `beagle.sh` on everything
4. ran `../check-runs.sh >finished.runs 2>unfinished.runs` to find which contigs need to be finished
5. finished these with `parallel -j 6 ../beagle.sh ::: $(cat unfinished.runs | awk ' {print $1".vcf.gz"} ')`
6. ran `../add-sort-info.sh` to sort them and add an addition sort key column at end

