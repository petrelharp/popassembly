
# these are the ids we have in the data
# from: zcat 03.POPRES_chr22.beagle.phased.gz | head -n 1 | cut -f 3- -d ' '  | tr ' ' '\n' > ../xinzhou/ids
ids <- scan("ids",skip=2)
# this is phs000145.v4.pht000659.v2.p2.c1.POPRES_v1_v2_Subject_Phenotypes.GRU.txt.gz
# BUT WITH EXTRA TAB ADDED AT END OF HEADER LINE
euro <- read.table("sample-info-2.txt",skip=10,header=TRUE,sep='\t')
euro <- subset(euro, euro$SUBJID %in% ids  )


# List of countries, with 0 or 1 if in europe or not
in.Europe <- read.table("dbgap_In_Europe.out",as.is=T)
europe.countries<-in.Europe[in.Europe$V1==1,2]

# Combine certain countries/locations together.
combine.country<-list()
# a finer grouping
combine.country <- c( combine.country, list(c("Russia","USSR")) )
combine.country <- c( combine.country, list(c("Netherlands","Holland")) )


# Make new COUNTRY_SELF that agrees with grandparents
# record original information
euro$ORIG_COUNTRYSELF<-euro$COUNTRY_SELF
gfolx <- as.matrix( euro[,c("COUNTRY_MGM","COUNTRY_MGF","COUNTRY_PGM","COUNTRY_PGF")] )
##For individuals with no granf
par.country <- as.matrix( euro[,c("COUNTRY_FATHER","COUNTRY_MOTHER")] )

euro$MIXEDGFOLX <- apply(gfolx, 1, function (x) length(unique(x))>1 )
num.grandpar<-apply(gfolx,1,function(grand){4-sum(grand=="")}) 


# Is the individual "mixed"?
euro$MIXED <- euro$GROUPING_PCA_LABEL1 == "Mix" | euro$COUNTRY_SELF=="Europe" | euro$MIXEDGFOLX
euro$COUNTRY_SELF <- NA
euro$COUNTRY_GFOLX <- NA
euro$COUNTRY_SELF[!euro$MIXED] <- euro$COUNTRY_GFOLX[!euro$MIXED] <- gfolx[!euro$MIXED,1]
## If no grandpar info. available use country self
euro$COUNTRY_SELF[num.grandpar==0]<-as.character(euro$ORIG_COUNTRYSELF[num.grandpar==0]) 

# combine countries
for (i in 1:length(combine.country)) {
    euro$COUNTRY_SELF[euro$COUNTRY_SELF %in% combine.country[[i]]]<-combine.country[[i]][1]
}

# split countries by language
# Switzerland
euro$COUNTRY_SELF[ euro$COUNTRY_SELF=="Switzerland" & euro$PRIMARY_LANGUAGE=="French" ] <- "Swiss French"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF=="Switzerland" & euro$PRIMARY_LANGUAGE=="German" ] <- "Swiss German"
# Balkans:
# with( subset(euro,COUNTRY_SELF %in% c("Yugoslavia","Albania","Serbia", "Macedonia", "Montenegro", "Croatia", "Kosovo", "Bosnia" )), table( COUNTRY_SELF, droplevels(PRIMARY_LANGUAGE) ) )
# COUNTRY_SELF    Albanian Bosnian Croatian French Hungarian Kosovan Macedonian Romanian Serbian Serbo-Croatian Yugoslavian
#   Albania     0        3       0        0      0         0       0          0        0       0              0           0
#   Bosnia      0        0       4        0      0         0       0          0        0       1              4           0
#   Croatia     0        0       0        7      0         0       0          0        0       0              1           0
#   Kosovo      1       11       0        0      0         0       2          0        0       0              2           1
#   Macedonia   0        0       0        0      0         0       0          4        0       0              0           0
#   Montenegro  0        0       0        0      0         0       0          0        0       1              0           0
#   Serbia      0        1       0        0      0         1       0          0        0       6              4           0
#   Yugoslavia  1        5       0        1      1         0       1          0        1       2              3           4
balkans <- c( "Albania", "Bosnia", "Croatia", "Kosovo", "Macedonia", "Montenegro", "Serbia", "Yugoslavia" )
balkan.langs <- c( "Albanian", "Bosnian", "Croatian", "Kosovan", "Macedonian", "Serbian", "Serbo-Croatian" )
# notes: "kosovan" probably = "albanian"
#  "serbo-croatian" includes serbian, croatian, bosnian and probably = "yugoslavian"
#   macedonian "forms a continuum" of south slavic languages with bulgarian and serbo-croatian
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia","Serbia") & euro$PRIMARY_LANGUAGE=="Albanian" ] <- "Albania"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia") & euro$PRIMARY_LANGUAGE=="Croatian" ] <- "Croatia"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia") & euro$PRIMARY_LANGUAGE=="Kosovan" ] <- "Kosovo"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia") & euro$PRIMARY_LANGUAGE=="Serbian" ] <- "Serbia"

write.table(euro,file="sample-info.csv",sep=",",row.names=FALSE)
