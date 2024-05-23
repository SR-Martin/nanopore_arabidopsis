# run methylkit package

# ==== install methylkit ====
#install.packages('devtools')
#install.packages('seqinr')
# 
#install_github("al2na/methylKit", build_vignettes=FALSE, 
#               repos=BiocManager::repositories(),
#                dependencies=TRUE)
#library(devtools)
library(methylKit)
library(seqinr)
library(stringr)


meth_tsv <- file.path('../data/methylation/ATRun1/ATRun1_CpG_5mC_counts_all.methylCall')
BS_tsv <- file.path('../data/methylation/BS/BS_CpG_frac.txt')

# TO ADD 2 TO POSITION.
nano=methRead(location = meth_tsv,sample.id="Nanopore ATRun1",assembly="mm9",header=TRUE, context="CpG", resolution="base", mincov = 10,
         pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=8 )
)

BS=methRead(location = BS_tsv,sample.id="Bisulfite",assembly="mm9",header=TRUE, context="CpG", resolution="base", mincov = 8,
         pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=8 )
)

getMethylationStats(nano,plot=T,both.strands=F)
getCoverageStats(nano,plot=TRUE,both.strands=FALSE)

getMethylationStats(BS,plot=T,both.strands=F)
getCoverageStats(BS,plot=TRUE,both.strands=FALSE)

meth_list <- as(list(BS,nano),"methylRawList")
meth_merge <- methylKit::unite(meth_list, destrand=T)
head(meth_merge)
getCorrelation(meth_merge, plot=TRUE)

