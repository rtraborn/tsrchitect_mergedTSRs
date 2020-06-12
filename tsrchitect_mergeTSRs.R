library(TSRchitect)

setwd("/home/ssnyde11/scratch/tsrchitect_tssthresh5")

BAMDIR  <- c("/home/ssnyde11/scratch/tsrchitect_tssthresh5/BAMDIR/")
PdGFF3  <- c("/home/ssnyde11/scratch/genomes/pulex_genome_files/PA42_v4_files/PA42.4.0.gff")

pdownstream  <- 150
pupstream    <- 150

useClustDist <- 30
TSSthreshold <- 5

load(file= "PdSTRIPE_complete.RData")

PdSTRIPE <- mergeSampleData(experimentName=PdSTRIPE, n.cores=1, tagCountThreshold=5)

PdSTRIPE <- determineTSR(experimentName=PdSTRIPE, n.cores=1, tssSetType="merged", tssSet="all", tagCountThreshold=TSSthreshold, clustDist=useClustDist, writeTable=TRUE, mixedorder=TRUE)


#do all samples in tsrSet
PdSTRIPE <- addAnnotationToTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=10, upstreamDist=pupstream, downstreamDist=pdownstream, feature="gene", featureColumn="ID", writeTable=TRUE)

#writeTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=1, tsrLabel="TSR_", mixedorder=TRUE, fileType="bed")
