library(TSRchitect)

setwd("/home/ssnyde11/scratch/tsrchitect_tssthresh5") #point to your project directory

#cp all of your processed BAM files in a single directory called 'BAMDIR'
BAMDIR  <- c("/home/ssnyde11/scratch/tsrchitect_tssthresh5/BAMDIR/") #provide a full path to BAMDIR

PdGFF3  <- c("/home/ssnyde11/scratch/genomes/pulex_genome_files/PA42_v4_files/PA42.4.0.gff") #provide a full path to the Dp annotation

#this is the only tricky part of the entire script
#change the sampleNames to sampleA-repN where A is the stage code and N is the replicate number
PdSTRIPE <- loadTSSobj(experimentTitle="DpDevel_TSS",
  inputDir=BAMDIR, n.cores=4, isPairedBAM=TRUE,
  sampleNames=c("sample1-rep1", "sample1-rep2","sample2-rep1",
  "sample2-rep2"), replicateIDs=c(1,1,2,2)) #datasets 1-2 and 3-4 are replicates

pdownstream  <- 150
pupstream    <- 150

useClustDist <- 30
TSSthreshold <- 5

load(file= "PdSTRIPE_complete.RData")

PdSTRIPE <- mergeSampleData(experimentName=PdSTRIPE, n.cores=1, tagCountThreshold=5)

PdSTRIPE <- determineTSR(experimentName=PdSTRIPE, n.cores=1, tssSetType="merged", tssSet="all", tagCountThreshold=TSSthreshold, clustDist=useClustDist, writeTable=TRUE, mixedorder=TRUE)


#do all samples in tsrSet
PdSTRIPE <- addAnnotationToTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=10, upstreamDist=pupstream, downstreamDist=pdownstream, feature="gene", featureColumn="ID", writeTable=TRUE)

writeTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=1, tsrLabel="TSR_A", mixedorder=TRUE, fileType="bed")
writeTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=2, tsrLabel="TSR_B", mixedorder=TRUE, fileType="bed")
writeTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=3, tsrLabel="TSR_C", mixedorder=TRUE, fileType="bed")
writeTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=4, tsrLabel="TSR_D", mixedorder=TRUE, fileType="bed")
writeTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=5, tsrLabel="TSR_E", mixedorder=TRUE, fileType="bed")
writeTSR(experimentName=PdSTRIPE, tsrSetType="merged", tsrSet=6, tsrLabel="TSR_F", mixedorder=TRUE, fileType="bed")
