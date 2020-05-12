
#BSgenome.Mmusculus.UCSC.mm10

library(SummarizedExperiment)
library(MultiAssayExperiment)
library(CAGEr)

inputFiles = c(
"CAGE_E1_ACC_STAR/CAGE_E1_ACC.bam",
"CAGE_E2_CAC_STAR/CAGE_E2_CAC.bam",
"CAGE_E3_AGT_STAR/CAGE_E3_AGT.bam",
"CAGE_A1_ACG_STAR/CAGE_A1_ACG.bam",
"CAGE_A2_GCT_STAR/CAGE_A2_GCT.bam",
"CAGE_C1_GCG_STAR/CAGE_C1_GCG.bam",
"CAGE_C2_ATG_STAR/CAGE_C2_ATG.bam",
"CAGE_C3_TAC_STAR/CAGE_C3_TAC.bam" )

ce <- CAGEexp( genomeName = "BSgenome.Mmusculus.UCSC.mm10", 
               inputFiles = inputFiles, 
			   inputFilesType = "bam",
			   sampleLabels   =  c("E1", "E2", "E3", "A1", "A2", "C1", "C2", "C3")
	         )
getCTSS(ce, correctSystematicG = TRUE, removeFirstG = FALSE )
#
library(rtracklayer)

mm10_vM22_gene <- import.gff("~/Project/mm10/Annotation/GENCODE_vM22/gencode.vM22.annotation.gtf", format="gtf")

annotateCTSS(ce, mm10_vM22_gene) 
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
CTSScoordinatesGR(ce)
plotAnnot(ce, "counts")

### this step take long, ~30 min
#corr.m <- plotCorrelation2( ce, samples = "all", tagCountThreshold = 10, applyThresholdBoth = FALSE, method = "pearson")

mergeSamples(ce, mergeIndex = c(1,1,1,2,2,3,3,3),  mergedSampleLabels = c("E14", "PWWP2ABKO_A1", "PWWP2ABKO_C7"))
annotateCTSS(ce, mm10_vM22_gene)
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
plotAnnot(ce, "counts")

print("######################")
librarySizes(ce)
#[1] 73131612 35027566 55407033
plotReverseCumulatives(ce, fitInRange = c(5, 10000), onePlot = TRUE)
normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 10000), alpha = 1.05, T = 10^7)

####### clustering
clusterCTSS(object = ce, threshold=1, thresholdIsTpm = TRUE, nrPassThreshold=1, method = "distclu", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)
#clusterCTSS(object = ce, threshold=1, thresholdIsTpm = TRUE, nrPassThreshold=1, method = "distclu", maxDist = 50, removeSingletons = TRUE, keepSingletonsAbove = 10)

#head(tagClusters(ce, sample = "E14"))
#head(tagClusters(ce, sample = "PWWP2ABKO_A1"))
#head(tagClusters(ce, sample = "PWWP2ABKO_C7"))

##### 
cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = T)
quantilePositions(ce, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
##### Tag clusters and their interquantile width can be retrieved
tagClustersGR( ce, "E14", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

###
plotInterquantileWidth(ce, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
aggregateTagClusters(ce, tpmThreshold = 2, qLow = 0.1, qUp = 0.9, maxDist = 100)

print("#############################")

ce$outOfClusters / ce$librarySizes * 100
############ consensus Cluster
consensusClustersGR(ce)
annotateConsensusClusters(ce, mm10_vM22_gene)
cumulativeCTSSdistribution(ce, clusters = "consensusClusters", useMulticore = TRUE)
quantilePositions(ce, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = TRUE)

annotateCTSS(ce, mm10_vM22_gene)
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]

write.table(file="E14_CAGE_vM22.xls", consensusClustersGR(ce, sample = "E14", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9),sep="\t", quote=FALSE, row.names=F, col.names=T)
write.table(file="A1_CAGE_vM22.xls",  consensusClustersGR(ce, sample = "PWWP2ABKO_A1", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9), sep="\t", quote=FALSE, row.names=F, col.names=T)
write.table(file="C7_CAGE_vM22.xls",  consensusClustersGR(ce, sample = "PWWP2ABKO_C7", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9), sep="\t", quote=FALSE, row.names=F, col.names=T)

print("test3")
### output 
#exportCTSStoBedGraph(ce, values = "normalized", format = "bedGraph", oneFile = TRUE)
#exportCTSStoBedGraph(ce, values = "normalized", format = "BigWig")
#exportToBed(object = ce, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = TRUE)

#### GeneExpression profiling
cumulativeCTSSdistribution(ce, clusters = "consensusClusters")

