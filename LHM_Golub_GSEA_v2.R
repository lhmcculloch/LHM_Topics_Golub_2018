#LHM
#9/26/18
#GSEA for Golub et al. 1999 Training Data

library(golubEsets)
library(fgsea)
library(hu6800.db)
library(EnrichmentBrowser)
library(reactome.db)
library(globaltest)
library(GO.db)
library(GSEABase)

data(Golub_Train) #Did this for the train data, but could do it for the merged dataset with data(Golub_Merge)

#Getting rid of all of the AFFX controls, then the Alu control
GTrainNoCtrl <- Golub_Train[!grepl("AFFX", rownames(Golub_Train)),]
GTrainNoCtrl <- GTrainNoCtrl[!grepl("hum_alu_at", rownames(GTrainNoCtrl)),]

#Pulling out the ALL and AML samples:
allSamples <- c()
amlSamples <- c()
for (i in 1:length(GTrainNoCtrl$Samples)) {
  if (GTrainNoCtrl$ALL.AML[i] == "ALL") {
    allSamples <- c(allSamples, GTrainNoCtrl$Samples[i])
  }
  else if (GTrainNoCtrl$ALL.AML[i] == "AML") {
    amlSamples <- c(amlSamples, GTrainNoCtrl$Samples[i])
  }
}

#Pulling out expression values by tumor type:
gTrainExp <- exprs(GTrainNoCtrl)
allExp <- gTrainExp[,allSamples]
amlExp <- gTrainExp[,amlSamples]

#Expression means across all samples, to be used for computing a ratio:
allExpMeans <- rowMeans(allExp)
amlExpMeans <- rowMeans(amlExp)

#Means ratio for ALL/AML. I used this to come up with the "ranks" list that fsgea called for,
#but I'm not sure how it compares to how they expect rank lists to be generated.
allOverAmlRatio <- allExpMeans/amlExpMeans

#Sorting the list by expression level. I saw some contradictory info about whether it should be ascending or
#descending, or whether it mattered. The ranked list seemed to be ascending (negative to positive), so I used
#that approach first. Testing both ways, the results seemed to be about the same.
allOverAmlRatioAsc <- sort(allOverAmlRatio)

#Converting the gene IDs from the hu6800 system used in early AFFY data to the Entrez IDs used by fgsea
hu6800Entrez <- as.list(hu6800ENTREZID)
mapped_probes <- mappedkeys(hu6800Entrez)
probesList <- as.list(hu6800Entrez[mapped_probes])

#Pulled out the Entrez IDs and the corresponding gene expression values and paired them. Not the most efficient way
#of doing this, but it seemed to work.
probesRepNames <- c()
probesRepData <- c()

for (i in 1:length(allOverAmlRatioAsc)) {
  for (j in 1:length(probesList)) {
    if (names(probesList[j]) == names(allOverAmlRatioAsc[i])) {
      probesRepNames <- c(probesRepNames, probesList[j])
      probesRepData <- c(probesRepData, allOverAmlRatioAsc[i])
    }
  }
}

#Checking that the Entrez IDs and expression data vectors are 1:1 in length before using them to make a new named
#num vector. Also manually checked the conversions for a few probes here (not shown).
length(probesRepNames) == length(probesRepData)

#Named, ranked list of gene expression mean ratios with Entrez IDs, in a format that fgsea appears to like:
names(probesRepData) <- probesRepNames

#Obtaining a gene set to look at. This is a spot where I wasn't sure what was best to look at, so I tried this
#for now. Could try other gene sets instead.
gs <- get.go.genesets(org="hsa", mode="biomart")

#Getting rid of the sample with a ratio of infinity to avoid errors in fgsea (it's the last gene in the ascending
#6014 gene set)
probesRepDataNoInf <- probesRepData[1:6013]

#Ran fgsea. Used the conditions recommended by the package tutorial as a starting point.
fgseaRes <- fgsea(pathways = gs, 
                  stats = probesRepDataNoInf,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)

#Looking at the results:
head(fgseaRes[order(pval), ])

#Looking at a couple of measures of significance:
sum(fgseaRes[, padj<0.05]) #Upon further examination of the data, this adjusted p value reached its minimum at ~0.8,
#so nothing came close to reaching this cutoff.
sum(fgseaRes[, pval<0.05]) #Using raw values, I could get 29-31 pathways (seemed to vary depending on the run)
#meeting this threshold

#Graphing the top 10 upregulated and downregulated genes in ALL vs. AML (based on fgsea tutorial):
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gs[topPathways], probesRepDataNoInf, fgseaRes, gseaParam = 0.5)

#Trying something similar with the reactomes data. There were a couple of warning messages this time, which were not an issue with the other dataset.
pathways <- reactomePathways(names(probesRepDataNoInf))
fgseaRes <- fgsea(pathways, probesRepDataNoInf, nperm=10000, maxSize=500)
head(fgseaRes[order(pval), ])

#Looking at the results
sum(fgseaRes[, pval<0.05]) #107 genes met this threshold (but none met a similar threshold for padj)

#Given the lack of interesting gene findings by this method, and the fact that these results seemed surprising
#given the expected differences in the dataset, I decided to try another method that looks at all samples across
#the dataset (rather than a list of ranked means) to see how the results compared. The package globaltest appears
#to recommend additional normalization on their data (they even show some examples for a couple of Golub genes), but
#for consistency with the rest of our efforts, I'm leaving our starting data as-is.

eg <- as.list(hu6800ENTREZID)
golubGtGO <- gtGO(ALL.AML, Golub_Train, probe2entrez = eg, annotation = "org.Hs.eg.db")
head(golubGtGO) #This gives more significant p-values for the data. A lot of the terms are fairly vague, but it does
#point out a lot of differentially regulated pathways.

#To see if I can get some additional gene information, trying the same method with the Broad gene sets.
#The gene set data includes information for genes not present in my dataset, so some of the results are NA:
broad <- getBroadSets("msigdb_v6.2.xml") #include path to local directory
broadRes <- gtBroad(ALL.AML, Golub_Train, collection = broad)
broadRes[1:25] #There are a lot of significant pathways for this data, though some are fairly vague (e.g.,
#"MODULE_53", which appears at the top of the list - it is a set of cell-line expressed, cancer-associated genes,
#but the meaning isn't readily understandable from the name).


