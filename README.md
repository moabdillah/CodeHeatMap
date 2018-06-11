# CodeHeatMap

library(WGCNA);
library(flashClust)

install.packages("flashClust")

suppressMessages(library(WGCNA))
allowWGCNAThreads()
## Allowing multi-threading with up to 4 threads.
suppressMessages(library(cluster)) 

options(stringsAsFactors = FALSE)  
# Read in the female liver data set 
DataCF = read.csv("FPKM_27>0.1.csv",header=TRUE, sep = ";");
dim(DataCF)
# the first 8 columns contain gene information 
names(DataCF)[1:9]
# Remove gene information and transpose the expression data 
datExprCF = as.data.frame(t(DataCF[, -c(1:5)])) 
names(datExprCF) = DataCF$Gene 
rownames(datExprCF) = names(DataCF)[-c(1:5)]


# Now we read in the physiological trait data 
traitData = read.csv("ClinicalTraitsMethylCF_27.csv",header=TRUE, sep = ";");
dim(traitData)
names(traitData)
# use only a subset of the columns 
#allTraits = traitData[, c(2:31)]
allTraits = traitData[, -c(1,3,4,15,27,28,29,30,31)] 

names(allTraits)

# Order the rows of allTraits so that they match those of datExprCF: 
CodeCF = rownames(datExprCF) 
traitRows = match(CodeCF, allTraits$CodeCF) 
datTraits = allTraits[traitRows, -1] 
rownames(datTraits) = allTraits[traitRows, 1] 
# show that row names agree 
table(rownames(datTraits) == rownames(datExprCF))

# sample network based on squared Euclidean distance note that we transpose the data 
A = adjacency(t(datExprCF), type = "distance") 
# this calculates the whole network connectivity 
k = as.numeric(apply(A, 2, sum)) - 1 
# standardized connectivity 
Z.k = scale(k)
print(Z.k)

# Designate samples as outlying if their Z.k value is below the threshold 
thresholdZ.k = -5  # often -2.5 

# the color vector indicates outlyingness (red) 
outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")  
# calculate the cluster tree using flahsClust or hclust 

sampleTree = flashClust(as.dist(1 - A), method = "average") 

#sampleTree = hclust(as.dist(1 - A), method = "average") 
# Convert traits to a color representation: where red indicates high values 
traitColors = data.frame(numbers2colors(datTraits, signed = TRUE)) 
dimnames(traitColors)[[2]] = paste(names(datTraits), sep = "") 
datColors = data.frame(outlier = outlierColor, traitColors) 
# Plot the sample dendrogram and the colors underneath. 
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,      
                    main = "Sample dendrogram and trait heatmap")


# Remove outlying samples from expression and trait data 
remove.samples = Z.k < thresholdZ.k | is.na(Z.k) 
# the following 2 lines differ from what is written in the book 
datExprCF = datExprCF[!remove.samples, ] 
datTraits = datTraits[!remove.samples, ] 
# Recompute the sample network among the remaining samples 
A = adjacency(t(datExprCF), type = "distance") 
# Let's recompute the Z.k values of outlyingness 
k = as.numeric(apply(A, 2, sum)) - 1 
Z.k = scale(k)

# Choose a set of soft thresholding powers 
powers = c(1:20)  # in practice this should include powers up to 20. 
# choose power based on SFT criterion 
#sft = pickSoftThreshold(datExprCF, powerVector = powers)
sft=pickSoftThreshold(datExprCF,powerVector=powers, networkType = "signed")



# Plot the results: 
par(mfrow = c(1, 2)) 
# SFT index as a function of different powers 
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],     
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", 
     main = paste("Scale independence")) 
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "red") 
# this line corresponds to using an R^2 cut-off of h 
abline(h = 0.9, col = "red") 
# Mean connectivity as a function of different powers 
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)",      ylab = "Mean Connectivity", main = paste("Mean connectivity")) 
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")


mergingThresh = 0.25 
net = blockwiseModules(datExprCF, corType = "bicor", maxBlockSize = 16000,      
                       networkType = "signed", power = 12, minModuleSize = 30, mergeCutHeight = mergingThresh,      
                       numericLabels = TRUE, saveTOMs = TRUE, pamRespectsDendro = FALSE, 
                       saveTOMFileBase = "CFdataTOM") 
moduleLabelsAutomatic = net$colors 
table(net$colors)

# Convert labels to colors for plotting 
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)  
# A data frame with module eigengenes can be obtained as follows 
MEsAutomatic = net$MEs 

# this is the body FEV1 
FEV1 = as.data.frame(datTraits$FEV1_percent) 
names(FEV1) = "FEV1" 
# Next use this trait to define a gene significance variable 
GS.FEV1 = as.numeric(cor(datExprCF, FEV1, use = "p")) 
# This translates the numeric values into colors 
GS.FEV1Color = numbers2colors(GS.FEV1, signed = T) 
blocknumber = 1 
datColors = data.frame(moduleColorsAutomatic, GS.FEV1Color)[net$blockGenes[[blocknumber]],      ] 


# Plot the dendrogram and the module colors underneath 
plotDendroAndColors(net$dendrograms[[blocknumber]], colors = datColors, 
                    groupLabels = c("Module colors",      "GS.FEV1"), 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05)
####12.2.3 Blockwise module detection for large networks
bwnet = blockwiseModules(datExprCF, corType = "bicor", maxBlockSize = 16000,      
                         networkType = "signed", power = 12, minModuleSize = 30, 
                         mergeCutHeight = mergingThresh,      
                         numericLabels = TRUE, 
                         saveTOMs = TRUE, 
                         pamRespectsDendro = FALSE, 
                         saveTOMFileBase = "CFdataTOM-blockwise",      
                         verbose = 3)

# Relabel blockwise modules so that their labels match those from our previous analysis 
moduleLabelsBlockwise = matchLabels(bwnet$colors, moduleLabelsAutomatic) 
# Convert labels to colors for plotting 
moduleColorsBlockwise = labels2colors(moduleLabelsBlockwise)  
# measure agreement with single block automatic procedure 
mean(moduleLabelsBlockwise == moduleLabelsAutomatic)

blockNumber = 1
# Plot the dendrogram for the chosen block 
plotDendroAndColors(bwnet$dendrograms[[blockNumber]], 
                    moduleColorsBlockwise[bwnet$blockGenes[[blockNumber]]],     
                    "Module colors", 
                    main = paste("Dendrogram and module colors in block", blockNumber),      
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05)

######12.2.4 Manual, stepwise module detection
# We now calculate the weighted adjacency matrix, using the power 7: 
A = adjacency(datExprCF, power = 12, type='signed') 
# Digression: to define a signed network choose A = adjacency(datExprCF, power = 12, type='signed')

# define a dissimilarity based on the topological overlap 
dissTOM = TOMdist(A)

# hierarchical clustering 
geneTree = flashClust(as.dist(dissTOM), method = "average") 
# here we define the modules by cutting branches 
moduleLabelsManual1 = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",      
                                    deepSplit = 2, 
                                    pamRespectsDendro = F, 
                                    minClusterSize = 30)
# TOMplot: Network heatmap plot, all genes###########
plotTOM = dissTOM^12;
# # Set diagonal to NA for a nicer plot

diag(plotTOM) = NA;
# # Call the plot function

sizeGrWindow(9,9)
TOMplot(dissim=dissTOM^12,dendro=geneTree,colors=moduleColorsCF, main = "Network heatmap plot, all genes")

# Relabel the manual modules so that their labels match those from our previous analysis 

moduleLabelsManual2 = matchLabels(moduleLabelsManual1, moduleLabelsAutomatic) 
# Convert labels to colors for plotting 
moduleColorsManual2 = labels2colors(moduleLabelsManual2)
# Calculate eigengenes 
MEList = moduleEigengenes(datExprCF, colors = moduleColorsManual2) 
MEs = MEList$eigengenes

# Add the weight to existing module eigengenes 
MET = orderMEs(cbind(MEs)) 
# Plot the relationships among the eigengenes and the trait 
plotEigengeneNetworks(MET, "", marDendro = c(0, 4, 1, 2), 
                      marHeatmap = c(3,      4, 1, 2), 
                      cex.lab = 0.8, xLabelsAngle = 90)

# automatically merge highly correlated modules 
merge = mergeCloseModules(datExprCF, moduleColorsManual2, cutHeight = mergingThresh)


# resulting merged module colors 
moduleColorsManual3 = merge$colors 
# eigengenes of the newly merged modules: 
MEsManual = merge$newMEs 

# Show the effect of module merging by plotting the original and merged module colors below the tree 

datColors = data.frame(moduleColorsManual3, moduleColorsAutomatic, 
                       moduleColorsBlockwise,      
                       GS.FEV1Color) 
plotDendroAndColors(geneTree, colors = datColors, groupLabels = c("manual hybrid",      
                                                                  "single block", "2 block", "GS.FEV1"), 
                    dendroLabels = FALSE, 
                    hang = 0.03,      
                    addGuide = TRUE, 
                    guideHang = 0.05)



# check the agreement between manual and automatic module labels 
mean(moduleColorsManual3 == moduleColorsAutomatic)

###Relating modules to physiological traits Choose a module assignment 

moduleColorsCF = moduleColorsAutomatic 
# Define numbers of genes and samples 
nGenes = ncol(datExprCF) 
nSamples = nrow(datExprCF) 

# Recalculate MEs with color labels 
MEs0 = moduleEigengenes(datExprCF, moduleColorsCF)$eigengenes 
MEsFC = orderMEs(MEs0) 
modTraitCor = cor(MEsFC, datTraits, use = "p") 
modTraitP = corPvalueStudent(modTraitCor, nSamples)


# names (colors) of the modules
modNames = substring(names(MEsFC), 3)
geneModuleMembership = as.data.frame(cor(datExprCF, MEsFC, use = "p"));

###Modules MEsFC utilisés#########
MMPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
write.table(MMPvalue1, "~/Resultat-RDataCF/GenModule-p.MM1.txt", sep="\t")


# Since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table. We 
# color code each association by the correlation value: Will display correlations and their p-values 

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")",      sep = "") 
dim(textMatrix) = dim(modTraitCor) 

#par(mar = c(6, 8.5, 3, 3))
par(mar = c(6, 13, 3, 3))

# Display the correlation values within a heatmap plot 
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = names(MEsFC),      
               ySymbols = names(MEsFC), colorLabels = FALSE, colors = greenWhiteRed(50),      
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, 
               zlim = c(-1,          1), main = paste("Module-trait relationships"))


#######12.4.1 Connectivity-, TOM-, and MDS plots

# # WARNING: On some computers, this code can take a while to run (20 minutes??).  
# # I suggest you skip it. 
# # Set the diagonal of the TOM disscimilarity to NA 

diag(dissTOM) = NA 

# # # Transform dissTOM with a power to enhance visibility 
TOMplot(dissim=dissTOM^12,dendro=geneTree,colors=moduleColorsCF, main = "Network heatmap plot, all genes")

##### Another methode for TOMplot #####

# ##Multidimensional scaling (MDS) A la place Du Plot (Tree, dendrogramme; cluster) qui bloque Le cluster###
cmd1 = cmdscale(as.dist(dissTOM), 2) 
par(mfrow = c(1, 1)) 
plot(cmd1, col = moduleColorsCF, 
      main = "MDS plot", xlab = "Scaling Dimension 1",      ylab = "Scaling Dimension 2")


####Gene relationship to trait and important modules: Gene Significance and Module

#calculate the module membership values (aka. module eigengene based # connectivity kME): 
datKME = signedKME(datExprCF, MEsFC)

#Intramodular analysis: identifying genes with high GS and MM

#######################################
FEV1 = as.data.frame(datTraits$FEV1_percent) 
names(FEV1) = "FEV1" 
# Next use this trait to define a gene significance variable 
GS.FEV1 = as.numeric(cor(datExprCF, FEV1, use = "p")) 
# This translates the numeric values into colors 
GS.FEV1Color = numbers2colors(GS.FEV1, signed = T) 
blocknumber = 1 
datColors = data.frame(moduleColorsAutomatic, GS.FEV1Color)[net$blockGenes[[blocknumber]],      ] 

######For FEV1#############
colorOfColumn = substring(names(datKME), 4) 
par(mfrow = c(2, 2))
selectModules = c("lightcyan", "cyan", "salmon", "royalblue") 

### greenyellow
par(mfrow = c(2, length(selectModules)/2)) 
for (module in selectModules) {     
  column = match(module, colorOfColumn)     
  restModule = moduleColorsCF == module     
  verboseScatterplot(datKME[restModule, column], 
                     GS.FEV1[restModule], 
                     xlab = paste("Module Membership ", module, "module"), 
                     ylab = "GS.FEV1_percent", main = paste("kME.", module, "vs. GS"), 
                     col = module) }

########## For Diabetes ###########
#calculate the module membership values (aka. module eigengene based # connectivity kME): 
datKME = signedKME(datExprCF, MEsFC)
Diabetes = as.data.frame(datTraits$diabetes) 
names(Diabetes) = "Diabetes" 

# Next use this trait to define a gene significance variable 
GS.Diabetes = as.numeric(cor(datExprCF, Diabetes, use = "p")) 

# This translates the numeric values into colors 
GS.DiabetesColor = numbers2colors(GS.Diabetes, signed = T) 
blocknumber = 1 
datColors1 = data.frame(moduleColorsAutomatic, GS.DiabetesColor)[net$blockGenes[[blocknumber]],      ] 
colorOfColumn1 = substring(names(datKME), 4) 
par(mfrow = c(2, 2))
selectModules1 = c("red", "lightyellow", "lightgreen","salmon") 
par(mfrow = c(2, length(selectModules1)/2)) 
for (module in selectModules1) {     
  column = match(module, colorOfColumn)     
  restModule = moduleColorsCF == module     
  verboseScatterplot(datKME[restModule, column], 
                     GS.Diabetes[restModule], 
                     xlab = paste("Module Membership ", module, "module"), 
                     ylab = "GS.Diabetes", main = paste("kME.", module, "vs. GS"), 
                     col = module) }



########## For BMI ###########
#calculate the module membership values (aka. module eigengene based # connectivity kME): 
datKME = signedKME(datExprCF, MEsFC)
BMI = as.data.frame(datTraits$BMI) 
names(BMI) = "BMI" 

# Next use this trait to define a gene significance variable 
GS.BMI = as.numeric(cor(datExprCF, BMI, use = "p")) 

# This translates the numeric values into colors 
GS.BMIColor = numbers2colors(GS.BMI, signed = T) 
blocknumber = 1 
datColors2 = data.frame(moduleColorsAutomatic, GS.BMIColor)[net$blockGenes[[blocknumber]],      ] 
colorOfColumn2 = substring(names(datKME), 4) 
par(mfrow = c(2, 2))
selectModules1 = c("lightcyan", "cyan", "grey60","salmon") 
par(mfrow = c(2, length(selectModules1)/2)) 
for (module in selectModules1) {     
  column = match(module, colorOfColumn)     
  restModule = moduleColorsCF == module     
  verboseScatterplot(datKME[restModule, column], 
                     GS.BMI[restModule], 
                     xlab = paste("Module Membership ", module, "module"), 
                     ylab = "GS.BMI", main = paste("kME.", module, "vs. GS"), 
                     col = module) }
####liste des gènes de chaque module avec leurs KME respectifs


#We define a data frame containing he module membership (MM) values for each module. 
#In the past, we called the module membership values kME.

datKME1=signedKME(datExprCF, MEsFC, outputColumnName="MM.")
head(datKME)
write.table(datKME1, "~/Resultat-RDataCF/GenModule-KME.txt", sep="\t")

#################################################################
 ####################### Choose Cytoscape or VisANT #############
 ################################################################
          ######### Cytoscape plot and software #############
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExprCF, power = 12);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotationCF.csv");
# Select modules
# select modules 
modules = c("lightgreen") 

# Select module probes 
probes = names(datExprCF)
inModule=is.finite(match(moduleColorsCF,modules)) 
modProbes=probes[inModule] 
match1=match(modProbes,annot$Description) 
modGenes=annot$Gene[match1] 
# Select the corresponding Topological Overlap 
modTOM = TOM[inModule, inModule] 
dimnames(modTOM) = list(modProbes, modProbes) 
# Export the network into edge and node list files for Cytoscape 
cyt = exportNetworkToCytoscape(modTOM,   
                               edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),   
                               nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),   
                               weighted = TRUE, threshold = 0.05,nodeNames=modProbes,   
                               altNodeNames = modGenes, 
                               nodeAttr = moduleColorsCF[inModule])


################## WARNING: VisANT Does not Work weell : Use Cytoscape #######################
############# VisANT plot and software ##################
## Recalculate topological overlap 
TOM = TOMsimilarityFromExpr(datExprCF, power=12)
# #Read in the anotation fle
GeneAnnotationCF = read.csv("GeneAnnotationCF.csv",header=TRUE, sep = ";") 
# ###Module For FEV1######################
# #Select module
module = "greenyellow" 
# # Select module probes 
probes = names(datExprCF)
inModule = (moduleColorsCF==module)
modProbes = probes[inModule] 

#Select the corresponding Topological Overlap 
modTOM = TOM[inModule, inModule] 
dimnames(modTOM) = list(modProbes, modProbes)

# # Export the network into an edge list file VisANT can read

# vis = exportNetworkToVisANT(modTOM[selectHubs,selectHubs],   
                            file=paste("VisANTInput-", module, ".xml", sep="\t"),   
                            weighted=TRUE,
                            threshold = 0.10, 
                            probeToGene=   data.frame(GeneAnnotationCF$,Gene
                                                      GeneAnnotationCF$Description))
# # Because the module is rather large, 
# # we focus on the 30 top intramodular hub genes
# # intramodular connectivity 
nTopHubs = 30
kIN = softConnectivity(datExprCF[, modProbes]) 
selectHubs = (rank (-kIN) <= nTopHubs) 
vis30 = exportNetworkToVisANT(modTOM[selectHubs,selectHubs],   
                              file=paste("VisANTInput-", module, "-top30.xml", sep="\t"),   
                              weighted=TRUE,
                               threshold = 0.10, 
                             probeToGene=   data.frame(GeneAnnotationCF$Gene,
                                                         GeneAnnotationCF$Description))


##############################################################################
############################# YES WE CAN #####################################
##############################################################################

