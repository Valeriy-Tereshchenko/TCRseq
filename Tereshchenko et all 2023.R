# TCR sequencing in mouse models of solid organ transplantation unveils the features of allorecognition and rejection 


# install and load

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

BiocManager::install("NOISeq")

BiocManager::install("baySeq")

install.packages('immunarch')
install.packages('assertthat')
install.packages('eulerr')
install.packages("readr")
install.packages("igraph")
install.packages("rgl")
install.packages("viridis")
install.packages("wesanderson")
install.packages("divo")
install.packages("abdiv")

library(eulerr)
library(BiocManager)
library(edgeR)
library(baySeq)
library(NOISeq)
if(require("parallel")) cl <- makeCluster(8) else cl <- NULL
library(immunarch)
library(assertthat)
library(readr)
library(igraph)
library(viridis)
library(wesanderson)
library(stringr)
library(abdiv)
library(divo)

browseVignettes("baySeq")

# loading data without non-coding and out of frame sequencies

file_path = "C:\\Tereshchenko\\Гугл-диск\\с компа\\Работа\\Treg памяти\\Секвенирование TCR\\Data"
CorrectData = repLoad(file_path)
print(CorrectData$meta, n = 50)
CorrectData


#  filter invivo and invitro data

InVitroData = repFilter(CorrectData, .method = "by.meta", .query = list( inv = include("in vitro")))
print(InVitroData$meta, n = 50)
InVivoData = repFilter(CorrectData, .method = "by.meta", .query = list( inv = include("in vivo")))
print(CorrectData$meta, n = 50)

# Create matrix with clonotype counts for tmm-normalization
ClonesMatrix = function(imdata) {
  LongestElement = which.max(lapply(imdata, nrow))
  ClonesMatrix = matrix(1, nrow = length(imdata[[LongestElement]]$Clones) , ncol = length(imdata))
  for (i in 1:length(imdata)) {
    ClonesMatrix[,i] = c(imdata[[i]]$Clones, rep.int(0, (length(imdata[[LongestElement]]$Clones)-length(imdata[[i]]$Clones))))
  }
  colnames(ClonesMatrix) = names(imdata)
  return(ClonesMatrix)
}

ClonesMatrix(InVitroData$data)#[10000:10100,]
ClonesMatrix(InVivoData$data)#[22000:22100,]

# tmm- normalization of matrix with clonotype counts
NInVitroMatrix = tmm(ClonesMatrix(InVitroData$data)) #invitro
NInVitroMatrix#[6000:6100,]

NInVitroMatrix[NInVitroMatrix < 1 & NInVitroMatrix > 0] = 1
NInVitroMatrix


NInVivoMatrix = tmm(ClonesMatrix(InVivoData$data)) #invivo
NInVivoMatrix[22000:22100,]

NInVivoMatrix[NInVitroMatrix < 1 & NInVitroMatrix > 0] = 1
NInVivoMatrix


# replace original counts with tmm-normalized counts 

Normalized_Replacement = function(imdata, normalized_table) {
  normalized_table[normalized_table < 1 & normalized_table > 0] = 1
  for (i in 1:length(imdata)) {
    clones = as.integer(normalized_table[,i])
    length(clones) = length(imdata[[i]]$Clones)
    imdata[[i]]$Clones = clones
  }
  return(imdata)
}

NInVitroData = Normalized_Replacement(InVitroData$data, NInVitroMatrix) #invitro
NInVitroData = list(NInVitroData, InVitroData$meta)
names(NInVitroData) = c("data", "meta")
NInVitroData
print(InVitroData$meta, n = 50)

NInVivoData = Normalized_Replacement(InVivoData$data, NInVivoMatrix) #invivo
NInVivoData = list(NInVivoData, InVivoData$meta)
names(NInVivoData) = c("data", "meta")
NInVivoData
print(NInVivoData$meta, n = 50)

NCorrectData = c(NInVitroData$data, NInVivoData$data)
NCorrectData = list(NCorrectData, CorrectData$meta)
names(NCorrectData) = c("data", "meta")
NCorrectData



# tmm check
# Clones and Clonotypes Count 
ClonesAndClonotypes = function(imdata) {
  SampleName = c()
  ClonesNum = c()
  ClonotypeNum = c()
  for (i in 1:length(imdata$data)) {
    SampleName =  names(imdata$data)
    ClonesNum = append(ClonesNum, sum(imdata$data[[i]]$Clones))
    ClonotypeNum = append(ClonotypeNum, length(imdata$data[[i]]$Clones))
  }
  print(data.frame(SampleName, ClonesNum, ClonotypeNum, ClonesNum/ClonotypeNum))
}
ClonesAndClonotypes(NInVitroData)
ClonesAndClonotypes(InVitroData)

ClonesAndClonotypes(NInVivoData)
ClonesAndClonotypes(InVivoData)

write.csv(ClonesAndClonotypes(NCorrectData),  "C:\\Users\\sirius\\Desktop\\repertoires.csv") # Supplementary data 1(repertoires)

# join the samples into groups of stimulation
#in vitro
# NBaL0 = repSample(NBaL0, .method = "sample" ,.n = 17952)
# NBaL0 = filter(NBaL0, !duplicated(NBaL0$CDR3.aa, fromLast = F))

NBaL0 = rbind(NInVitroData$data$"Ba-L-05", NInVitroData$data$"Ba-L-09", NInVitroData$data$"Ba-L-12", NInVitroData$data$"Ba-L-13")
NBaL0 = NBaL0[order(NBaL0$Clones, decreasing = TRUE),]
NBaL0


NBaLBa = rbind(NInVitroData$data$"BaLBa-05", NInVitroData$data$"BaLBa-07", NInVitroData$data$"BaLBa-09", NInVitroData$data$"BaLBa-12") 
NBaLBa = NBaLBa[order(NBaLBa$Clones, decreasing = TRUE),]
NBaLBa

NBa = rbind(NBaL0, NBaLBa) 
NBa = NBa[order(NBa$Clones, decreasing = TRUE),]
NBa

NBaLBl = rbind(NInVitroData$data$"BaLBl-05", NInVitroData$data$"BaLBl-07", NInVitroData$data$"BaLBl-09", NInVitroData$data$"BaLBl-10", NInVitroData$data$"BaLBl-12")
NBaLBl = NBaLBl[order(NBaLBl$Clones, decreasing = TRUE),]
NBaLBl


NBaLHum = rbind(NInVitroData$data$"BaLHum-05", NInVitroData$data$"BaLHum-07", NInVitroData$data$"BaLHum-09", NInVitroData$data$"BaLHum-10")
NBaLHum = NBaLHum[order(NBaLHum$Clones, decreasing = TRUE),]
NBaLHum

NBlL0 = rbind(NInVitroData$data$"Bl-L-05", NInVitroData$data$"Bl-L-07", NInVitroData$data$"Bl-L-10", NInVitroData$data$"Bl-L-12", NInVitroData$data$"Bl-L-13")
NBlL0 = NBlL0[order(NBlL0$Clones, decreasing = TRUE),]
NBlL0

NBlLBl = rbind(NInVitroData$data$"BlLBl-05", NInVitroData$data$"BlLBl-07", NInVitroData$data$"BlLBl-10", NInVitroData$data$"BlLBl-12", NInVitroData$data$"BlLBl-13")
NBlLBl = NBlLBl[order(NBlLBl$Clones, decreasing = TRUE),]
NBlLBl

NInVitroList = list(NBaL0, NBaLBa, NBaLHum, NBaLBl, NBlL0, NBlLBl)
names(NInVitroList) = c("BaL0", "BaLBa", "BaLHum", "BaLBl", "BlL0", "BlLBl")
NInVitroList

#in vivo

NInt = rbind(NInVivoData$data$"int-05", NInVivoData$data$"int-06", NInVivoData$data$"int-09")
NInt = NInt[order(NInt$Clones, decreasing = TRUE),]
NInt

NTrans = rbind(NInVivoData$data$"trans-01", NInVivoData$data$"trans-02", NInVivoData$data$"trans-03",NInVivoData$data$"trans-04")
NTrans = NTrans[order(NTrans$Clones, decreasing = TRUE),]
NTrans

NVacBl = rbind(NInVivoData$data$"VacBl-01", NInVivoData$data$"VacBl-02", NInVivoData$data$"VacBl-03", NInVivoData$data$"VacBl-04")
NVacBl = NVacBl[order(NVacBl$Clones, decreasing = TRUE),]
NVacBl

NVacDC = rbind(NInVivoData$data$"VacDC-01", NInVivoData$data$"VacDC-02", NInVivoData$data$"VacDC-03", NInVivoData$data$"VacDC-04", NInVivoData$data$"VacDC-05")
NVacDC = NVacDC[order(NVacDC$Clones, decreasing = TRUE),]
NVacDC

NVacHum = rbind(NInVivoData$data$"VacHum-01", NInVivoData$data$"VacHum-03", NInVivoData$data$"VacHum-04", NInVivoData$data$"VacHum-05")
NVacHum = NVacHum[order(NVacHum$Clones, decreasing = TRUE),]
NVacHum

NInVivoList = list(NInt, NVacHum, NVacBl, NVacDC, NTrans)
names(NInVivoList) = c("Int", "VacHum", "VacBl", "VacDC",  "Trans")
NInVivoList

NList = c(NInVitroList, NInVivoList)
NList


# Creating matrix for fold change expression analysis
FoldChangesMatrix = function(imdata, int){
  FoldChangesMatrix = list()
  for (i in 1:length(imdata)){
    ClonesInt = int[match(imdata[[i]]$CDR3.aa, int$CDR3.aa, nomatch = 0)[match(imdata[[i]]$CDR3.aa, int$CDR3.aa, nomatch = 0)!=0],]$Clones
    ClonesImdata = filter(imdata[[i]], match(imdata[[i]]$CDR3.aa, int$CDR3.aa, nomatch = 0) > 0)$Clones
    FoldChangesMatrix[i] = list(matrix(c(ClonesInt, ClonesImdata), ncol = 2))
    colnames(FoldChangesMatrix[[i]]) = c('Int', names(imdata[i]))
    rownames(FoldChangesMatrix[[i]]) = int$CDR3.aa[match(imdata[[i]]$CDR3.aa, int$CDR3.aa, nomatch = 0)]
  }
  names(FoldChangesMatrix) = names(imdata)
  return(FoldChangesMatrix)
}
NInVitroList
InVitroFCM = FoldChangesMatrix(NInVitroList, NBa) #in vitro
InVitroFCM = InVitroFCM[3:6]
names(InVitroFCM)
InVitroFCM
InVitroFCM$BaLHum

InVivoFCM = FoldChangesMatrix(NInVivoList, NInt) #in vivo
InVivoFCM = InVivoFCM[2:5]
names(InVivoFCM)


InVitroDataFCM = FoldChangesMatrix(NInVitroData$data, NBa)
names(InVitroDataFCM)


# differential expression analysis by Bayesian posterior likelihood estimation

if(require("parallel")) cl <- makeCluster(8) else cl <- NULL
BayesFunction = function(FoldChangesData, Nimdata, x) {
  BayesResults = list()
  Result = list()
  for (i in 1:length(FoldChangesData)) {
    CD <- new("countData", data = FoldChangesData[[i]], replicates = c("Int", "names(FoldChangesData[[i]])"), groups =  list(NDE = c(1,1), DE = c(1,2)))
    libsizes(CD) <- getLibsizes(CD)
    CD <- getPriors.NB(CD, samplesize = length(FoldChangesData[[i]][,1]), estimation = "QL", cl = cl)
    CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = FALSE)
    BayesResults[[i]] = topCounts(CD, group = "DE", number = length(FoldChangesData[[i]][,1]))
    BayesResults[[i]] = cbind(BayesResults[[i]], rownames(FoldChangesData[[i]]), as.vector(FoldChangesData[[i]][,2]/FoldChangesData[[i]][,1]))
    BayesResults[[i]] = subset(BayesResults[[i]], BayesResults[[i]][,8]>= x & BayesResults[[i]][,3]>=0.95)
    Result[[i]] = filter(Nimdata[[i]], Nimdata[[i]]$CDR3.aa %in% BayesResults[[i]][,7])
    Result[[i]] = subset(Result[[i]], Result[[i]]$Clones>= x )
  }
  names(BayesResults) = names(FoldChangesData)
  names(Result) = names(Nimdata)
  return(Result)
}

fold_increase_generation = function() {
  InVitroBayes2 <<- BayesFunction(InVitroFCM, NInVitroList[3:6], 2)
  InVitroBayes4 <<- BayesFunction(InVitroFCM, NInVitroList[3:6], 4)
  InVitroBayes8 <<- BayesFunction(InVitroFCM, NInVitroList[3:6], 8)
  InVitroBayes16 <<- BayesFunction(InVitroFCM, NInVitroList[3:6], 16)
  InVitroBayes32 <<- BayesFunction(InVitroFCM, NInVitroList[3:6], 32)
  
  InVivoBayes2 <<- BayesFunction(InVivoFCM, NInVivoList[2:5], 2)
  InVivoBayes4 <<- BayesFunction(InVivoFCM, NInVivoList[2:5], 4)
  InVivoBayes8 <<- BayesFunction(InVivoFCM, NInVivoList[2:5], 8)
  InVivoBayes16 <<- BayesFunction(InVivoFCM, NInVivoList[2:5], 16)
  InVivoBayes32 <<- BayesFunction(InVivoFCM, NInVivoList[2:5], 32)
}

fold_increase_generation()

# Nonexpanded lists

Non_expanded_lists = function(immdata_list, bayes_data){
  for (i in 1:length(bayes_data)) {
    immdata_list[[i+2]]= filter(immdata_list[[i+2]], !immdata_list[[i+2]]$CDR3.aa %in% bayes_data[[i]]$CDR3.aa)
  }
  return(immdata_list[3:6])
}


non_exp_invitro_2x = Non_expanded_lists(NInVitroList, InVitroBayes2)

non_exp_invivo_2x = Non_expanded_lists(NInVivoList, InVivoBayes2)




# Overlap
repOverlap(NCorrectData$data, .method = "overlap", .verbose = F)%>%  #  FIG. S1 
  vis("heatmap2", .by = c("Mice"), .meta = NCorrectData$meta)


Fig1A = repOverlap(NInVitroList, .method = "morisita", .verbose = F)%>%  #  Fig1A 
  vis(.title = "",.axis.text.size = 30, .text.size = 10, .signif.digits = 2, .legend = F) + scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
Fig1A

Fig1B = repOverlap(NInVivoList, .method = "morisita", .verbose = F)%>%  #  Fig 1B 
  vis(.title = "",.axis.text.size = 30, .text.size= 10, .signif.digits = 2, .legend = F ) + scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

Fig1A + Fig1B


FigS3_InVitro2 = repOverlap(InVitroBayes2, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2,.title = "") + scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVitro4 = repOverlap(InVitroBayes4, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "" )+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVitro8 = repOverlap(InVitroBayes8, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "" )+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVitro16 = repOverlap(InVitroBayes16, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "" )+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVitro32 = repOverlap(InVitroBayes32, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2 ,.title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))


FigS3_InVivo2 = repOverlap(InVivoBayes2, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVivo4 = repOverlap(InVivoBayes4, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVivo8 = repOverlap(InVivoBayes8, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVivo16 = repOverlap(InVivoBayes16, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))
FigS3_InVivo32 = repOverlap(InVivoBayes32, .method = "morisita", .verbose = F)%>%  #  Fig S3
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

(Fig1A+Fig1B)/ (FigS3_InVitro2 + FigS3_InVivo2)/ (FigS3_InVitro4 + FigS3_InVivo4)/(FigS3_InVitro8 + FigS3_InVivo8)/ (FigS3_InVitro16 + FigS3_InVivo16)/(FigS3_InVitro32 + FigS3_InVivo32) #Fig 1

FigS2_1x= repOverlap(NList, .method = "morisita", .verbose = F)%>%  #  Fig S2 
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "") + scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

FigS2_2x = repOverlap(c(InVitroBayes2, InVivoBayes2), .method = "morisita", .verbose = F)%>%  #  Fig S2 
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "") + scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

FigS2_4x = repOverlap(c(InVitroBayes4, InVivoBayes4), .method = "morisita", .verbose = F)%>%  #  Fig S2 
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

FigS2_8x = repOverlap(c(InVitroBayes8, InVivoBayes8), .method = "morisita", .verbose = F)%>%  #  Fig S2 
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

FigS2_16x = repOverlap(c(InVitroBayes16, InVivoBayes16), .method = "morisita", .verbose = F)%>%  #  Fig S2 
  vis(.legend = F,  .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

FigS2_32x = repOverlap(c(InVitroBayes32, InVivoBayes32), .method = "morisita", .verbose = F)%>%  #  Fig S2 
  vis(.legend = F, .axis.text.size = 30, .text.size = 11, .signif.digits = 2, .title = "")+ scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))


(FigS2_1x+FigS2_2x+FigS2_4x) / (FigS2_8x+FigS2_16x+FigS2_32x) # Fig S3

Fig2SA = repOverlap(non_exp_invitro_2x, .method = "morisita", .verbose = F)%>%    
  vis(.title = "In vitro", .axis.text.size = 30, .text.size = 10, .signif.digits = 2, .legend = F) + scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))


Fig2SB = repOverlap(non_exp_invivo_2x, .method = "morisita", .verbose = F)%>%  #  Fig 1B 
  vis(.title = "In vivo",.axis.text.size = 30, .text.size= 10, .signif.digits = 2, .legend = F ) + scale_fill_gradient2(low="#053061" ,mid="#f7f7f7"  , high="#67001f", midpoint=0.2, limits=c(0, 0.4))

Fig2SA + Fig2SB # FigS2

#igraph and GLIPH2

# preparing data for graph
immdata_graph = function(imdata){
  imdata = filter(imdata, !duplicated(imdata$CDR3.aa, fromLast = F))
  imdata_matrix = matrix(data = c("aa", "aa"), ncol =2, nrow = 1)
  for (i in 1:length(imdata$Clones)){
    for (j in 1:length(imdata$Clones)) {
      if (nchar(imdata$CDR3.aa[[i]]) == nchar(imdata$CDR3.aa[[j]])) {
        if (sum(!(strsplit(imdata$CDR3.aa[[i]], "")[[1]] == strsplit(imdata$CDR3.aa[[j]], "")[[1]]) , na.rm = TRUE) == 1 & !(imdata$CDR3.aa[[j]] %in% imdata_matrix[,1])) {
          imdata_matrix = rbind(imdata_matrix, c(imdata$CDR3.aa[[i]], imdata$CDR3.aa[[j]]))}
      }
    }
  }

  return(imdata_matrix[-1,])
}

# functionalize 

plot_clustered_graph = function(imdata,imgraph) {
  node_size = log2(filter(imdata, imdata$CDR3.aa %in% V(imgraph)$name)$Clones)
  colors <- hcl.colors(max(clusters(imgraph)$membership), palette = "spectral")
  set.seed(20)
  return(plot(imgraph,
              vertex.color=colors[membership(clusters(imgraph))], 
              vertex.size=node_size, 
              vertex.label = NA, 
              main = substitute(imdata), 
              layout = layout.fruchterman.reingold(imgraph)))
}



create_GLIPH_input = function(imdata, file_name){
  imdata = imdata %>%  mutate(across('J.name', str_replace, "\\*.*", ''))
  imdata = imdata %>%  mutate(across('V.name', str_replace, "\\*.*", ''))
  GLIPH2_TABLE = matrix(data = c(imdata$CDR3.aa, 
                                 imdata$V.name, 
                                 imdata$J.name, 
                                 rep.int(NA, times = length(imdata$CDR3.aa)),
                                 rep.int(paste0(substitute(file_name), ": Public"), times = length(imdata$CDR3.aa)),
                                 imdata$Clones),
                        nrow = length(imdata$CDR3.aa), ncol = 6)
  GLIPH2_TABLE = as.data.frame(GLIPH2_TABLE)
  write_delim(GLIPH2_TABLE, file = paste0("C:\\Users\\sirius\\Desktop\\GLIPH2\\", substitute(file_name), ".txt"), col_names = F, delim = "\t")
  print("Written to C:\\Users\\sirius\\Desktop\\GLIPH2\\")
}



get_GLIPH_motifs = function(file_name){
  GLIPH2_Result = read.csv(file = paste0("C:\\Users\\sirius\\Desktop\\GLIPH2\\",  substitute(file_name),".csv"))
  GLIPH2_Result=table(GLIPH2_Result$pattern)
  GLIPH2_motifs = names(GLIPH2_Result[GLIPH2_Result >=3])
  GLIPH2_motifs = str_replace_all(GLIPH2_motifs, "%", ".")
  GLIPH2_motifs =GLIPH2_motifs[!GLIPH2_motifs == "single"]
  GLIPH2_motifs =GLIPH2_motifs[str_length(GLIPH2_motifs) >= 4] 
  return(GLIPH2_motifs)
}


plot_GLIPH_clusters = function(imdata,imgraph, GLIPH_motifs){
  GLIPH_clusters = motif_clusters(imgraph, GLIPH_motifs)
  motif_colors <- c("#f7f7f7", hcl.colors(max(GLIPH_clusters$membership), palette = "spectral"))
  node_size = log2(filter(imdata, imdata$CDR3.aa %in% V(imgraph)$name)$Clones)
  set.seed(20)
  plot(imgraph,
       vertex.color=motif_colors[membership(GLIPH_clusters)], 
       vertex.size=node_size, 
       vertex.label = NA, 
       main = substitute(imdata) ,
       layout = layout.fruchterman.reingold(imgraph))
  
}

motif_clusters = function(imgraph, motif_vector){
  clusters_vector = rep_len(1, length.out = length(V(imgraph)))
  names(clusters_vector) = V(imgraph)$name
  for (i in 1:length(motif_vector)) {
    clusters_vector[which(str_detect(names(V(imgraph)), motif_vector[i]))] = i+1
  }
  motif_communities = make_clusters(imgraph, membership = clusters_vector)
  return(motif_communities)
}

stimulation_clusters = function(imgraph, imdata1, imdata2, imdata3){
  clusters_vector = rep_len(1, length.out = length(V(imgraph)))
  names(clusters_vector) = V(imgraph)$name
  for (i in 1:length(imdata1$CDR3.aa)) {
    clusters_vector[which(str_detect(names(V(imgraph)), imdata1$CDR3.aa[i]))] = 1
  }
  for (i in 1:length(imdata2$CDR3.aa)) {
    clusters_vector[which(str_detect(names(V(imgraph)), imdata2$CDR3.aa[i]))] = 2
  }
  for (i in 1:length(imdata3$CDR3.aa)) {
    clusters_vector[which(str_detect(names(V(imgraph)), imdata3$CDR3.aa[i]))] = 3
  }
  motif_communities = make_clusters(imgraph, membership = clusters_vector)
  return(motif_communities)
}


plot_stim_clusters = function(imdata, imgraph, stim_clusters){
  stim_colors <- c('#A71B4B', "#F9BE56", '#584B9F' )
  node_size = log2(filter(imdata, imdata$CDR3.aa %in% V(imgraph)$name)$Clones)
  set.seed(20)
  plot(imgraph,
       vertex.color=stim_colors[membership(stim_clusters)], 
       vertex.size=node_size, 
       vertex.label = NA, 
       main = substitute(imdata) ,
       layout = layout.fruchterman.reingold(imgraph))
}


# Graphs_plotting
#Direct_graphs

Direct_df = rbind(InVitroBayes2$BlL0 ,InVivoBayes2$VacDC, InVivoBayes2$Trans)
Direct_df
Direct_Graph_Data = immdata_graph(Direct_df)
Direct_Graph_Data
Direct_Graph = graph_from_edgelist(Direct_Graph_Data, directed = F)
Direct_Graph
direct_graph_5 <- delete_vertices(Direct_Graph, V(Direct_Graph)[clusters(Direct_Graph)$membership %in% which(clusters(Direct_Graph)$csize <= 5)])
direct_graph_5

plot_clustered_graph(Direct_df, direct_graph_5) # Fig 2A

create_GLIPH_input(Direct_df, "Direct")

direct_motifs = get_GLIPH_motifs(Direct)
direct_motifs_с = direct_motifs[!(direct_motifs %in% in_direct_motifs)]
direct_motifs_с
write.csv(direct_motifs,  "C:\\Users\\sirius\\Desktop\\GLIPH2\\direct_motifs.csv") # Supplementary data 2


plot_GLIPH_clusters( Direct_df,direct_graph_5, direct_motifs) # Fig 2B

direct_stim_cluster = stimulation_clusters(direct_graph_5, InVitroBayes2$BlL0, InVivoBayes2$VacDC, InVivoBayes2$Trans)

plot_stim_clusters(Direct_df, direct_graph_5 ,direct_stim_cluster) # Fig 2C
direct_clusters = motif_clusters(direct_graph_5, direct_motifs)
stim_colors <- c('#A71B4B', "#F9BE56", '#584B9F' )
text(1.2, seq(0.5,-0.5, length=3), labels = c("BlL0", "VacDC", "Trans"), 
       col = stim_colors, cex = 2 )



# indirect graphs
InDirect_df = rbind(InVitroBayes2$BaLBl ,InVivoBayes2$VacBl, InVivoBayes2$Trans)
InDirect_df

InDirect_Graph_Data = immdata_graph(InDirect_df)
InDirect_Graph_Data

InDirect_Graph = graph_from_edgelist(InDirect_Graph_Data, directed = F)
InDirect_Graph
InDirect_Graph5 = delete_vertices(InDirect_Graph, V(InDirect_Graph)[clusters(InDirect_Graph)$membership %in% which(clusters(InDirect_Graph)$csize <= 5)])
InDirect_Graph5

plot_clustered_graph(InDirect_df , InDirect_Graph5) # Fig 2D

create_GLIPH_input(InDirect_df, "InDirect")
in_direct_motifs = get_GLIPH_motifs(InDirect)
in_direct_motifs_с = in_direct_motifs[!(in_direct_motifs %in% direct_motifs)]
in_direct_motifs_с
write.csv(in_direct_motifs,  "C:\\Users\\sirius\\Desktop\\GLIPH2\\in_direct_motifs.csv") # Supplementary data 2
plot_GLIPH_clusters(InDirect_df,InDirect_Graph5, in_direct_motifs) # Fig 2E

indirect_stim_cluster = stimulation_clusters(InDirect_Graph5, InVitroBayes2$BaLBl, InVivoBayes2$VacBl, InVivoBayes2$Trans)

plot_stim_clusters(InDirect_df, InDirect_Graph5 ,indirect_stim_cluster) # Fig 2F
text(1.2, seq(0.5,-0.5, length=3), labels = c("BaLBl", "VacBl", "Trans"), 
     col = stim_colors, cex = 2 )

# int graph
Int_df = rbind(NBaL0 ,NBaLBa, NInt)
Int_df = filter(Int_df, !duplicated(Int_df$CDR3.aa))
Int_df = Int_df[1:3200,] 
Int_df


Int_Graph_Data = immdata_graph(Int_df)
Int_Graph_Data

Int_Graph = graph_from_edgelist(Int_Graph_Data, directed = F)
Int_Graph


Int_Clusters_above_5 <- delete_vertices(Int_Graph, V(Int_Graph)[clusters(Int_Graph)$membership %in% which(clusters(Int_Graph)$csize <= 5)])
clusters(Int_Clusters_above_5)$csize
Int_Clusters_above_5 

plot_clustered_graph(Int_df , Int_Clusters_above_5) # Fig 2G

create_GLIPH_input(Int_df, "Int")

int_motifs = get_GLIPH_motifs(Int)
length(int_motifs)
int_motifs
write.csv(int_motifs,  "C:\\Users\\sirius\\Desktop\\GLIPH2\\int_motifs.csv") # Supplementary data 2

plot_GLIPH_clusters(Int_df,Int_Clusters_above_5, int_motifs) # Fig 2H

int_stim_cluster = stimulation_clusters(Int_Clusters_above_5, NBaL0 ,NBaLBa, NInt)


plot_stim_clusters(Int_df, Int_Clusters_above_5 ,int_stim_cluster) # Fig 2I
text(1.2, seq(0.5,-0.5, length=3), labels = c("BaL0", "BaLBa", "Int"), 
     col = stim_colors, cex = 2 )



# samples_in_clusters
samples_in_clusters = function(motifs, imdata_common, imdata1, imdata2, imdata3){
  sources = 1:length(motifs)
  for (i in 1:length(motifs)) {
    clonotypes = imdata_common[str_which(imdata_common$CDR3.aa, motifs[i]),]$CDR3.nt
    sources[i] = as.numeric(any(clonotypes %in% imdata1$CDR3.nt)) + as.numeric(any(clonotypes %in% imdata2$CDR3.nt)) + as.numeric(any(clonotypes %in% imdata3$CDR3.nt))
  }
  return(table(sources))
}
samples_in_clusters(direct_motifs_с, Direct_df, InVitroBayes2$BlL0, InVivoBayes2$VacDC, InVivoBayes2$Trans)
samples_in_clusters(in_direct_motifs_с, InDirect_df, InVitroBayes2$BaLBl, InVivoBayes2$VacBl, InVivoBayes2$Trans)
samples_in_clusters(int_motifs, Int_df, NInVitroList$BaL0, NInVitroList$BaLBa, NInVivoList$Int)

# venn

plot(
  venn(
    list(direct= direct_motifs,indirect = in_direct_motifs, unstimulated = int_motifs)) , 
  fills = c('#A71B4B', "#F9BE56", '#584B9F'), labels = list(col = "black", fontsize = 20), quantities = list(col = "black", fontsize = 20))
) # Fig.4S


#pathway clonotypes

get_pathaway_clonotypes = function(motifs, imdata) {
  numbers = c()
  for (i in 1:length(motifs)) {
    numbers = append(numbers, str_which(imdata$CDR3.aa, motifs[i]))
  }
  numbers = unique(numbers)
  
  return(imdata[numbers,])
}


direct_clonotypes = get_pathaway_clonotypes(direct_motifs_с, Direct_df)
direct_clonotypes
in_direct_clonotypes = get_pathaway_clonotypes(in_direct_motifs_с, InDirect_df)
in_direct_clonotypes
int_clonotypes = get_pathaway_clonotypes(int_motifs, Int_df)
int_clonotypes

# alloreactive Clonotypes from bayes

get_pathaway_clonotypes_from_samples = function(motifs_direct, motifs_indirect, imdata) {
  new_imdata = list()
  
  numbers = c()
  for (i in 1:length(motifs_direct)){
    numbers = append(numbers, str_which(imdata[[3]]$CDR3.aa, motifs_direct[i]))
  }
  numbers = unique(numbers)
  new_imdata[[1]] = imdata[[3]][numbers,]
  numbers = c()
  
  for (i in 1:length(motifs_direct)){
    numbers = append(numbers, str_which(imdata[[7]]$CDR3.aa, motifs_direct[i]))
  }
  numbers = unique(numbers)
  new_imdata[[2]] = imdata[[7]][numbers,]
  numbers = c()
  
  for (i in 1:length(motifs_direct)){
    numbers = append(numbers, str_which(imdata[[8]]$CDR3.aa, motifs_direct[i]))
  }
  numbers = unique(numbers)
  new_imdata[[3]] = imdata[[8]][numbers,]
  numbers = c()
  
  
  for (i in 1:length(motifs_indirect)){
      numbers = append(numbers, str_which(imdata[[2]]$CDR3.aa, motifs_indirect[i]))
}
    numbers = unique(numbers)
    new_imdata[[4]] = imdata[[2]][numbers,]
    numbers = c()
    
  for (i in 1:length(motifs_indirect)){
      numbers = append(numbers, str_which(imdata[[6]]$CDR3.aa, motifs_indirect[i]))
    }
    numbers = unique(numbers)
    new_imdata[[5]] = imdata[[6]][numbers,]
    numbers = c()
    
  for (i in 1:length(motifs_indirect)){
    numbers = append(numbers, str_which(imdata[[8]]$CDR3.aa, motifs_indirect[i]))
    }
    numbers = unique(numbers)
    new_imdata[[6]] = imdata[[8]][numbers,]
    numbers = c()
    
  names(new_imdata) = c("BlL0", "VacDC","Trans direct", "BaLBl", "VacBl","Trans indirect")
  return(new_imdata)
}  

get_pathaway_clonotypes_from_int = function() {
  new_imdata = list()
  
  numbers = c()
  for (i in 1:length(direct_motifs_с)){
    numbers = append(numbers, str_which(NInVivoList[[1]]$CDR3.aa, direct_motifs_с[i]))
  }
  new_imdata[[1]] = NInVivoList[[1]][numbers,]
  
  numbers = c()
  for (i in 1:length(in_direct_motifs_с)){
    numbers = append(numbers, str_which(NInVivoList[[1]]$CDR3.aa, in_direct_motifs_с[i]))
  }
  new_imdata[[2]] = NInVivoList[[1]][numbers,]
  names(new_imdata) = c("Int direct","Int indirect")
  return(new_imdata)
}
pathway_int =get_pathaway_clonotypes_from_int()
pathway_int

#violin plots

violin_plot = function(joined_sample, pathway_sample){
  plot = matrix(c(joined_sample$Proportion[joined_sample$Clones>1], pathway_sample$Proportion, rep(1, times = length(joined_sample$Proportion[joined_sample$Clones>1])), rep(2, times = length(pathway_sample$Proportion))) , nrow = (length(joined_sample$Proportion[joined_sample$Clones>1])+length(pathway_sample$Proportion)), ncol=2)
  colnames(plot) = c("freq", "sample")
  plot = as.data.frame(plot)
  plot$sample = as.factor(plot$sample)
  plot = ggplot(plot, aes(x= sample, y= freq, fill= sample)) + 
    geom_violin()+scale_y_continuous(trans='log10')+ 
    scale_x_discrete(labels = c('Whole','Allo'))+
    stat_summary(fun.y=median, geom="point", size=10, color='#A71B4B') + 
    scale_fill_manual(values=c("#F9BE56", '#584B9F')) + 
    theme_classic()+ 
    theme(legend.position="none", text = element_text(size = 70)) +
    labs(title=match.call()[[3]],x="", y = "log10(frequency)")
  return(plot)
}



(violin_plot(NInVitroList$BlL0, pathway_bayes2$BlL0) + violin_plot(NInVivoList$VacDC, pathway_bayes2$VacDC) + violin_plot(NInVivoList$Trans, pathway_bayes2$`Trans direct`))/
  (violin_plot(NInVitroList$BaLBl, pathway_bayes2$BaLBl) + violin_plot(NInVivoList$VacBl, pathway_bayes2$VacBl) + violin_plot(NInVivoList$Trans, pathway_bayes2$`Trans indirect`))

violin_plot(NInVivoList$Int, pathway_int$`Int direct`[pathway_int$`Int direct`$Clones >1,]) + violin_plot(NInVivoList$Int, pathway_int$`Int indirect`[pathway_int$`Int indirect`$Clones >1,])

(violin_plot(NInVitroList$BlL0, pathway_bayes2$BlL0) + violin_plot(NInVitroList$BaLBl, pathway_bayes2$BaLBl))/
  (violin_plot(NInVivoList$VacDC, pathway_bayes2$VacDC)+ violin_plot(NInVivoList$VacBl, pathway_bayes2$VacBl))/
  (violin_plot(NInVivoList$Trans, pathway_bayes2$`Trans direct`) + violin_plot(NInVivoList$Trans, pathway_bayes2$`Trans indirect`))/
  (violin_plot(NInVivoList$Int, pathway_int$`Int direct`[pathway_int$`Int direct`$Clones >1,]) + violin_plot(NInVivoList$Int, pathway_int$`Int indirect`[pathway_int$`Int indirect`$Clones >1,])) # Fig 3 B
  
# Fig 3 B statistics
wilcox.test(NInVitroList$BlL0$Proportion, pathway_bayes2$BlL0$Proportion)

wilcox.test(NInVitroList$BaLBl$Proportion, pathway_bayes2$BaLBl$Proportion)

wilcox.test(NInVivoList$VacDC$Proportion, pathway_bayes2$VacDC$Proportion) 

wilcox.test(NInVivoList$Trans$Proportion, pathway_bayes2$`Trans direct`$Proportion) 

wilcox.test(NInVivoList$VacBl$Proportion, pathway_bayes2$VacBl$Proportion) 

wilcox.test(NInVivoList$Trans$Proportion, pathway_bayes2$`Trans indirect`$Proportion) 

wilcox.test(NInVivoList$Int[NInVivoList$Int$Clones>1, ]$Proportion, pathway_int$`Int direct`[pathway_int$`Int direct`$Clones >1,]$Proportion)

wilcox.test(NInVivoList$Int[NInVivoList$Int$Clones>1, ]$Proportion, pathway_int$`Int indirect`[pathway_int$`Int indirect`$Clones >1,]$Proportion) 
  

(median(pathway_bayes2$BlL0$Proportion)/median(NInVitroList$BlL0$Proportion) +
median(pathway_bayes2$VacDC$Proportion)/median(NInVivoList$VacDC$Proportion)+
median(pathway_bayes2$`Trans direct`$Proportion)/median(NInVivoList$Trans$Proportion))/3



(median(pathway_bayes2$BaLBl$Proportion)/median(NInVitroList$BaLBl$Proportion)+
median(pathway_bayes2$VacBl$Proportion)/median(NInVivoList$VacBl$Proportion)+
median(pathway_bayes2$`Trans indirect`$Proportion)/median(NInVivoList$Trans$Proportion))/3



hist = c(sum(pathway_bayes2$BlL0$Proportion)/sum(NInVitroList$BlL0$Proportion),
         sum(pathway_bayes2$VacDC$Proportion)/sum(NInVivoList$VacDC$Proportion),
         sum(pathway_bayes2$`Trans direct`$Proportion)/sum(NInVivoList$Trans$Proportion),
         
         
         sum(pathway_bayes2$BaLBl$Proportion)/sum(NInVitroList$BaLBl$Proportion),
         sum(pathway_bayes2$VacBl$Proportion)/sum(NInVivoList$VacBl$Proportion),
         sum(pathway_bayes2$`Trans indirect`$Proportion)/sum(NInVivoList$Trans$Proportion))


hist*100  # Fig 3 A 



# VDJ usage



c(InVitroBayes2, InVivoBayes2)
pathway_bayes2 = get_pathaway_clonotypes_from_samples(direct_motifs_с, in_direct_motifs_с, c(InVitroBayes2, InVivoBayes2))
geneUsage(pathway_bayes2, "musmus.trbv", .norm = T) %>%   vis() 

pathway_bayes4 = get_pathaway_clonotypes_from_samples(direct_motifs_с, in_direct_motifs_с, c(InVitroBayes4, InVivoBayes4))
geneUsage(pathway_bayes4, "musmus.trbv", .norm = T) %>%   vis() 

pathway_bayes8 = get_pathaway_clonotypes_from_samples(direct_motifs_с, in_direct_motifs_с, c(InVitroBayes8, InVivoBayes8))
geneUsage(pathway_bayes8, "musmus.trbv", .norm = T) %>%   vis() 

pathway_bayes16 = get_pathaway_clonotypes_from_samples(direct_motifs_с, in_direct_motifs_с, c(InVitroBayes16, InVivoBayes16))
geneUsage(pathway_bayes16, "musmus.trbv", .norm = T) %>%   vis() 





# Fig 4 rework
# direct
(geneUsage(NInVitroList$BlL0, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40)))) /
  geneUsage(NInVivoList$VacDC, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))/
  geneUsage(NInVivoList$Trans, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))) 


#geneUsage(int_clonotypes, "musmus.trbv", .norm = T) %>%  subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%   subset(Clones > 0.01) %>%    vis()/
(geneUsage(pathway_bayes2$BlL0, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))/
geneUsage(pathway_bayes2$VacDC, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))/
geneUsage(pathway_bayes2$`Trans direct`, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40)))))

# indirect
geneUsage(NInVitroList$BaLBl, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))/
  geneUsage(NInVivoList$VacBl, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))/
  geneUsage(NInVivoList$Trans, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))


#geneUsage(NInVivoList$Int, "musmus.trbv", .norm = T) %>%  subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%   subset(Clones > 0.01) %>%    vis()/
geneUsage(pathway_bayes2$BaLBl, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))/
geneUsage(pathway_bayes2$VacBl, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))/
geneUsage(pathway_bayes2$`Trans indirect`, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40))))

geneUsage(pathway_bayes2, "musmus.trbv", .norm = T) %>% geneUsageAnalysis( .method = "cosine", .verbose = F) %>%   vis(.axis.text.size = 40, .text.size = 20, .labs = c(" ", " "), .title = "Gene usage cosine similarity", .legend = F) # Fig S5 rework

# int_pathway??????
geneUsage(pathway_int$`Int direct`, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40)))) /
  geneUsage(pathway_int$`Int indirect`, "musmus.trbv", .norm = T)%>% subset( Names %in% c('TRBV13-1', 'TRBV13-2','TRBV13-3',"TRBV19", "TRBV31", "TRBV5" ) )  %>%  subset(Clones > 0.01)  %>%   vis(.add.layer = list(scale_y_continuous(limits = c(0, 0.35)),  theme(text = element_text(size = 40)))) 
  
# fig S5B ??????
repOverlap(pathway_bayes2, .method = "cosine", .verbose = F, .col = "aa") %>% vis(.axis.text.size = 40, .text.size = 20, .labs = c(" ", " "), .title = "", .legend = F)

table(pathway_bayes2$`Trans direct`$CDR3.aa %in% pathway_bayes2$BlL0$CDR3.aa)
pathway_bayes2$`Trans direct`$CDR3.aa %in% pathway_bayes2$`Trans indirect`$CDR3.aa
pathway_bayes2$`Trans direct`$CDR3.aa


#nonexpanded
 
geneUsage(non_exp_invitro_2x, "musmus.trbv", .norm = T) %>%  vis()
geneUsage(non_exp_invivo_2x, "musmus.trbv", .norm = T) %>%  vis()


#pathwway list
PathWayList = list(int_clonotypes, direct_clonotypes, in_direct_clonotypes)
names(PathWayList) = c("unstim", "direct", "indirect")
PathWayList

# Supplementary data 3
write.csv2(PathWayList$direct, file = "C:\\Users\\sirius\\Desktop\\GLIPH2\\Directly alloreactive clonotypes.csv")
write.csv2(PathWayList$indirect, file = "C:\\Users\\sirius\\Desktop\\GLIPH2\\Indirectly alloreactive clonotypes.csv")

#spectratype
PathWayList$direct = PathWayList$direct  %>%  mutate(across('V.name', str_replace, "\\*.*", '')) 
PathWayList$indirect = PathWayList$indirect  %>%  mutate(across('V.name', str_replace, "\\*.*", '')) 
PathWayList$unstim = PathWayList$unstim  %>%  mutate(across('V.name', str_replace, "\\*.*", ''))
NInVivoList$Int = NInVivoList$Int  %>%  mutate(across('V.name', str_replace, "\\*.*", ''))



Fig.5A = spectratype(PathWayList$direct, .quant = "id", .col = "aa+v") %>%  vis() 
Fig.5A

colSums(Fig.3C[Fig.3C$Length > 14,3])/colSums(Fig.3C[Fig.3C$Length < 15,3])


Fig.5B = spectratype(PathWayList$indirect, .quant = "id", .col = "aa+v") %>% vis(.title = "D. Indirect alloreactive clonotypes spectrotype")
Fig.5B

Fig.3E = spectratype(PathWayList$unstim, .quant = "id", .col = "aa+v") %>% vis(.title = "E. Spectrotype of unstimulated clonotypes")
Fig.3E

Fig.5C = spectratype(NInVivoList$Int, .quant = "id", .col = "aa+v") # %>% vis(.title = "F. Spectrotype of intact mice repertoire")
Fig.5C = subset(Fig.5C, 9 < Fig.3F$Length & Fig.3F$Length < 19)
Fig.5C = vis(Fig.5C, .title = "F. Spectrotype of intact mice repertoire", .add.layer = theme(text = element_text(size = 40)))
Fig.5C
colSums(Fig.5C[Fig.5C$Length > 14,3])/colSums(Fig.5C[Fig.5C$Length < 15,3])

Fig.5A + Fig.5B + Fig.5C


#k-mers and motifs

getKmers(PathWayList$direct, 14) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq") /   # Fig 4 A-C
    getKmers(PathWayList$indirect, 14) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")/
     getKmers(NInVivoList$Int, 14) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq") #/
  getKmers(PathWayList$unstim, 14) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")
  
getKmers(PathWayList$direct, 14) %>% kmer_profile(.method = "prob")
getKmers(PathWayList$indirect, 14) %>% kmer_profile(.method = "prob")
getKmers(NInVivoList$Int, 14) %>% kmer_profile(.method = "prob")
  
# for 14 aa CDR3 length
(getKmers(PathWayList$direct[ str_length(PathWayList$direct$CDR3.aa) == 14,], 14) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")/
  getKmers(PathWayList$indirect[ str_length(PathWayList$indirect$CDR3.aa) == 14,], 14) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq"))/
  getKmers(NInVivoList$Int[str_length(NInVivoList$Int$CDR3.aa) == 14,], 14) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")

getKmers(PathWayList$direct[ str_length(PathWayList$direct$CDR3.aa) == 14,], 14) %>% kmer_profile(.method = "prob")
getKmers(PathWayList$indirect[ str_length(PathWayList$indirect$CDR3.aa) == 14,], 14) %>% kmer_profile(.method = "prob")
getKmers(NInVivoList$Int[str_length(NInVivoList$Int$CDR3.aa) == 14,], 14) %>% kmer_profile(.method = "prob")

(getKmers(PathWayList$direct[ str_length(PathWayList$direct$CDR3.aa) == 13,], 13) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")/
    getKmers(PathWayList$indirect[ str_length(PathWayList$indirect$CDR3.aa) == 13,], 13) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq"))/
  getKmers(NInVivoList$Int[str_length(NInVivoList$Int$CDR3.aa) == 13,], 13) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")

(getKmers(PathWayList$direct[ str_length(PathWayList$direct$CDR3.aa) == 15,], 15) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")/
    getKmers(PathWayList$indirect[ str_length(PathWayList$indirect$CDR3.aa) == 15,], 15) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq"))/
  getKmers(NInVivoList$Int[str_length(NInVivoList$Int$CDR3.aa) == 15,], 15) %>% kmer_profile(.method = "prob") %>% vis(.plot = "seq")


#database search


vdjdb = dbLoad("C:\\Tereshchenko\\Гугл-диск\\с компа\\Работа\\Treg памяти\\Секвенирование TCR\\Статья\\db\\VDJdb.tsv", "vdjdb")
vdjdb


McPAS = dbLoad("C:\\Tereshchenko\\Гугл-диск\\с компа\\Работа\\Treg памяти\\Секвенирование TCR\\Статья\\db\\McPAS.csv", "McPAS")
McPAS



# the clonotypes found in dbs

annotate_imdata = function(imdata, vdj_db = vdjdb, McPAS_db =McPAS) {
  CDR3_vdjdb = imdata[imdata$CDR3.aa %in% vdj_db$CDR3,]$CDR3.aa
  epitope_vdjdb = c()
  antigen_vdjdb = c()
  species_vdjdb = c()
  vdjdb_matrix = matrix(NA, nrow = length(CDR3_vdjdb), ncol = 5)
  vdjdb_matrix[,1] = CDR3_vdjdb
  vdjdb_matrix[,2] = imdata[imdata$CDR3.aa %in% vdj_db$CDR3,]$CDR3.nt
  for (i in 1:length(CDR3_vdjdb)){
    epitope_vdjdb = append(epitope_vdjdb, (vdj_db[which(vdj_db$CDR3 == CDR3_vdjdb[i]), 9]$Epitope)[1])
    antigen_vdjdb = append(antigen_vdjdb, (vdj_db[which(vdj_db$CDR3 == CDR3_vdjdb[i]), 10]$`Epitope gene`)[1])
    species_vdjdb = append(species_vdjdb, (vdj_db[which(vdj_db$CDR3 == CDR3_vdjdb[i]), 11]$`Epitope species`)[1])
    }
  vdjdb_matrix[,3] = epitope_vdjdb
  vdjdb_matrix[,4] = antigen_vdjdb
  vdjdb_matrix[,5] = species_vdjdb
  
  CDR3_mcpas = imdata[imdata$CDR3.aa %in% McPAS_db$CDR3.beta.aa,]$CDR3.aa
  epitope_mcpas = c()
  antigen_mcpas = c()
  species_mcpas = c()
  mcpas_matrix = matrix(NA, nrow = length(CDR3_mcpas), ncol = 5)
  mcpas_matrix[,1] = CDR3_mcpas
  mcpas_matrix[,2] = imdata[imdata$CDR3.aa %in% McPAS_db$CDR3.beta.aa,]$CDR3.nt
  for (i in 1:length(CDR3_mcpas)){
    epitope_mcpas = append(epitope_mcpas, (McPAS_db[which(McPAS_db$CDR3.beta.aa == CDR3_mcpas[i]), 13]$Epitope.peptide)[1])
    antigen_mcpas = append(antigen_mcpas, (McPAS_db[which(McPAS_db$CDR3.beta.aa == CDR3_mcpas[i]), 11]$Antigen.protein)[1])
    species_mcpas = append(species_mcpas, (McPAS_db[which(McPAS_db$CDR3.beta.aa == CDR3_mcpas[i]), 5]$Pathology)[1])
  }
  mcpas_matrix[,3] = epitope_mcpas
  mcpas_matrix[,4] = antigen_mcpas
  mcpas_matrix[,5] = species_mcpas
  
  for (i in 1:length(CDR3_vdjdb)){
    if (!CDR3_vdjdb[i] %in% CDR3_mcpas){
      mcpas_matrix = rbind(mcpas_matrix, vdjdb_matrix[i,])
    }
  }
  colnames(mcpas_matrix) = c("CDR3.aa", "CDR3.nt" ,"Epitope", "Antigen", "Pathology")
  mcpas_matrix[,5] = str_replace(mcpas_matrix[,5], "Systemic Lupus Erythematosus [:punct:]SLE[:punct:]", "SLE")
  mcpas_matrix[,5] = str_replace(mcpas_matrix[,5], "MCMV", "mCMV")
  mcpas_matrix[,5] = str_replace(mcpas_matrix[,5], "PlasmodiumBerghei", "Plasmodium berghei")
  mcpas_matrix[,5] = str_replace(mcpas_matrix[,5], "Vesicular stomatitis virus", "VSV")
  return( mcpas_matrix)
}

direct_annotation = annotate_imdata(PathWayList$direct)
direct_annotation
write.csv(direct_annotation, "C:\\Tereshchenko\\Гугл-диск\\с компа\\Работа\\Treg памяти\\Секвенирование TCR\\Статья\\annotation\\direct_clonotypes.csv") # Supplementary data 3


indirect_annotation = annotate_imdata(PathWayList$indirect)
write.csv(indirect_annotation, "C:\\Tereshchenko\\Гугл-диск\\с компа\\Работа\\Treg памяти\\Секвенирование TCR\\Статья\\annotation\\indirect_clonotypes.csv") # Supplementary data 3

unstim_annotation = annotate_imdata(PathWayList$unstim)
write.csv(unstim_annotation, "C:\\Tereshchenko\\Гугл-диск\\с компа\\Работа\\Treg памяти\\Секвенирование TCR\\Статья\\annotation\\unstim_annotation.csv") # Supplementary data 3
unstim_annotation

int_annotation = annotate_imdata(NInVivoList$Int)
write.csv(int_annotation, "C:\\Tereshchenko\\Гугл-диск\\с компа\\Работа\\Treg памяти\\Секвенирование TCR\\Статья\\annotation\\int_annotation.csv") # Supplementary data 3




pie(table(direct_clonotypes$CDR3.aa %in% direct_annotation[,1]), labels = c("Not found", "Found"), main = "Direct

    Found/Not found = 0.079", clockwise = F, col = c("#F9BE56", '#584B9F'), init.angle = 90, radius = 1) #Fig 5 A
table(direct_clonotypes$CDR3.aa %in% direct_annotation[,1]) # 79/999 = 0.079 

pie(table(in_direct_clonotypes$CDR3.aa %in% indirect_annotation[,1]), labels = c("Not found", "Found"), main = "Indirect 

     Found/Not found = 0.061", clockwise = F, col = c("#F9BE56", '#584B9F'), init.angle = 90, radius = 1) #Fig 5 B
table(in_direct_clonotypes$CDR3.aa %in% indirect_annotation[,1]) # 111/1828 = 0.061

pie(table(int_clonotypes$CDR3.aa %in% unstim_annotation[,1]), labels = c("Not found", "Found"), main = "Unstimulated 

    Found/Not found = 0.041", clockwise = F, col = c("#F9BE56", '#584B9F'), init.angle = 90, radius = 1) #Fig 5 C

pie(table(NInVivoList$Int$CDR3.aa %in% int_annotation[,1]), labels = c("Not found", "Found"), main = "Intact mice repertoire 

    Found/Not found = 0.011", clockwise = F, col = c("#F9BE56", '#584B9F'), init.angle = 90, radius = 1)  #Fig 5 H
table(NInVivoList$Int$CDR3.aa %in% int_annotation[,1]) # 836/72273 = 0.011 


pie(table(direct_annotation[,4]), clockwise = T, col = hcl.colors(length(table(direct_annotation[,4])), palette = "spectral"), main = "Direct", cex = 1.5)  #Fig 5 E
pie(table(indirect_annotation[,4]), clockwise = T, col = hcl.colors(length(table(indirect_annotation[,4])), palette = "spectral"), main = "Indirect", cex = 1.5) #Fig 5 F
pie(table(unstim_annotation[,4]), clockwise = T, col = hcl.colors(length(table(unstim_annotation[,4])), palette = "spectral"), main = "Unstimulated", cex = 1.5) #Fig 5 G
pie(table(int_annotation[,4]), clockwise = T, col = hcl.colors(length(table(int_annotation[,4])), palette = "spectral"), main = "Intact mice repertoire", 
labels = c(NA,NA, "Diabetes Type 1", NA,NA,NA,NA, "Influenza", NA,NA, "mCMV", NA,NA,NA,NA, "Plasmodium berghei", "RSV", NA,NA, "SLE    ", NA, "
                                                                                                                                                                     Tumor", NA, "     West Nile virus"), cex = 1.5) #Fig 5 H


