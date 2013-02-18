# Sample code MsigDB Curation with RObject
# Dataset <- add Layer
# project id <- add  Dataset, etc 
###############################################################################
## load the synapse client and login
library(synapseClient)
library(GSA)
library(affy)
synapseLogin()

## set up a project
myName1 <- "C1"
Ann1<-"positional gene sets"
C1 <- GSA.read.gmt("~/CurationPathway/MSigDB/C1/c1.all.v3.1.symbols.gmt")
C1.ALL<-C1[[1]]
names(C1.ALL)<-C1[[2]]

myName2 <- "C2"
Ann2<-"curated gene sets"
C2 <- GSA.read.gmt("~/CurationPathway/MSigDB/C2/c2.all.v3.1.symbols.gmt")
C2.ALL<-C2[[1]]
names(C2.ALL)<-C2[[2]]

myName2.1 <- "C2.CGP"
Ann2.1<-"chemical and genetic perturbations"
CGP <- GSA.read.gmt("~/CurationPathway/MSigDB/C2/CGP/c2.cgp.v3.1.symbols.gmt")
C2.CGP<-CGP[[1]]
names(C2.CGP)<-CGP[[2]]

myName2.2 <- "C2.CP"
Ann2.2<-"canonical pathways"
CP <- GSA.read.gmt("~/CurationPathway/MSigDB/C2/CP/c2.cp.v3.1.symbols.gmt")
C2.CP<-CP[[1]]
names(C2.CP)<-CP[[2]]

myName2_1 <- "C2.CP.KEGG"
Ann2_1<-"KEGG"
KEGG <- GSA.read.gmt("~/CurationPathway/MSigDB/C2/CP/KEGG/c2.cp.kegg.v3.1.symbols.gmt")
C2.CP.KEGG<-KEGG[[1]]
names(C2.CP.KEGG)<-KEGG[[2]]

myName2_2 <- "C2.CP.BIOCARTA"
Ann2_2<-"BIOCARTA"
BIOCARTA <- GSA.read.gmt("~/CurationPathway/MSigDB/C2/CP/BIOCARTA/c2.cp.biocarta.v3.1.symbols.gmt")
C2.CP.BIOCARTA<-BIOCARTA[[1]]
names(C2.CP.BIOCARTA)<-BIOCARTA[[2]]

myName2_3 <- "C2.CP.REACTOME"
Ann2_3<-"REACTOME"
REACTOME <- GSA.read.gmt("~/CurationPathway/MSigDB/C2/CP/REACTOME/c2.cp.reactome.v3.1.symbols.gmt")
C2.CP.REACTOME<-REACTOME[[1]]
names(C2.CP.REACTOME)<-REACTOME[[2]]

myName3 <- "C3"
Ann3<-"motif gene sets"
C3 <- GSA.read.gmt("~/CurationPathway/MSigDB/C3/c3.all.v3.1.symbols.gmt")
C3.ALL<-C3[[1]]
names(C3.ALL)<-C3[[2]]

myName3.1 <- "C3.MIR"
Ann3.1<-"microRNA targets"
MIR <- GSA.read.gmt("~/CurationPathway/MSigDB/C3/MIR/c3.mir.v3.1.symbols.gmt")
C3.MIR<-MIR[[1]]
names(C3.MIR)<-MIR[[2]]

myName3.2 <- "C3.TFT"
Ann3.2<-"transcription factor targets"
TFT <- GSA.read.gmt("~/CurationPathway/MSigDB/C3/TFT/c3.tft.v3.1.symbols.gmt")
C3.TFT<-TFT[[1]]
names(C3.TFT)<-TFT[[2]]


myName4 <- "C4"
Ann4<-"computational gene sets"
C4 <- GSA.read.gmt("~/CurationPathway/MSigDB/C4/c4.all.v3.1.symbols.gmt")
C4.ALL<-C4[[1]]
names(C4.ALL)<-C4[[2]]

myName4.1 <- "C4.CGN"
Ann4.1<-"cancer gene neighborhoods"
CGN <- GSA.read.gmt("~/CurationPathway/MSigDB/C4/CGN/c4.cgn.v3.1.symbols.gmt")
C4.CGN<-CGN[[1]]
names(C4.CGN)<-CGN[[2]]

myName4.2 <- "C4.CM"
Ann4.2<-"cancer modules"
CM <- GSA.read.gmt("~/CurationPathway/MSigDB/C4/CM/c4.cm.v3.1.symbols.gmt")
C4.CM<-CM[[1]]
names(C4.CM)<-CM[[2]]

myName5 <- "C5"
Ann5<-"gene ontology(GO) gene sets"
C5 <- GSA.read.gmt("~/CurationPathway/MSigDB/C5/c5.all.v3.1.symbols.gmt")
C5.ALL<-C5[[1]]
names(C5.ALL)<-C5[[2]]

myName5.1 <- "C5.GO_BP"
Ann5.1<-"GO terms biological process"
GO_BP <- GSA.read.gmt("~/CurationPathway/MSigDB/C5/GO_BP/c5.bp.v3.1.symbols.gmt")
C5.GO_BP<-GO_BP[[1]]
names(C5.GO_BP)<-GO_BP[[2]]

myName5.2 <- "C5.GO_CC"
Ann5.2<-"GO terms cellular components"
GO_CC <- GSA.read.gmt("~/CurationPathway/MSigDB/C5/GO_CC/c5.cc.v3.1.symbols.gmt")
C5.GO_CC<-GO_CC[[1]]
names(C5.GO_CC)<-GO_CC[[2]]

myName5.3 <- "C5.GO_MF"
Ann5.3<-"GO terms molecular functions"
GO_MF <- GSA.read.gmt("~/CurationPathway/MSigDB/C5/GO_MF/c5.mf.v3.1.symbols.gmt")
C5.GO_MF<-GO_MF[[1]]
names(C5.GO_MF)<-GO_MF[[2]]


myName6 <- "C6"
Ann6<-"oncogenic signature gene sets"
C6 <- GSA.read.gmt("~/CurationPathway/MSigDB/C6/c6.all.v3.1.symbols.gmt")
C6.ALL<-C6[[1]]
names(C6.ALL)<-C6[[2]]



projName <- "Pathway Database" #as.character(gsub("-",".",Sys.Date()))

myProj <- createEntity(Project(list(name=projName)))

myFolder <- createEntity(Folder(list(name = paste("MSigDB ",as.character(gsub("-",".",Sys.Date())),sep=""), parentId = myProj$properties$id)))

newEntity <- Data(list(name= "All with Gene Symbol IDs", parentId = myFolder$properties$id))
newEntity<-createEntity(newEntity)          
addObject(newEntity, C1.ALL)
addObject(newEntity, C2.ALL)
addObject(newEntity, C2.CGP)
addObject(newEntity, C2.CP)
addObject(newEntity, C2.CP.KEGG)
addObject(newEntity, C2.CP.BIOCARTA)
addObject(newEntity, C2.CP.REACTOME)
addObject(newEntity, C3.ALL)
addObject(newEntity, C3.MIR)
addObject(newEntity, C3.TFT)
addObject(newEntity, C4.ALL)
addObject(newEntity, C4.CGN)
addObject(newEntity, C4.CM)
addObject(newEntity, C5.ALL)
addObject(newEntity, C5.GO_BP)
addObject(newEntity, C5.GO_CC)
addObject(newEntity, C5.GO_MF)
addObject(newEntity, C6.ALL)
storeEntity(newEntity)
