#the purpose of this script is to format the toxcast assay results
#(if a compond was active or not) and the compond metadata
#it depends on the output from Reactome processing script (Reactome2Network/ReactomeClassv2.R) ReactomePathways2UniProt.txt
#linking the data back to rat for use with TGGates rat data depends on another scripts output (ReactomeXSpecies.R) ReactomePathwayName2IDxSpec.txt
#bring in the assay results. This has a 1 if hit 0 if no, NA if not tested

hit<-read.csv("../Data/NewData/AllResults_hitc_Matrix_141121.csv", header=TRUE, row.names=1)
hit$code<-row.names(hit)
#this gives more information about the chemical. to join on column "code"
chemSum<-read.csv("../Data/NewData/Chemical_Summary_141121.csv", header=TRUE)
#this is additonal info on the assays 
#assaySum$aenm is a subset of column headings in "hit" (code isnt present)
assaySum<-read.csv("../Data/NewData/Assay_Summary_141121.csv", header=TRUE)
# #this is the cytotox data
# cyto<-read.csv("../Data/NewData/AllResults_cyto_dist_141121.csv", header=TRUE)
#assay information
assayInfo<-read.csv("../Data/NewData/ToxCast Assay Annotation Assay_Target_Info_20141021.csv", header=TRUE)
#because they didnt export the file correctly. grrrr
assayInfo<-assayInfo[,2:ncol(assayInfo)]
#note that the aenm is assay_component_endpoint_name in assayInfo and one of the assays appears to be mislabeled
assay<-merge(unique(assaySum[,c(1:4,12,13)]), unique(assayInfo[,c(2:7,9:13)]),by.x="aeid", by.y="aeid", all=TRUE)


#need to generate a matrix where i have a table of chemcial and asssay pairs
#WIP
#This processes the score data
sc<-t(subset(hit, select =-code))
#-1 values indicates that too few conc were used to fit
sc[sc==-1]<-0
sc[is.na(sc)]<-0
n<-matrix(nrow=nrow(sc), ncol=ncol(sc))
for(i in 1:nrow(sc)){
  #this puts the assay name in for wach time a chemical had a "hit"
  n[i,]<-gsub('1',rownames(sc)[i], sc[i,])
}

n[n==0]<-NA
n2<-as.data.frame(n)
acdat<-NULL
acnam<-colnames(sc)
#for(i in 1:2){
for(i in 1:ncol(n2)){
  tmp<-na.omit(n2[,i])
  tmp2<-cbind(Code=rep(acnam[i], length(tmp)), aenm=as.character(tmp))
  acdat<-rbind(acdat, tmp2)
}
acdat<-as.data.frame(acdat); acdat$Type<-'TCA2TCC'; colnames(acdat)<-c('V1','V2','Type')

#need to address the alternate id. for now the alt id is going to be the assay name
assayS<-unique(assay[,c(3,3,8,15)]); colnames(assayS)<-c("ID", "altID", "ClassUpper","ClassLower"); assayS<-cbind(assayS, Group="Toxcast")
chem<-chemSum[,c(1,4)]; colnames(chem)<-c("ID", "altID"); chem<-cbind(chem, Group="Chem", ClassUpper="Chem", ClassLower="Chem")
labels<-rbind(assayS, chem)

#to expand the assay out to link them to other data
#specifically go from assay target to a biological process
library(biomaRt)
#note that intended target sub fam appears to be specific
as1<-assay[,c(3,16,5,8,15)]
#because only human stuff is easy at this point
as1<-na.omit(subset(as1, organism =='human'))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


test<-getBM(attributes= c("entrezgene",'ensembl_gene_id', 'uniprot_swissprot'), mart=ensembl)
test1<-subset(test, entrezgene %in% unique(as1$intended_target_gene_id))

#creating mapping to Reactome
#bring in the data
path<-unique(read.table ("/Reactome/Output53/ReactomePathways2UniProt.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE))
path1<-unique(merge(path[,c(1,2,4,7:9)],unique(test1[,c(1,3)]), by.x="UniProtID", by.y="uniprot_swissprot"))
pathAssay<-merge(as1, path1, by.x="intended_target_gene_id", by.y="entrezgene", all.x=TRUE)
A2R<-unique(subset(pathAssay, is.na(ReactomeID) != TRUE)[,c(2,7,6)]); colnames(A2R)<-c('V1','V2','UniProtID')
A2R$Type<-"Toxcast2Reactome"
#mapping the TC assays to the biological processes etc for those with no reactome pathway
x<-setdiff(unique(as.character(as1$aenm)), as.character(unique(A2R$V1)))
A2P<-subset(pathAssay, aenm %in% x)[,c(2,5,6)]; colnames(A2P)<-c('V1','V2','UniProtID')
A2P$Type<-"Toxcast2Biospace";
#get the nodes
Anodes<-rbind(A2R[,c(1,2,3,4)], A2P[,c(1,2,4,3)])
#write.table(Anodes, file="../Output/Assay2Path20150127.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
#ALL nodes
acdat$UniProtID<-NA
nodes<-rbind(acdat, Anodes)
write.table(nodes, file="../Output/ToxCastChemAssayPathNodes20150703.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
#made the labels file
#we want ALL the labels
x<-c(as.character(nodes$V1), as.character(nodes$V2))
tmp<-unique(subset(path, ReactomeID %in% x)[,c(1,2,8,9)])
colnames(tmp)<-c("ID", "altID","ClassUpper","ClassLower"); tmp$Group<-"Reactome"
#the other are from the bioprocess annotations used as nodes
tmp2<-unique(subset(pathAssay, intended_target_family_sub %in% x)[,c(5,5,4,5)])
colnames(tmp2)<-c("ID", "altID","ClassUpper","ClassLower"); tmp2$Group<-"BioProc"
#and get the ed
labs<-rbind(labels, tmp, tmp2)
subset(labs, is.na(ClassUpper))
labs$ClassUpper<-as.character(labs$ClassUpper)
labs$ClassUpper[is.na(labs$ClassUpper)]<-'Tox21'
# labs$ClassLower<-as.character(labs$ClassLower)
# labs$ClassLower[is.na(labs$ClassLower)]<-'Tox21'

write.table(labs, file="../Output/ChemAssayBioHitLabels20150703.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
######################################
######################################
#because the TGGates stuff is all in Rat space, I need to move the reactomeIDs to that domain
#bring in the mapping file
name2ID<-read.table("/Reactome/Output53/ReactomePathwayName2IDxSpec.txt", sep='\t',comment.char='',quote="\"", header=TRUE)
colnames(name2ID)
#note that I just want columnsfor rat and human
#this means both rat and human have data
name2ID2<-subset(name2ID, is.na(name2ID$ReactomeID.Homo.sapiens) == FALSE & is.na(name2ID$ReactomeID.Rattus.norvegicus) == FALSE)[,c("Name","ReactomeID.Homo.sapiens","ReactomeID.Rattus.norvegicus")]
#note that all the reactome lables are in V2
A2Rr<-subset(A2R, Type =="Toxcast2Reactome")
A2Rr<-merge(name2ID2, A2Rr, by.x="ReactomeID.Homo.sapiens", by.y="V2", all=TRUE)
colnames(A2Rr)<-gsub("ReactomeID.Homo.sapiens", "V2",colnames(A2Rr))
A2Rr<-subset(A2Rr, is.na(V1) == FALSE & is.na(V2) ==FALSE)
#update type to reflect the odd mapping
A2Rr$Type<-"TCHum2ReactRat"
colnames(A2Rr)<-c('V2H',"Name","V2","V1","UniProtID","Type")
#removing the parts where there is no rat equ
A2Rr<-subset(A2Rr, is.na(V2) ==FALSE)
Anodes2<-rbind(A2Rr[,c(4,3,6,5)], A2P[,c(1,2,4,3)])
#ALL nodes
nodes2<-unique(rbind(acdat, Anodes2))
write.table(nodes2, file="../Output/ToxCastChemAssayRatPathNodes20150703.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
