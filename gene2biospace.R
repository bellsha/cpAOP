#this is ment to get the gene expression (changes) based on chemcial
#as the dataset requires descretization it will be based on array/treatment (treatment is easier)
#and dose_level='Control' will not be included in final results as all the values will be zero
#libraries
library(doBy)
library(igraph);library(arules)
#get the differential experss data
deDesc<-read.table("../Output/NetRun1501/Desc2FC-T1-0.75-FC-2-Liver-NetRun1501.txt", sep='\t', header=TRUE) #DE data from https://github.com/bellsha/TGGATESProc/blob/master/CombineArrayDE.R
deDA<-abs(deDesc)
#############################33
#get the experimental mappings
probe<-read.table("../Files/rat2302.probe.entrez.go_20150515.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=FALSE) #from: https://github.com/bellsha/TGGATESProc/blob/master/ProbeAnnotation.R
colnames(probe)<-c("ProbeID","ENTREZID","GOID", "Evi","GOprocess", "UniprotID", "GOTerm", "GODef")
tmp2<-subset(probe, ProbeID %in% row.names(deDesc))
path<-read.table("../Files/Reactome/ReactomePathways2UniProtRAT.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE) #from: https://github.com/bellsha/Reactome2Network/blob/master/ReactomeClassv2.R
#remove the NA
path<-subset(path, UniProtID != "NA")
#fix the whitespace issue
########################
# Abstracting from probe to "biological space"
#there are a few different levels to explore.
#Merging everything to keep that into the options availabl
################
pathMap<-merge(tmp2[,c(1:3,6)], path[,c(1,2,4,7,9)], by.x='UniprotID', by.y='UniProtID')  

############################################
#lets build the pathway enrichment, with each row corresponding to a "community" 
m<-matrix(nrow=nrow(deDA), ncol=ncol(deDA))
for(i in 1:nrow(deDA)){
	m[i,]<-gsub('1',rownames(deDA)[i], deDA[i,])
}
m[m==0]<-NA
m2<-as.data.frame(m)
cdat<-NULL
cnam<-colnames(deDA)
for(i in 1:ncol(m2)){
	tmp<-na.omit(m2[,i])
	tmp2<-cbind(Treat=rep(cnam[i], length(tmp)), Probe=as.character(tmp))
	cdat<-rbind(cdat, tmp2)
}
cdat<-as.data.frame(cdat)
#Calculate the number of de probes for a given treatment
#Note that we only care about the probes for which we have pathway information
cdat2<-subset(cdat, Probe %in% unique(pathMap$ProbeID))
sizes<-summaryBy(Probe~Treat,data=cdat2, FUN=length, keep.names=TRUE) 
cdata<-merge(cdat2, sizes, by.x="Treat", by.y="Treat")
colnames(cdata)<-c('Treat','Probe','Size')
#merging using uniprotID as the identifier to go back to pathMap
cdataM<-merge(cdata, unique(pathMap[,1:2]), by.x="Probe", by.y="ProbeID")
sdataM<-na.omit(cdataM)
###########
#for reactome
source("../Code/HyperEnrich.R")

#by reactome pathway...order is 1st column=clusterID and second column=featureID and thrid column= size
#need to go by id
anno<-ClusReport(sdataM[,c(2,4,3)], background=sdataM$UniprotID, RefData=pathMap[,c(5,1)])
#anno<-ClusReport(sdataM[,c(2,4,3)], background=sdataM$UniprotID, RefData=pathMap[,c(6,1)])
#now make a dataframe for the items
pwdf<-matrix(nrow=length(unique(anno$CLID)), ncol=length(unique(anno$Label)))
ID<-unique(anno$CLID)
labID<-unique(anno$Label)
colnames(pwdf)<-labID
rownames(pwdf)<-ID

for(i in 1:length(ID)){
	tmp<-subset(anno, CLID ==ID[i])
	pwdf[i,]<-labID %in% tmp$Label
}
pwdf2<-as.data.frame(apply(pwdf,2,as.numeric))
pwdf2$treatment<-as.factor(row.names(pwdf))
#adding in the name to the anno file (bc enrichment function doesnt have it)
anno2<-merge(anno, unique(path[,1:2]), by.x="Label", by.y="ReactomeID", all.x=TRUE)
#all the enriched pathways
write.table(anno2,file="../Output/NetRun1501/Update/PathEnrichReactomeAnnoTable.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
write.table(pwdf2,file="../Output/NetRun1501/Update/PathEnrichReactomeDataFrame.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )


############
# #for GO
# #libraries
# source("../Code/OldCode/readGOorgSB130625.R")
# #this is code for the GO analysis
# source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/GOHyperGAll.txt")
# library(doBy)
# library(GOstats)
# library(GO.db)
# library(WGCNA)
# library(qvalue)
# #save the old wd to come back to it
# oldDir<-getwd()
# #create a new wd to host the GO stuff
# dir.create("../Files/GO")
# setwd("../Files/GO")
# annFile<-"../rat2302.probe.entrez.go_20150113.txt"
# #input order is goid, geneid, gocat
# #using uniprotID
# readGOorgvSB(myfile = annFile, colno = c(3,6,5), geneList=NULL)
# bkground <-NULL
# #gene2GOlistSB(rootUK=T)
# if (is.null(bkground)){
#   bkgd<-probe$EntrezID #all the probes considered
# }else{
#   bkgd<-bkground
# }
# gene2GOlistSB(rootUK=T, background=bkgd)
# 
# GO_BP_DF <- read.table(file="GO_BP_DF", header=T, colClasses = "character")
# go_df <- read.table(file="go_df", header=T, colClasses = "character")
# ###
# #go analysis...1st column=ProbeID and 2nd column=clusterID and 3rd column= size
# annoGO<- GOCluster_ReportSB(CL_DF=sdataM[,c(4,2,3)], background=sdataM$UniprotID, method="simplify", id_type="gene", CLSZ=4, cutoff=0.05, gocats="BP")
# #because we want at least 2 genes in the cluster to have the term
# annoGO<-subset(annoGO, SampleMatch >1 & Term != 'NA' & Padj <= 0.05)
# 
# #now make a dataframe for the items
# pwdfg<-matrix(nrow=length(unique(annoGO$CLID)), ncol=length(unique(annoGO$GOID)))
# ID<-unique(annoGO$CLID)
# labID<-unique(annoGO$GOID)
# colnames(pwdfg)<-labID
# rownames(pwdfg)<-ID
# 
# for(i in 1:length(ID)){
#   tmp<-subset(annoGO, CLID ==ID[i])
#   pwdfg[i,]<-labID %in% tmp$GOID
# }
# pwdfGO2<-as.data.frame(apply(pwdfg,2,as.numeric))
# pwdfGO2$treatment<-as.factor(row.names(pwdfg))
# #all the enriched pathways
# setwd(oldDir)
# write.table(annoGO,file="../Output/NetRun1501/PathEnrichGoAnnoTable.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
# write.table(pwdfGO2,file="../Output/NetRun1501/PathEnrichGODataFrame.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
######
# #generate the go to reactome rules...
# #getting the rules
# #combine the datasets
# mDataPath<-merge(pwdf2, pwdfGO2, by.x='treatment',by.y='treatment', all=TRUE)
# dataBin<-subset(mDataPath, select=-treatment)
# rownames(dataBin)<-mDataPath$treatment
# #to avoid infinity OR values
# tmp<-rep(c(0,1), c(ncol(pwdf2)-1, ncol(pwdfGO2)-1))
# tmp2<-rep(c(1,0), c(ncol(pwdf2)-1, ncol(pwdfGO2)-1))
# 
# dataBin2<-rbind(tmp,tmp2, dataBin)
# 
# #remove Rows/Columns which are really sparse...requiring more than 2 observations
# dataBinS<-dataBin2[!apply(dataBin2,1,function(x) sum(abs(x), na.rm=TRUE)<=2),]
# dataBinS<-dataBinS[,!apply(dataBinS,2,function(x) sum(abs(x), na.rm=TRUE)<=2)]
# 
# dataBinS<-as.data.frame(apply(dataBinS,2,as.factor))
# bin.data<-as(dataBinS, "transactions")
# #support is set for 2 occurances in the dataset. Note this only applys to cross type measure
# #for withinin associations manual filtering will be needed
# rulesL<-apriori(bin.data, parameter =list(support = 2/nrow(dataBinS), confidence = 0.1, minlen=2, maxlen=2))
# quality(rulesL)<-cbind(quality(rulesL), oddsRatio = interestMeasure(rulesL, method="oddsRatio", bin.data))
# 
# #
# #inspect(head(sort(rulesL, by = "lift"), n = 3))
# reactItemsL<-intersect(colnames(bin.data), paste(colnames(pwdf2),1, sep='='))
# goItemsL<-intersect(colnames(bin.data), paste(colnames(pwdfGO2),1, sep='='))
# #asking what is the probability of RHS given LHS
# #note this does impact the values....so may want to reconsider later
# rulesL.gr<-subset(rulesL, subset =lhs %in% goItemsL & lift >1 & rhs %in% reactItemsL)
# df.gr<-as(rulesL.gr, "data.frame")
# #currently rules is a factor, setting up to remove redundent rules
# #these will have different C values
# df.gr$rules<-as.character(df.gr$rules)
# temp<-sapply(strsplit(df.gr$rules, split="[{}]"), unlist)
# df.gr$lhs<-temp[2,]
# df.gr$rhs<-temp[4,]
# #stripping the '=1' for ease in xref
# #for merged networks I need to strip the '=1' from the interactions
# df.gr$lhs<-gsub('=1','',df.gr$lhs)
# df.gr$rhs<-gsub('=1','',df.gr$rhs)
# #add in the additional interest measure that is support*oddsratio
# df.gr$interest<-df.gr$support*df.gr$oddsRatio
# write.table(df.gr,file="../Output/NetRun1501/LiverGOReactEnrichRulesGraph20150115.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )





