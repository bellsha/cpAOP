#this is ment to get the gene expression (changes) based on chemcial
#as the dataset requires descretization it will be based on array/treatment (treatment is easier)
#and dose_level='Control' will not be included in final results as all the values will be zero
#libraries
library(doBy)
library(igraph);library(arules)
#get the differential experss data
deData<-read.table(file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/DiffExprDatav2-T1-0-FC-2-Liver-NetRun1410.txt", sep='\t', header=TRUE, comment.char='')
#dataLabels<-read.table("/wikihomedir/sbell/TG-Gates/Data/Microarray_Data_RatLabels_20141030.txt", sep='\t', header=TRUE, na.strings='', quote="\"")
info<-read.table("/wikihomedir/sbell/TG-Gates/Output/NetRun1410/RunInfo-Liver-NetRun1410.txt", sep='\t', header=FALSE, na.strings='', quote="\"")
colnames(info)<-c("runID", "fileName","SPECIES","ORGAN_ID", "CHEMICAL", "DOSE_LEVEL", "SACRIFICE_PERIOD","TEST-TYPE","SINGLE_REPEAT_TYPE", "folder", "treatment", "experiment")
#descritize the data to decrease the search space
deDesc<-deData
#using a cutoff of 2 fold change
deDesc[deDesc <= -1]<--1
deDesc[deDesc > -1]<-0
deDesc[deData >= 1]<-1
#because I think that we will improve resoluition at the pathway level but just looking at changes
#and not direction
deDA<-abs(deDesc)
write.table(deDesc, file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/Desc2FCDiffExprRun1410.txt", sep='\t', col.names=TRUE, row.names=TRUE,quote=FALSE )
#deDesc<-read.table("/home/sbell/Desktop/wikihomedir/sbell/TG-Gates/network/TGGatesServeOut/Desc2FCDiffExprNetRun.txt", sep='\t', header=TRUE)
#deDA<-abs(deDesc)
#############################33
#get the experimental mappings
probe<-read.table("/wikihomedir/sbell/TG-Gates/Files/rat2302.probe.entrez.go_20150113.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=FALSE)
colnames(probe)<-c("ProbeID","ENTREZID","GOID", "Evi","GOprocess", "UniprotID", "GOTerm", "GODef")
tmp2<-subset(probe, ProbeID %in% row.names(deData))
path<-read.table("/wikihomedir/sbell/TG-Gates/Files/ReactomePathways2UniProtRATv51.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE)
#remove the NA
path<-subset(path, UniProtID != "NA")
#fix the whitespace issue
#expanding to make associaitons based on minor pathways
pathMap<-merge(tmp2[,c(1,2,6)], path[,c(1,2,4,9)], by.x='UniprotID', by.y='UniProtID')
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
cdataM<-merge(cdata, unique(pathMap[,1:2]), by.x="Probe", by.y="ProbeID")


source("/wikihomedir/sbell/TG-Gates/Code/HyperEnrich.R")
sdataM<-na.omit(cdataM)
anno<-ClusReport(sdataM[,c(2,4,3)], background=sdataM$UniprotID, RefData=pathMap[,c(6,1)])
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
#all the enriched pathways
write.table(anno,file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/PathEnrichMinorAnnoTable.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )

#####################################################
#####################
#now bring in pathology phenotype and the clinical chem data
labSum<-read.table(file="/wikihomedir/sbell/TG-Gates/Data/Lab/LiverLabDescData201411.txt", sep='\t', header=TRUE)
LiverSum<-read.table(file="/wikihomedir/sbell/TG-Gates/Data/Pathology/open_tggates_pathologyDESC.txt", sep='\t', header=TRUE)
#combined pheno data
phenoData<-merge(labSum, LiverSum[,2:ncol(LiverSum)], by.x='treatment',by.y='treatment', all=TRUE)
#pheno plus pathway data
mDataPath<-merge(pwdf2, phenoData, by.x='treatment',by.y='treatment', all=TRUE)

#############################################3

#getting the rules
dataBin<-subset(mDataPath, select=-treatment)
rownames(dataBin)<-mDataPath$treatment
#to avoid infinity OR values
tmp<-rep(c(0,1), c(ncol(pwdf2)-1, ncol(phenoData)-1))
tmp2<-rep(c(1,0), c(ncol(pwdf2)-1, ncol(phenoData)-1))

dataBin2<-rbind(tmp,tmp2, dataBin)

#remove Rows/Columns which are really sparse...requiring more than 2 observations
dataBinS<-dataBin2[!apply(dataBin2,1,function(x) sum(abs(x), na.rm=TRUE)<=2),]
dataBinS<-dataBinS[,!apply(dataBinS,2,function(x) sum(abs(x), na.rm=TRUE)<=2)]

dataBinS<-as.data.frame(apply(dataBinS,2,as.factor))
#added
#dataBinS[dataBinS==0]<-NA
bin.data<-as(dataBinS, "transactions")
#support is set for 2 occurances in the dataset. Note this only applys to cross type measure
#for reactome:reactome manual filtering will be needed
rulesL<-apriori(bin.data, parameter =list(support = 2/nrow(dataBinS), confidence = 0.1, minlen=2, maxlen=2))
quality(rulesL)<-cbind(quality(rulesL), oddsRatio = interestMeasure(rulesL, method="oddsRatio", bin.data))

#
#inspect(head(sort(rulesL, by = "lift"), n = 3))
geneItemsL<-intersect(colnames(bin.data), paste(colnames(pwdf2),1, sep='='))
phenoItemsL<-intersect(colnames(bin.data), paste(colnames(phenoData),1, sep='='))
#asking what is the probability of RHS given LHS
rulesL.gp<-subset(rulesL, subset =lhs %in% geneItemsL & lift >1 & rhs %in% phenoItemsL)
df.gp<-as(rulesL.gp, "data.frame")
#currently rules is a factor, setting up to remove redundent rules
#these will have different C values
df.gp$rules<-as.character(df.gp$rules)
temp<-sapply(strsplit(df.gp$rules, split="[{}]"), unlist)
df.gp$lhs<-temp[2,]
df.gp$rhs<-temp[4,]
#stripping the '=1' for ease in xref
#for merged networks I need to strip the '=1' from the interactions
df.gp$lhs<-gsub('=1','',df.gp$lhs)
df.gp$rhs<-gsub('=1','',df.gp$rhs)
#add in the additional interest measure that is support*oddsratio
df.gp$interest<-df.gp$support*df.gp$oddsRatio
write.table(df.gp,file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomePhenoEnrichRulesGraph111314.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )

######################################3
#adding in chemicals->pathways (quick)
treat2chems<-merge(unique(info[,c(5:7,11,9)]),unique(anno[,c(1,4,8,9)]), by.x='treatment', by.y='CLID')
write.table(treat2chems, file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomeEnrichChemicalEdges111314.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
################
########################################################
#getting additonal associaitons
#asking what is the probability of RHS given LHS
rulesL.gg<-subset(rulesL, subset =lhs %in% geneItemsL & lift >1 & rhs %in% geneItemsL)
df.gg<-as(rulesL.gg, "data.frame")
#currently rules is a factor, setting up to remove redundent rules
#these will have different C values
df.gg$rules<-as.character(df.gg$rules)
temp<-sapply(strsplit(df.gg$rules, split="[{}]"), unlist)
df.gg$lhs<-temp[2,]
df.gg$rhs<-temp[4,]
#stripping the '=1' for ease in xref
#for merged networks I need to strip the '=1' from the interactions
df.gg$lhs<-gsub('=1','',df.gg$lhs)
df.gg$rhs<-gsub('=1','',df.gg$rhs)
#remove those with less than 3 for support, 
#becaucse  line of support was a factor added for the reactome x pathology
df.gg<-subset(df.gg, support >= 3/nrow(dataBinS))
#now to remove redundent edges (bc the OR doesnt really matter here, just support)
library(igraph)
g<-graph.data.frame(df.gg[,c('lhs','rhs','support')], directed=FALSE, vertices=NULL)
#make sure that the names are transfered bc apparently they are not
#V(g)$name<-colnames(pcc2)
#removes unconnected nodes(vertices)
#new_g <- delete.vertices(g, which(degree(g) < 1))
edgeG<-as.data.frame(get.edgelist(g))
#double check output expected
#edgeG[1,]
edgeG$support<-E(g)$support
edgeG<-unique(edgeG)
write.table(edgeG,file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomeReactomeEnrichRulesGraph11132014.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
#################################################################################
#generate the combined edge and node files for visualization in cytoscape
#bring in the data files with the edges
p2c<-read.table("/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverPathologyChemicalEdges201411.txt", sep='\t', header=TRUE, quote="\"")
colnames(p2c)<-c("treatment", "CHEMICAL","SACRIFICE_PERIOD","DOSE_LEVEL","SINGLE_REPEAT_TYPE","Out")
r2c<-read.table("/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomeEnrichChemicalEdges111314.txt", sep='\t', header=TRUE, quote="\"")
p2r<-read.table("/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomePhenoEnrichRulesGraph111314.txt", sep='\t', header=TRUE, quote="\"")
c2l<-read.table("/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverLabChemEdged201411.txt", sep='\t', header=TRUE, quote="\"")
#short exposure
p2cS<-subset(p2c, SINGLE_REPEAT_TYPE =='Single')
r2cS<-subset(r2c, SINGLE_REPEAT_TYPE =='Single')
prc.single<-unique(c(as.character(p2cS$CHEMICAL), as.character(p2cS$Out), as.character(r2cS$CHEMICAL), as.character(r2cS$Label)))
p2rS<-subset(p2r, lhs %in% prc.single & rhs %in% prc.single)
#long exposure
p2cR<-subset(p2c, SINGLE_REPEAT_TYPE =='Repeat')
r2cR<-subset(r2c, SINGLE_REPEAT_TYPE =='Repeat')
prc.repeat<-unique(c(as.character(p2cR$CHEMICAL), as.character(p2cR$Out), as.character(r2cR$CHEMICAL), as.character(r2cR$Label)))
p2rR<-subset(p2r, lhs %in% prc.repeat & rhs %in% prc.repeat)

#so to get all edges (I think I want a seperate label if it is short vs long exposure)
p2rS$SINGLE_REPEAT_TYPE <-'Single'
p2rR$SINGLE_REPEAT_TYPE <-'Repeat'
p2r2<-rbind(p2rS, p2rR)
p2r2$TYPE<-"pheno2reactome"
#write it
#write.table(p2r2, "/home/sbell/Desktop/wikihomedir/sbell/TG-Gates/network/TGGatesServeOut/LiverRPnoEnrichRulesEDITGraph041014.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
r2r<-read.table("/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomeReactomeEnrichRulesGraph11132014.txt", sep='\t', header=TRUE, quote="\"")
#r2r<- edgeG
r2rS<-subset(r2r, V1 %in% prc.single & V2 %in% prc.single)
r2rR<-subset(r2r, V1 %in% prc.repeat & V2 %in% prc.repeat)
#so to get all edges (I think I want a seperate label if it is short vs long exposure)
r2rS$SINGLE_REPEAT_TYPE <-'Single'
r2rR$SINGLE_REPEAT_TYPE <-'Repeat'
r2r2<-rbind(r2rS, r2rR)
r2r2$TYPE<-"reactome2reactome"
#write it
#write.table(r2r2, "G:/NetworkFiles/network/TGGatesServeOut/LiverReactomeReactomeEnrichRulesEDITGraph0414.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
#x<-c(colnames(r2r),colnames(p2c),colnames(r2c),colnames(p2r2),colnames(c2l))
p2c$TYPE<-"pheno2chem"
r2c$TYPE<-"reactome2chem"
c2l$TYPE<-"pheno2chem"
#all networks combine
#standardize relevent column names and adding add'l blanks
p2r2v2<-p2r2[,c('lhs','rhs','SINGLE_REPEAT_TYPE','TYPE','support','interest')]
colnames(p2r2v2)<-c('V1','V2','SINGLE_REPEAT_TYPE','TYPE','support','interest')
p2r2v2$Padj<-NA; p2r2v2$EntrezIDs<-NA;p2r2v2$treatment<-NA; p2r2v2$DOSE_LEVEL<-NA; p2r2v2$SACRIFICE_PERIOD<-NA;p2r2v2$Change<-NA 
#
colnames(p2c)<-c('treatment','V1','SACRIFICE_PERIOD','DOSE_LEVEL','SINGLE_REPEAT_TYPE','V2','TYPE')
p2c$support<-NA;p2c$interest<-NA;p2c$Padj<-NA; p2c$EntrezIDs<-NA; p2c$Change<-NA
#
colnames(r2c)<-c('treatment','V1','DOSE_LEVEL','SACRIFICE_PERIOD','SINGLE_REPEAT_TYPE','V2','Padj','EntrezIDs','TYPE')
r2c$support<-NA;r2c$interest<-NA; r2c$Change<-NA
#
colnames(c2l)<-c('treatment','V1','SINGLE_REPEAT_TYPE','SACRIFICE_PERIOD','DOSE_LEVEL','V2','Change','TYPE')
c2l$support<-NA;c2l$interest<-NA; c2l$Padj<-NA; c2l$EntrezIDs<-NA

#
r2r2$Padj<-NA; r2r2$EntrezIDs<-NA;r2r2$treatment<-NA; r2r2$interest<-NA; r2r2$DOSE_LEVEL<-NA; r2r2$SACRIFICE_PERIOD<-NA; r2r2$Change<-NA
#merge
allEdges<-rbind(p2r2v2, p2c, r2c, c2l, r2r2)
allEdges$ID<-paste("A", c(1:nrow(allEdges)), sep='')

#create a short exposure network (edges that show up at 24hrs or less)
shortEdges<-subset(allEdges, SINGLE_REPEAT_TYPE =='Single')
#shortEdges$ID<-paste("S", c(1:nrow(shortEdges)), sep='')
#create a network of long/repeat exposure (after the first 24hrs)
longEdges<-subset(allEdges, SINGLE_REPEAT_TYPE =='Repeat')
#longEdges$ID<-paste("L", c(1:nrow(longEdges)), sep='')
#because cytoscape is incapable of actually importing things correctly,
#a node file and a seperate edge file
allNode<-allEdges[,c("V1","V2","ID")]
shortNode<-shortEdges[,c("V1","V2","ID")]
longNode<-longEdges[,c("V1","V2","ID")]
#save the files
write.table(allEdges, "/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomePhenoChemGraph111414.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
#
write.table(allNode, "/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomePhenoChemNode111414.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
#write.table(shortNode, "/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomePhenoChemSHORTNode111414.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
#write.table(longNode, "/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverReactomePhenoChemLONGNode111414.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
##
###################################################
#Generate the label files
## ID  Group  ClassUpper ClassLower
chemMap<-read.table("/wikihomedir/sbell/TG-Gates/Files/Chem2MechMapping1013.txt", sep='\t', header=TRUE)
#get in the chemical mapping file (for class)
chemMap$CHEMICAL<-tolower(chemMap$NAME)
uChemMap<-unique(subset(chemMap, CATEGORY=='MECHANISM')[,c('CHEMICAL','VALUE')])
CM<-cbind(ID=uChemMap$CHEMICAL, Group='Chem',ClassUpper="chem", ClassLower=as.character(uChemMap$VALUE))
noMech<-setdiff(unique(as.character(info$CHEMICAL)), unique(uChemMap$CHEMICAL))
CMno<-cbind(ID=noMech, Group='Chem', ClassUpper="chem", ClassLower='ChemUnk')
#
#get pathology mapping
tmp<-cbind(ID=colnames(LiverSum[,2:ncol(LiverSum)])[colnames(phenoData) !='treatment'], Group='Pheno', ClassUpper='Path', ClassLower='Path')
tmp2<-cbind(ID=colnames(labSum)[colnames(labSum) !='treatment'], Group='Pheno', ClassUpper="Lab", ClassLower='Lab')
PD<-rbind(tmp, tmp2)
#
#get pathway/reactome mapping
#RPW<-cbind(unique(path[,c(2,8,9)]), Group='Reactome'); colnames(RPW)<-c('ID','ClassUpper','ClassLower','Group')
RPW<-cbind(unique(path[,c(9,8,9)]), Group='Reactome'); colnames(RPW)<-c('ID','ClassUpper','ClassLower','Group')
#RPW<-cbind(unique(path[,5:7]), Group='Reactome'); colnames(RPW)<-c('ID','Parent','ClassLower','Group')
RPW<-RPW[,c('ID','Group','ClassUpper', 'ClassLower')]
labels<-as.data.frame(rbind(CM,CMno, PD, RPW))
#labels
write.table(labels,file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverPhenoChemReactomeLabels11132014.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
