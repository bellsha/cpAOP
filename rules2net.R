#generating the assocation rules and edge files

library(doBy)
library(igraph);library(arules)
#
info<-read.table("../Output/NetRun1501/RunInfo-Liver-NetRun1501.txt", sep='\t', header=FALSE, na.strings='', quote="\"")
colnames(info)<-c("runID", "fileName","SPECIES","ORGAN_ID", "CHEMICAL", "DOSE_LEVEL", "SACRIFICE_PERIOD","TEST-TYPE","SINGLE_REPEAT_TYPE", "folder", "treatment", "experiment")

#the annotation dataframe
anno<-read.table("../Output/NetRun1501/Update/PathEnrichReactomeAnnoTable.txt", sep='\t', header=TRUE) #from https://github.com/bellsha/cpAOP/blob/master/gene2biospace.R
path<-read.table("../Files/Reactome/ReactomePathways2UniProtRAT.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE)

#the reactome pathway mapping
#need tohave check.names=FALSE to avoid issues
pwdf2<-read.table("../Output/NetRun1501/Update/PathEnrichReactomeDataFrame.txt", sep='\t', header=TRUE, check.names=FALSE)
#####################################################
#####################
#now bring in pathology phenotype and the clinical chem data
labSum<-read.table(file="../Data/Lab/LiverLabDescData201411.txt", sep='\t', header=TRUE)
LiverSum<-read.table(file="../Data/Pathology/open_tggates_pathologyDESC.txt", sep='\t', header=TRUE)
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
write.table(df.gp,file="../Output/NetRun1501/Update/LiverReactomePhenoEnrichRulesGraph.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )

######################################3
#adding in chemicals->pathways (quick)
treat2chems<-merge(unique(info[,c(5:7,11,9)]),unique(anno[,c(2,1,8,10,9)]), by.x='treatment', by.y='CLID')
write.table(treat2chems, file="../Output/NetRun1501/Update/LiverReactomeEnrichChemicalEdges.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
################
########################################################
#getting additonal associations between the different reactome groups
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
write.table(edgeG,file="../Output/NetRun1501/Update/LiverReactomeReactomeEnrichRulesGraph.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
#################################################################################
#getting additonal associations between the different phenotypes
#asking what is the probability of RHS given LHS
rulesL.pp<-subset(rulesL, subset =lhs %in% phenoItemsL & lift >1 & rhs %in% phenoItemsL)
df.pp<-as(rulesL.pp, "data.frame")
df.pp$rules<-as.character(df.pp$rules)
temp<-sapply(strsplit(df.pp$rules, split="[{}]"), unlist)
df.pp$lhs<-temp[2,]
df.pp$rhs<-temp[4,]
#stripping the '=1' for ease in xref
#for merged networks I need to strip the '=1' from the interactions
df.pp$lhs<-gsub('=1','',df.pp$lhs)
df.pp$rhs<-gsub('=1','',df.pp$rhs)
df.pp$interest<-df.pp$support*df.pp$oddsRatio
#remove those with less than 3 for support, 
#becaucse  line of support was a factor added for the reactome x pathology
df.pp<-subset(df.pp, support >= 3/nrow(dataBinS))
#now to remove redundent edges (bc the OR doesnt really matter here, just support)
g<-graph.data.frame(df.pp[,c('lhs','rhs','support')], directed=FALSE, vertices=NULL)
#make sure that the names are transfered bc apparently they are not
#V(g)$name<-colnames(pcc2)
#removes unconnected nodes(vertices)
#new_g <- delete.vertices(g, which(degree(g) < 1))
edgeG<-as.data.frame(get.edgelist(g))
#double check output expected
#edgeG[1,]
edgeG$support<-E(g)$support
edgeG<-unique(edgeG)
write.table(edgeG,file="../Output/NetRun1501/Update/LiverPhenoPhenoEnrichRulesGraph.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )

