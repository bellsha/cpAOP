#generate the combined edge and node files for visualization in cytoscape
#this is across all the datasets and data types
#this file includes, TG gates chemical, phenotype, and gene expression (mapped to Reactome)
#as well as the toxcast data and the phenotypes mapped to AO files
#format for the edge file
#V1  V2  Type  other attributes......
#format for the node files
#ID  altID  Group  UpperClass  LowerClass

#resource files
chemMap<-read.table("../Files/Chem2MechMapping1013.txt", sep='\t', header=TRUE) #totally not really helpful right now
info<-read.table("../Output/NetRun1501/RunInfo-Liver-NetRun1501.txt", sep='\t', header=FALSE, na.strings='', quote="\"")
colnames(info)<-c("runID", "fileName","SPECIES","ORGAN_ID", "CHEMICAL", "DOSE_LEVEL", "SACRIFICE_PERIOD","TEST-TYPE","SINGLE_REPEAT_TYPE", "folder", "treatment", "experiment")
path<-read.table("../Files/Reactome/ReactomePathways2UniProtRAT.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE)
#Need to expand so i can do something on a class lower bases...right now using the AO for that
#aop<-read.table("../Files/LiverPhenotypeLabelsV2.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE)
aop<-read.table("../Files/LiverPhenotypeLabels.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE)
#remove kidney 
#aop<-subset(aop, Source != "kidney")
#toxcast mapping file...not disable quoting here
tclab<-read.table("../Files/ToxCast/ChemAssayBioHitLabels20150703.txt", sep='\t',na.strings='', header=TRUE, quote="")
#need to get the labels for chemicals uniform until i get everything back to a cas number
tcctemp<-subset(tclab, Group=="Chem"); tcntemp<-subset(tclab, Group !="Chem")
tcctemp$altID<-tolower(tcctemp$altID)
tcctemp$altID<-gsub(" ","_",tcctemp$altID)
tclab<-rbind(tcctemp, tcntemp)
##only want the bioprocess (not reactome) and assay here)
#tclab<-subset(tclab, Group %in% c('Toxcast','BioProc'))
#need to get out the relevent chemicals?

#the reactome mapping file
#name2ID<-read.table("/home/sbell/Desktop/datadrive/Reactome/ReactomePathwayName2IDxSpecv51.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE)


labSum<-read.table(file="../Data/Lab/LiverLabDescData201411.txt", sep='\t', header=TRUE)
LiverSum<-read.table(file="../Data/Pathology/open_tggates_pathologyDESC.txt", sep='\t', header=TRUE)
# probe<-read.table("../Files/rat2302.probe.entrez.go_20150113.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=FALSE)
# colnames(probe)<-c("ProbeID","ENTREZID","GOID", "Evi","GOprocess", "UniprotID", "GOTerm", "GODef")

#bring in the data files with the edges
#pathology to chemical maping: tgpathologyprep201411.R
p2c<-read.table("../Output/NetRun1410/LiverPathologyChemicalEdges201411.txt", sep='\t', header=TRUE, quote="\"")
colnames(p2c)<-c("treatment", "CHEMICAL","SACRIFICE_PERIOD","DOSE_LEVEL","SINGLE_REPEAT_TYPE","Out")
#
#lab to chemical maping: TGlabDesc201411.R
c2l<-read.table("../Output/NetRun1410/LiverLabChemEdged201411.txt", sep='\t', header=TRUE, quote="\"")
colnames(c2l)<-c("treatment", "CHEMICAL","SINGLE_REPEAT_TYPE","SACRIFICE_PERIOD","DOSE_LEVEL","Out","Change")

##
#
#reactome to chemcial mapping: rules2net201501.R
r2c<-read.table("../Output/NetRun1501/Update/LiverReactomeEnrichChemicalEdges.txt", sep='\t', header=TRUE, quote="\"")
#phenotype (path and lab) to reactome mapping: rules2net201501.R
p2r<-read.table("../Output/NetRun1501/Update/LiverReactomePhenoEnrichRulesGraph.txt", sep='\t', header=TRUE, quote="\"")
r2r<-read.table("../Output/NetRun1501/Update/LiverReactomeReactomeEnrichRulesGraph.txt", sep='\t', header=TRUE, quote="\"")
p2p<-read.table("../Output/NetRun1501/Update/LiverPhenoPhenoEnrichRulesGraph.txt", sep='\t', header=TRUE, quote="\"")

#bring in the toxcast data (note that this is humanmapped to rat reactome pathways): 
tc<-read.table("../Files/ToxCast/ToxCastChemAssayRatPathNodes20150703.txt", sep='\t', header=TRUE, quote="\"")
colnames(tc)<-gsub("Type","TYPE",colnames(tc))
#TC to tg mapping
c2c<-read.table("../Output/NetRun1501/TC2TGgatesChemMap20150204.txt", sep='\t', header=TRUE, quote="\"")
#now to slim out the t2r to get rid of some of the unmapped chemicals
t2r<-subset(tc, TYPE != "TCA2TCC")
tmp<-subset(tc, V1 %in% c2c$V1)
t2r<-rbind(t2r, tmp)
# #removing the chemcial associations for now
# t2r<-subset(t2r, Type != "TCA2TCC")

#########
#reactome single & repeat
#single exposure (first 24hrs)
p2cS<-subset(p2c, SINGLE_REPEAT_TYPE =='Single')
r2cS<-subset(r2c, SINGLE_REPEAT_TYPE =='Single')
l2cS<-subset(c2l, SINGLE_REPEAT_TYPE =='Single')
prc.single<-unique(c(as.character(p2cS$CHEMICAL), as.character(p2cS$Out), as.character(l2cS$CHEMICAL), as.character(l2cS$Out),as.character(r2cS$CHEMICAL), as.character(r2cS$Label)))
p2rS<-subset(p2r, lhs %in% prc.single & rhs %in% prc.single)
#long exposure
p2cR<-subset(p2c, SINGLE_REPEAT_TYPE =='Repeat')
r2cR<-subset(r2c, SINGLE_REPEAT_TYPE =='Repeat')
l2cR<-subset(c2l, SINGLE_REPEAT_TYPE =='Repeat')
prc.repeat<-unique(c(as.character(p2cR$CHEMICAL), as.character(p2cR$Out), as.character(l2cR$CHEMICAL), as.character(l2cR$Out), as.character(r2cR$CHEMICAL), as.character(r2cR$Label)))
p2rR<-subset(p2r, lhs %in% prc.repeat & rhs %in% prc.repeat)
#so to get all edges (I think I want a seperate label if it is short vs long exposure)
p2rS$SINGLE_REPEAT_TYPE <-'Single'
p2rR$SINGLE_REPEAT_TYPE <-'Repeat'
p2r2<-rbind(p2rS, p2rR)
p2r2$TYPE<-"pheno2reactome"
#extracting the single/repeat dosing
r2rS<-subset(r2r, V1 %in% prc.single & V2 %in% prc.single)
r2rR<-subset(r2r, V1 %in% prc.repeat & V2 %in% prc.repeat)
#so to get all edges (I think I want a seperate label if it is short vs long exposure)
r2rS$SINGLE_REPEAT_TYPE <-'Single'
r2rR$SINGLE_REPEAT_TYPE <-'Repeat'
r2r2<-rbind(r2rS, r2rR)
r2r2$TYPE<-"reactome2reactome"
#
#extracting the single/repeat dosing
p2pS<-subset(p2p, V1 %in% prc.single & V2 %in% prc.single)
p2pR<-subset(p2p, V1 %in% prc.repeat & V2 %in% prc.repeat)
#so to get all edges (I think I want a seperate label if it is short vs long exposure)
p2pS$SINGLE_REPEAT_TYPE <-'Single'
p2pR$SINGLE_REPEAT_TYPE <-'Repeat'
p2p2<-rbind(p2pS, p2pR)
p2p2$TYPE<-"pheno2pheno"
#
#defining the phenotype to AO edges
#want catagory, label
a2p<-aop[,c(3,1)]; colnames(a2p)<-c("V1", "V2"); 

#specifying the edge type

p2c$TYPE<-"pheno2chem"
r2c$TYPE<-"reactome2chem"
c2l$TYPE<-"pheno2chem"
a2p$TYPE<-"ao2pheno"

#all networks combine
#standardize relevent column names and adding add'l blanks
p2r2v2<-p2r2[,c('lhs','rhs','SINGLE_REPEAT_TYPE','TYPE','support','interest')]
colnames(p2r2v2)<-c('V1','V2','SINGLE_REPEAT_TYPE','TYPE','support','interest')
p2r2v2$Padj<-NA; p2r2v2$UniProtID<-NA;p2r2v2$treatment<-NA; p2r2v2$DOSE_LEVEL<-NA; p2r2v2$SACRIFICE_PERIOD<-NA;p2r2v2$Change<-NA 
#
colnames(p2c)<-c('treatment','V1','SACRIFICE_PERIOD','DOSE_LEVEL','SINGLE_REPEAT_TYPE','V2','TYPE')
p2c$support<-NA;p2c$interest<-NA;p2c$Padj<-NA; p2c$UniProtID<-NA; p2c$Change<-NA
#
colnames(r2c)<-c('treatment','V1','DOSE_LEVEL','SACRIFICE_PERIOD','SINGLE_REPEAT_TYPE','V2','Padj','Name','UniProtID','TYPE')
r2c$support<-NA;r2c$interest<-NA; r2c$Change<-NA
r2c<-subset(r2c, select=-Name)
#
colnames(c2l)<-c('treatment','V1','SINGLE_REPEAT_TYPE','SACRIFICE_PERIOD','DOSE_LEVEL','V2','Change','TYPE')
c2l$support<-NA;c2l$interest<-NA; c2l$Padj<-NA; c2l$UniProtID<-NA
#
r2r2$Padj<-NA; r2r2$UniProtID<-NA;r2r2$treatment<-NA; r2r2$interest<-NA; r2r2$DOSE_LEVEL<-NA; r2r2$SACRIFICE_PERIOD<-NA; r2r2$Change<-NA
#
p2p2$Padj<-NA; p2p2$UniProtID<-NA;p2p2$treatment<-NA; p2p2$interest<-NA; p2p2$DOSE_LEVEL<-NA; p2p2$SACRIFICE_PERIOD<-NA; p2p2$Change<-NA
#
a2p$Padj<-NA; a2p$UniProtID<-NA;a2p$treatment<-NA; a2p$DOSE_LEVEL<-NA;a2p$SACRIFICE_PERIOD<-NA;a2p$Change<-NA;a2p$support<-NA;a2p$interest<-NA; a2p$SINGLE_REPEAT_TYPE<-NA 
#
t2r$Padj<-NA;t2r$treatment<-NA; t2r$DOSE_LEVEL<-NA;t2r$SACRIFICE_PERIOD<-NA;t2r$Change<-NA;t2r$support<-NA;t2r$interest<-NA; t2r$SINGLE_REPEAT_TYPE<-NA 
#
c2c$Padj<-NA;c2c$treatment<-NA; c2c$DOSE_LEVEL<-NA;c2c$SACRIFICE_PERIOD<-NA;c2c$Change<-NA;c2c$support<-NA;c2c$interest<-NA; c2c$SINGLE_REPEAT_TYPE<-NA; c2c$UniProtID<-NA; 

#merge
allEdges<-rbind(p2r2v2, p2c, r2c, c2l, r2r2,p2p2, a2p, t2r, c2c)
allEdges$ID<-paste("A", c(1:nrow(allEdges)), sep='')
#add in edges to join TCchems to tggated chems

# #create a short exposure network (edges that show up at 24hrs or less)
# shortEdges<-subset(allEdges, SINGLE_REPEAT_TYPE =='Single')
# #shortEdges$ID<-paste("S", c(1:nrow(shortEdges)), sep='')
# #create a network of long/repeat exposure (after the first 24hrs)
# longEdges<-subset(allEdges, SINGLE_REPEAT_TYPE =='Repeat')
# #longEdges$ID<-paste("L", c(1:nrow(longEdges)), sep='')
# #because cytoscape is incapable of actually importing things correctly,
# #a node file and a seperate edge file
# allNode<-allEdges[,c("V1","V2","ID")]
# shortNode<-shortEdges[,c("V1","V2","ID")]
# longNode<-longEdges[,c("V1","V2","ID")]
#save the files
write.table(allEdges, "../Output/NetRun1501/Update/LiverTCAReactomePhenoChemGraph.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)

#write.table(allNode, "../Output/NetRun1501/LiverGOReactomePhenoChemNode20150202.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
##
###################################################
#Generate the label files
## ID  Group  ClassUpper ClassLower
#get in the chemical mapping file (for class)
chemMap$CHEMICAL<-tolower(chemMap$NAME)
uChemMap<-unique(subset(chemMap, CATEGORY=='MECHANISM')[,c('CHEMICAL','VALUE')])
CM<-cbind(ID=uChemMap$CHEMICAL,altID=uChemMap$CHEMICAL, Group='Chem',ClassUpper="chem", ClassLower=as.character(uChemMap$VALUE))

noMech<-setdiff(unique(as.character(info$CHEMICAL)), unique(uChemMap$CHEMICAL))
CMno<-cbind(ID=noMech,altID=noMech, Group='Chem', ClassUpper="chem", ClassLower='ChemUnk')

#
# #get pathology mapping
tmp<-cbind(ID=colnames(LiverSum[,3:ncol(LiverSum)]),altID=colnames(LiverSum[,3:ncol(LiverSum)]), Group='Pheno', ClassUpper='Pathology', ClassLower='Pathology')
#he clinical chem mapping
tmp2<-cbind(ID=colnames(labSum)[colnames(labSum) !='treatment'],altID=colnames(labSum)[colnames(labSum) !='treatment'], Group='Pheno', ClassUpper="Lab", ClassLower='Lab')
PD<-as.data.frame(rbind(tmp, tmp2))
#
#get pathway/reactome mapping
RPW<-cbind(unique(path[,c(1:2,8:9)]), Group='Reactome'); colnames(RPW)<-c('ID','altID','ClassUpper','ClassLower','Group')
RPW<-RPW[,c('ID','altID','Group','ClassUpper', 'ClassLower')]
#
#Get the AO mapping to the phenotypes
# this is updated for the new input file
# colnames(aop)
# [1] "Label"       "Source"      "Category"    "Description" "Type" 
# get the ao mappings
# tmp<-unique(aop[,c("AO","class")]); tmp2<-cbind(ID=as.character(tmp$AO), altID = as.character(tmp$AO), Group= "AO", ClassUpper="AO", ClassLower=as.character(tmp$class))
# 
# aoMap<-cbind(ID=as.character(aop$Label),altID=as.character(aop$Description), Group="Pheno", ClassUpper=as.character(aop$Type), ClassLower=as.character(aop$class))
# aoMap<-as.data.frame(unique(rbind(tmp2, aoMap)))
# aoMap$ClassLower[is.na(aoMap$ClassLower)==TRUE]<-"AO" #for those I havent annotated

#Get the AO mapping to the phenotypes
tmp<-unique(aop$Category); tmp2<-cbind(ID=as.character(tmp), altID = as.character(tmp), Group= "AO", ClassUpper="AO", ClassLower="AO")
#this is mapping the AO to the phenotypes in a bit more detail
#filling out the description
nodesc<-subset(aop, is.na(Description)); desc<- subset(aop, !is.na(Description))
nodesc$Description<-nodesc$Label; aop2<-rbind(nodesc, desc)
aoMap<-cbind(ID=as.character(aop2$Label),altID=as.character(aop2$Description), Group="Pheno", ClassUpper=as.character(aop2$Type), ClassLower=as.character(aop2$Category))
#aoMap<-cbind(ID=as.character(aop$Label),altID=as.character(aop$Label), Group="Pheno", ClassUpper=as.character(aop$Type), ClassLower=as.character(aop$Type))
aoMap<-as.data.frame(unique(rbind(tmp2, aoMap)))


#to get the ones that for some reason arent in the file
x<-setdiff(PD$ID, aoMap$ID)
aoMap2<-rbind(aoMap, subset(PD, ID %in% x))
# the TC mapping are done on file input

labels<-as.data.frame(rbind(CM,CMno, aoMap2, RPW, tclab))
#a check
z<-setdiff(c(as.character(allEdges$V1), as.character(allEdges$V2)), as.character(labels$ID))

#labels
write.table(labels,file="../Output/NetRun1501/Update/LiverTCAPhenoChemReactomeLabels.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
