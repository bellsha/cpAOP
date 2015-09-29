#TG Gates chemical mapping...
#mapping from tg gates to toxcast

# tgdata<-read.table("../Data/Open-tggates_AllAttribute.txt", sep='\t', header=TRUE, na.strings='', quote="\"")
# tgChem<-unique(tgdata[,c(8:10)])
#to fix the TNFalpha
# tgChem$chem<-gsub("TNF*", 'TNFa',tgChem$COMPOUND_NAME)
# tgChem$chem<-tolower(tgChem$chem)
info<-read.table("../Output/NetRun1501/RunInfo-Liver-NetRun1501.txt", sep='\t', header=FALSE, na.strings='', quote="\"")
colnames(info)<-c("runID", "fileName","SPECIES","ORGAN_ID", "CHEMICAL", "DOSE_LEVEL", "SACRIFICE_PERIOD","TEST-TYPE","SINGLE_REPEAT_TYPE", "folder", "treatment", "experiment")
tgChem<-unique(as.data.frame(cbind(CHEMICAL=as.character(info[,5]), chem=tolower(info$CHEMICAL))))

tcChem<-unique(read.csv("../Files/ToxCast/Chemical_Summary_141121.csv", header=TRUE, na.strings='', quote="\""))
tcChem$chem<-tolower(tcChem$chnm)
tcChem$chem<-gsub(" ","_",tcChem$chem)
#overlap:
length(intersect(tgChem$chem, tcChem$chem))
#[1] 34
#create a mapping file for these in the form of an edgelist:
#V1  V2  Type  other attributes......
c2c<-merge(tcChem, tgChem)
c2c2<-c2c[,c(2,6)]; c2c2$TYPE<-"TCChem2TGChem"; colnames(c2c2)<-c("V1","V2", "TYPE")
write.table(c2c2, "../Output/NetRun1501/TC2TGgatesChemMap20150204.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
