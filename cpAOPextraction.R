#the aim of this script is to extract cpAOP networks
#this will include the extraction of the "fatty liver" cpAOPnet as well as chemical-specific aop

#setoutput directory:
outdir<-"../Output/NetRun1501/Paper"
indir<-"../Output/NetRun1501/Update"
#bring in the datasets...node lists
allEdges<-read.table(file.path(indir, "LiverTCAReactomePhenoChemGraph.txt"), sep='\t', header=TRUE, comment.char="" , strip.white=TRUE)
#labels
labels<-read.table(file.path(indir, "LiverTCAPhenoChemReactomeLabels.txt"), sep='\t', header=TRUE, comment.char="", quote="", strip.white=TRUE)

aop<-read.table("../Files/LiverPhenotypeLabelsV2.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE)
#remove kidney
#aop<-subset(aop, Source != "kidney")
####################
#this function will get out the neighbors for a graph
get_neighbors<-function(edgeList, target, bothNodes= FALSE){
  if(bothNodes==FALSE){
    tmp<-subset(edgeList, V1 %in% target | V2 %in% target); el2<-unique(c(as.character(tmp$V1), as.character(tmp$V2)))  
  }else{
    tmp<-subset(edgeList, V1 %in% target & V2 %in% target); el2<-unique(c(as.character(tmp$V1), as.character(tmp$V2)))
  }
  el2
}

#this adds in the human readable labels to the edges
temp<-merge(allEdges, unique(labels[,1:2]), by.x="V1", by.y="ID", all.x=TRUE)
colnames(temp)<-gsub("altID", "V1.name",colnames(temp))
temp2<-merge(temp, unique(labels[,1:2]), by.x="V2", by.y="ID", all.x=TRUE)
colnames(temp2)<-gsub("altID", "V2.name",colnames(temp2))
allEdges2<-temp2[,colnames(temp2)[c(2,14,1,15,3:13)]]
#write table
write.table(allEdges2, file=file.path(outdir,"LiverTCAReactomePhenoChemGraphv2.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )

#to generate the list of compounds that are connected to the outcomes of interest.
phen<-list(alteredMet=unique(subset(aop, Category =="altered metabolic activity")$Label), cholestasis=unique(subset(aop, Category =="cholestasis")$Label), lipidDys=unique(subset(aop, Category =="lipid dysregulation")$Label))
#chems<-subset(labels, Group=="Chem")
#paths<-subset(labels, Group=="Reactome")
#this is to get the first neightbor nodes for the different FL AO
FLpheno1<-get_neighbors(allEdges, target=phen[[1]])
FLpheno2<-get_neighbors(allEdges, target=phen[[2]])
FLpheno3<-get_neighbors(allEdges, target=phen[[3]])

#this is a list of all the nodes that connect to all 3 "aos" for fatty liver
FLpheno<-Reduce(intersect, list(FLpheno1, FLpheno2, FLpheno3))
#this is just the chemicals
FLChems<-FLpheno[FLpheno %in% as.character(subset(labels, Group=="Chem")$ID)]
FLChems
# [1] "acetaminophen"             "bendazac"                  "benzbromarone"            
# [4] "bortezomib"                "carbon_tetrachloride"      "chlormezanone"            
# [7] "colchicine"                "cycloheximide"             "dantrolene"               
# [10] "diazepam"                  "ethambutol"                "ethanol"                  
# [13] "ethinylestradiol"          "ethionamide"               "flutamide"                
# [16] "hexachlorobenzene"         "lomustine"                 "methapyrilene"            
# [19] "methylene_dianiline"       "monocrotaline"             "naphthyl_isothiocyanate"  
# [22] "naproxen"                  "nimesulide"                "nitrosodiethylamine"      
# [25] "papaverine"                "phenobarbital"             "phenylanthranilic_acid"   
# [28] "phorone"                   "promethazine"              "propylthiouracil"         
# [31] "puromycin_aminonucleoside" "rifampicin"                "tamoxifen"                
# [34] "terbinafine"               "thioacetamide"             "vitamin_A"                
# [37] "WY-14643"                           

#now to get out the carbon tet network : "carbon_tetrachloride" 

temp<-get_neighbors(allEdges2, target="carbon_tetrachloride")
temp2<-get_neighbors(allEdges2, target=temp)
#need to remove all the chemicals nodes except for the target one
temp3<-c(setdiff(temp2, subset(labels, Group=="Chem")$ID),"carbon_tetrachloride")
ccl4net<-subset(allEdges2, V1%in% temp3 & V2 %in% temp3)
write.table(ccl4net, file=file.path(outdir, "ccl4Net.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
#for a streamline network 
ccl4netSLIM<-unique(ccl4net[,c(2,4,5,6)])
write.table(ccl4netSLIM, file=file.path(outdir, "ccl4NetSLIM.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )

#######################################
#to get out the cpAOP network for lipid dysregulation
#use all the chemicals tied to the endpoints and look at the overlap

for(i in 1:length(FLChems)){
  #get the second neighbors of the seeding compound
  #first neighbors
  fls1<-get_neighbors(allEdges2, target=FLChems[i]);
  fls2<-get_neighbors(allEdges2, target=fls1)
  if(i==1){
    nodesFL<-fls2
  }else{
    nodesFL<-unique(intersect(nodesFL, fls2))
  }
}
#removing the chemicals that may be attached bc of shared phenotypes etc
nodesFL2<-setdiff(as.character(nodesFL), as.character(subset(labels, Group== "Chem")$ID))
#getting the network out, note that a lot of the atributes arent needed bc they just dont make sense in this context(the chemicals have been removed)
FLcpAOP<-subset(allEdges2, V1 %in% nodesFL2 & V2 %in% nodesFL2)[,c(1:8)]

#note that it is the "name" that you would want to map on bc that gets rid of some of the coding isses b/w datasets
write.table(FLcpAOP, file=file.path(outdir, "FLcpAOPchemisect.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
#to add in just the toxcast
tctemp<-subset(allEdges2, (V1 %in% nodesFL2 & TYPE =="TCHum2ReactRat") | (V2 %in% nodesFL2 & TYPE =="TCHum2ReactRat"))
#adding to the network
FLcpAOPtc<-rbind(FLcpAOP, tctemp[,1:8])
write.table(FLcpAOPtc, file=file.path(outdir,"FLcpAOPchemisectTC.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
