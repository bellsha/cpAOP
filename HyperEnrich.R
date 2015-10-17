#general (internal) group enrichment function
#Ref_DF is a dataframe with col1="label" and col2="feature"
#samp contains just features 
#background is the background list of features, ie those present in the dataset

HyperEnrich <- function(samp, background,Ref_DF, Reference, statsB=NULL,Nannot=2,m=m, n=n, ...) {
  #m = the # of features associated with the catagories for Reference
  if(is.null(statsB) | is.null(m)){
    statsB <- data.frame(Size=sapply(Reference, length))
    statsB <- data.frame(Label=row.names(statsB), statsB)
    row.names(statsB)<-1:nrow(statsB)
    m <- as.vector(statsB$Size)
  }
  
  #x = # of features in sample set for each Label in background/dataset
  #this will get all the ones for the node counted for how many times they appear in the sample
  #note that it will return a count of 1 greater than the number of samples, so subtract 1 from x
  statsS <- as.vector(sapply(Reference, function(f) sum(samp %in% unlist(f))))
  x <- statsS       
  
  ## (n): Obtain the number of total number of features with a Label in dataset
  if(is.null(n)){
    if(is.null(background)){
      n <- length(unique(Ref_DF[, 2]))  
    }else{
      #note this is to adjust the total numnber of "balls" according to the # of annotations
      #it is expected that the supplied input 
      n <- length(background[background %in% unique(unique(Ref_DF[, 2]))])
    }
  }

  ## (k): Obtain number of features in the sample (having a Label)
  #to accomidate for duplicates, provide background
  if(is.null(background)){
    k <- length(unique(Ref_DF[Ref_DF[,2] %in% samp, 2]))	
  } else{
    k <- length(background[background %in% unique(Ref_DF[Ref_DF[,2] %in% samp, 2])])
  }
  ###############
  ## Obtain features names matching Labels
  #match_key <-sapply(Reference, function(f) unlist(f)[unlist(f) %in% as.character(samp)])
  match_key <-sapply(Reference, function(f) f[unlist(f) %in% as.character(samp)])
  match_key <- sapply(match_key, function(f) { paste(f, collapse=" ") } )
  match_key <- as.vector(match_key)
  key <- match_key; key[key==""] <- "NA"
  ###

  ####################3  
  ## Apply phyper function (x=x-1 because we want P[X>= x]; ie over representation)
  phyp_v <- phyper(x-1, m, n-m , k, lower.tail = FALSE)
  
  ## P-value correction according to Bioinformatics, 20, 3710-3715
  Ncorrect <- table(Ref_DF[Ref_DF[,2] %in% samp, 1]) # Obtain the Labels with direct annotations from sample set 
  Ncorrect <- sum(Ncorrect >= Nannot) # Count only those that have 2 or more annotations from sample set
  if(Ncorrect<=1) { 
    adj_phyp_v <- phyp_v # no adjustment necessary if Ncorrect <= 1
  } else {
    adj_phyp_v <- phyp_v * Ncorrect # Calculates simple Bonferroni correction. 
    adj_phyp_v[adj_phyp_v >= 1] <- 1
  }
  
  ## Generate output data format
  result_df <- data.frame(NumAnn=k, statsB, SampleMatch=x, Phyper=phyp_v, Padj=adj_phyp_v, SampleKeys=key )
  result_df <- result_df[order(result_df$Phyper), ]
  result_df
}

#data is a dataframe with first column=clusterID and second column=featureID and thrid column= size
ClusReport <- function(data, background=NULL,Reference=NULL,RefData=NULL, MinAnn=2,MinSamp=2, pMax=0.05){
  if(is.null(RefData)){
    print("Must supply a reference dataframe")
    break
  }
  #will not work right if duplicates in RefData
  RefData<-unique(RefData)
  if(!is.null(background)){
    
    Ref_DF<-merge(RefData, as.data.frame(background), by.x=colnames(RefData)[2], by.y="background")
    Ref_DF<-unique(Ref_DF[,c(2,1)])
  }else{
    Ref_DF<-RefData
  }
  colnames(Ref_DF)<-c('label','feature')
  if(is.null(Reference)){
    lab<-unique(as.character(Ref_DF[,1]))
    Labs<-list()
    for(d in 1:length(lab)){
      #  print(lab[d])
      Labs[[lab[d]]]<-subset(Ref_DF, Ref_DF[,1]==lab[d])[,2]
    }
    Reference<-Labs
  }
  #just so they are not being calculated each itteration:
  statsB <- data.frame(Size=sapply(Reference, length))
  statsB <- data.frame(Label=row.names(statsB), statsB)
  row.names(statsB)<-1:nrow(statsB)
  m <- as.vector(statsB$Size)
  if(is.null(background)){
    n <- length(unique(Ref_DF[, 2]))  
  }else{
    #note this is to adjust the total numnber of "balls" according to the # of annotations
    #it is expected that the supplied input 
    n <- length(background[background %in% unique(Ref_DF[, 2])])
  }
  ###
  #note that the Reference should be matched to the sample for correct accounting, ie thi
  clus<-unique(as.vector(data[,1]))
  containerDF<-NULL
  for(i in clus){
    #cat("Processing group", i, "\n")
    test_sample<-as.vector(data[data[,1]==i,2])
    test_result<-HyperEnrich(samp=test_sample, background=background, Reference=Reference, Ref_DF=Ref_DF,MinAnn=MinAnn, statsB=statsB, m=m, n=n)
    tempDF<-data.frame(CLID=rep(i, times=nrow(test_result)), CLSZ=rep(unique(as.vector(data[data[,1]==i,3])), times=nrow(test_result)), test_result)
    containerDF <- rbind(containerDF, tempDF)
    
  }
  results<-subset(containerDF, Padj <= pMax & SampleMatch >= MinSamp & Size >=MinAnn)
  return(results)
}
