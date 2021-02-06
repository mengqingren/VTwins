### VTwins

* the R code is used to analyze high dimensional metagenomic data
* the function focused on enriched species or genus for phenotypes
* can be run on windows or linux

##### Install package dependency (Packages_version)

* tidyverse_1.3.0
* readxl_1.3.1
* qvalue_2.16.0
* fdrtool_1.2.15 

##### Usage

    source(file = "Pair_Find.R")
    pair_find(data=data,RAN_num=15,RAP_num=30,k="euclidean",rawoff=0.05)

##### Examples

    #source function
    source(file = "Pair_Find.R")
    
    #1. mock dataset
    Mock <- matrix(runif(8000),ncol = 100,nrow = 80) %>% data.frame() # must be a dataframe
    colnames(Mock) <- paste("Feature",1:100,sep = "")
    rownames(Mock) <- c(paste("Ctl",1:39,sep = ""),paste("Disease",1:41,sep = ""))
    
    #2. Metadata
    MetaData <- data.frame(SampleID = rownames(Mock),Group = rep(c("Ctl","Disease"),c(39,41)))
    
    filename = "Ctl-Disease-Pair.csv" # must provide a filename to write the output file in current directory
    groupnum = table(MetaData$Group) %>% data.frame() # phenotype samples count and must be corresponding to RAN(Control) and RAP(Disease), repsectively.
    #we should always keep the ctl ahead of disease manually
    res<-pair_find(data=Mock,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
    # the redundant pair number is 254
    write.csv(res,file = "Feature-PairWilcoxonRank.Res.csv",row.names = F)

