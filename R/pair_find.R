#' @title Perform phenotype-enriched features analysis
#'
#' @description Perform phenotype-enriched features analysis
#' @param data
#' @param phenodata
#' @param k
#' @param SavePath
#' @param ShuffleWstat
#' @param BoundarySample
#' @param BoundaryPair
#' @param BoundaryPair
#' @param ShuffleTime
#' @param Uppercent
#' @param DownPercent
#' @return Mid.Matrix
#' @examples
#' @export

pair_find<-function(data=data,phenodata=data.frame(),k="euclidean",SavePath = NULL,ShuffleWstat = NULL, BoundarySample = NULL,BoundaryPair=NULL,PvalueCutoff=0.05,ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8){ # colnames(phenodata) = c("id","grp"); grp1 = "Ctrl", grp2 = "Disease"
  suppressMessages(library(tidyverse))
  #suppressMessages(library(fdrtool))
  #suppressMessages(library(qvalue))
  phenodata <- phenodata %>% dplyr::arrange(grp) %>% data.frame()
  #data must be ordered as grp order
  groupnum <- table(phenodata$grp) %>% data.frame()
  RAN_num = groupnum[1,2]
  RAP_num = groupnum[2,2]
  RAN.Samples<-phenodata %>% dplyr::filter(grp=="grp1")
  RAP.Samples<-phenodata %>% dplyr::filter(grp=="grp2")
  RAN <- data[rownames(data) %in% RAN.Samples$id,] %>% data.frame()
  RAP <- data[rownames(data) %in% RAP.Samples$id,] %>% data.frame()
  n=dim(data)[1]
  num=floor(sqrt(RAN_num+RAP_num))
  results_whole=matrix(nrow=0,ncol=(num+1))
  for (i in 1:n) {
    new<-rbind(data[i,],data)
    new_dis<-dist(new,method = k)
    new_dis<-as.matrix(new_dis)
    wa_d<-sort(new_dis[,1])[3:(num+2)]
    new_row<-c(rownames(data)[i],wa_d)
    results_whole<-rbind(results_whole,new_row)
  }
  mean_whole=mean(as.numeric(as.vector(results_whole[,2:(num+1)])))
  sd_whole=sd(as.numeric(as.vector(results_whole[,2:(num+1)])))

  RAP_knn_num=floor(sqrt(RAP_num))
  RAN_knn_num=floor(sqrt(RAN_num))

  results_sample_pair<-matrix(nrow=0,ncol=(2*RAN_knn_num+2*RAP_knn_num+2))

  for (j in 1:(RAN_num+RAP_num)){
    sample_mean<-mean(as.numeric(results_whole[j,2:(num+1)]))
    if (sample_mean <= (mean_whole + sd_whole)){
      #print(results_whole[j,1])
      if (results_whole[j,1] %in% rownames(RAN)){
        new_RANN<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAN)
        #print(results_whole[j,1])
        new_RANN_dis<-dist(new_RANN,method = k)
        new_RANN_dis<-as.matrix(new_RANN_dis)
        wa_RANN_dis<-sort(new_RANN_dis[,1])[3:(RAN_knn_num+2)]
        knn_RANN_mean<-mean(as.numeric(wa_RANN_dis))
        knn_RANN_sd <- sd(as.numeric(wa_RANN_dis))

        new_RANP<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAP)
        new_RANP_dis<-dist(new_RANP,method = k)
        new_RANP_dis<-as.matrix(new_RANP_dis)
        wa_RANP_dis<-sort(new_RANP_dis[,1])[2:(RAP_knn_num+1)]
        knn_RANP_mean<-mean(as.numeric(wa_RANP_dis))
        knn_RANP_sd<-sd(as.numeric(wa_RANP_dis))

        name_NN<-c(rep(0,RAN_knn_num))
        name_NP<-c(rep(0,RAP_knn_num))

        if(knn_RANP_mean <= (knn_RANN_mean+knn_RANN_sd)){
          for (m in 1:RAN_knn_num){
            RANN_sample<-rownames(new_RANN_dis)[which(new_RANN_dis[,1] == wa_RANN_dis[m])]
            name_NN[m]<-RANN_sample
          }
          for (h in 1:RAP_knn_num){
            RANP_sample<-rownames(new_RANP_dis)[which(new_RANP_dis[,1] == wa_RANP_dis[h])]
            name_NP[h]<-RANP_sample
          }
          new.row<-c(as.character(results_whole[j,1]),rownames(data[which(rownames(data) == results_whole[j,1]),]),name_NN, as.vector(wa_RANN_dis), name_NP, as.vector(wa_RANP_dis))
          #print(new.row)
          results_sample_pair<-rbind(results_sample_pair,new.row)
        }
      } else {
        new_RAPN<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAN)
        new_RAPN_dis<-dist(new_RAPN,method = k)
        new_RAPN_dis<-as.matrix(new_RAPN_dis)
        wa_RAPN_dis<-sort(new_RAPN_dis[,1])[2:(RAN_knn_num+1)]
        knn_RAPN_mean<-mean(wa_RAPN_dis)

        new_RAPP<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAP)
        new_RAPP_dis<-dist(new_RAPP,method = k)
        new_RAPP_dis<-as.matrix(new_RAPP_dis)
        wa_RAPP_dis<-sort(new_RAPP_dis[,1])[3:(RAP_knn_num+2)]
        knn_RAPP_mean<-mean(wa_RAPP_dis)
        knn_RAPP_sd<-sd(wa_RAPP_dis)

        name_PN<-c(rep(0,RAN_knn_num))
        name_PP<-c(rep(0,RAP_knn_num))

        if (knn_RAPN_mean <= (knn_RAPP_mean+knn_RAPP_sd)){
          for (m in 1:RAN_knn_num){
            RAPN_sample<-rownames(new_RAPN_dis)[which(new_RAPN_dis[,1] == wa_RAPN_dis[m])]
            name_PN[m]<-RAPN_sample
          }
          for (h in 1:RAP_knn_num){
            RAPP_sample<-rownames(new_RAPP_dis)[which(new_RAPP_dis[,1] == wa_RAPP_dis[h])]
            name_PP[h]<-RAPP_sample
          }
          new.row<-c(as.character(results_whole[j,1]),rownames(data[which(rownames(data) == results_whole[j,1]),]), name_PN, as.vector(wa_RAPN_dis), name_PP, as.vector(wa_RAPP_dis))
          #print(new.row)
          results_sample_pair<-rbind(results_sample_pair,new.row)
        }
      }
    }
  }
  cat("KNN FINISHED\n")

  if (dim(results_sample_pair)[1] > 1) {
    res<-results_sample_pair[,1:(2*RAN_knn_num+2*RAP_knn_num+2)]
  }else{
    cat("Nothing Done \n")
    return(NULL)
  }

  if (is.null(SavePath)) {
    SavePath="./"
  }else{
    if(!dir.exists(SavePath)){
      dir.create(SavePath)
    }
  }

  #res<-res[,c(2:(2+RAN_knn_num),(2+1+2*RAN_knn_num):(2+2*RAN_knn_num+RAP_knn_num))] %>% data.frame() %>% remove_rownames()
  res<-res[,c(2:(2+2*RAN_knn_num+2*RAP_knn_num))] %>% data.frame() %>% remove_rownames()
  if (!is.null(BoundarySample)) {
    write.csv(res,paste(SavePath,"/",BoundarySample,".csv",sep = ''),row.names = F)
  }
  #View(res)
  #return(res)
  #sum(str_detect(res$X1,groupnum$Var1[1]%>% as.vector()))
  #View(res)
  #write.csv(res,file = filename)
  pairinfor = matrix(ncol =3,nrow = 0)
  for (i in 1:dim(res)[1]) {
    indexpo1 = which(phenodata[,1] == res[i,1])
    indexpo2 = which(phenodata[,1] == res[i,2])
    if (phenodata[indexpo1,2] == phenodata[indexpo2,2]) {
      for (mkl in 1:RAP_knn_num) {
        pairinfor <- rbind(pairinfor,c(res[i,1],res[i,1+2*RAN_knn_num+mkl],res[i,1+2*RAN_knn_num+mkl+RAP_knn_num]))
      }
    }else{
      for (pkl in 1:RAN_knn_num) {
        pairinfor <- rbind(pairinfor,c(res[i,pkl+1],res[i,1],res[i,pkl+1+RAN_knn_num]))
      }
    }
  }

  Newmatrix <<- matrix(ncol = 3,nrow = 0)
  Extract_Dist <- function(PairData,...){
    if (dim(PairData)[1] > 0) {
      Findline <- PairData[which(PairData[,3] == min(PairData[,3])),]
      Newmatrix <<- rbind(Newmatrix,Findline)
      #View(Newmatrix)
      NewMidData <- PairData %>% filter(!(Ctl %in% Findline$Ctl) & !(Disease %in% Findline$Disease))
      #cat(dim(NewMidData)[1])
      #cat("\n")
      return(Extract_Dist(NewMidData))
    }else{
      return(Newmatrix)
    }
  }

  #View(pairinfor)
  pairinfor <- pairinfor %>% data.frame() %>% dplyr::rename(Ctl=1,Disease=2,Distance=3) %>%
    mutate(Distance = as.numeric(as.character(Distance))) %>% dplyr::arrange(Ctl,Disease) %>% unique()
  #View(pairinfor)
  pairinfor <- Extract_Dist(pairinfor)
  pairinfor <- pairinfor %>% data.frame() %>% dplyr::select(-Distance)

  if (!is.null(BoundaryPair)) {
    write.csv(pairinfor,paste(SavePath,"/",BoundaryPair,".csv",sep = ''),row.names = F)
  }
  #write.csv(pairinfor,file = "test.csv",row.names = F)
  cat(paste("the redundant pair number is ",dim(pairinfor)[1],"\n",sep = ''))
  cat("PAIR FINISHED\n")
  #return(pairinfor)

  ###### Pair shuffle

  ShufflePair <- lapply(1:ShuffleTime, function(Shuffle){
    Random <- runif(1,DownPercent,Uppercent)
    n <- round(dim(pairinfor)[1]*Random,0)
    RandomIndex <- sample(1:dim(pairinfor)[1],n,replace = F)

    pairinforMiddata.mid1.1 <- pairinfor[-RandomIndex,] %>% data.frame()
    pairinforMiddata.mid2 <- pairinfor[RandomIndex,] %>% data.frame()
    pairinforMiddata.mid2.1 <- data.frame(Ctl = pairinforMiddata.mid2$Disease, Disease = pairinforMiddata.mid2$Ctl)
    pairinforMiddata.mid <- rbind(pairinforMiddata.mid1.1,pairinforMiddata.mid2.1) %>% data.frame()
  })

  cat("SHUFFLE PAIR FINISHED\n")

  ####### calculate W stat and rank
  FinalMatrix <- matrix(ncol = 2+ShuffleTime,nrow = 0) #Species, Original W,
  MeanData <- matrix(ncol = 3,nrow = 0) # record mean value of pair
  for (i in 1:dim(data)[2]) {
    ### original cohort
    former<-rep(NA,dim(pairinfor)[1])
    latter<-rep(NA,dim(pairinfor)[1])
    for (j in 1:dim(pairinfor)[1]) {
      index_former<-which(rownames(data) == pairinfor$Ctl[j])
      former[j]=as.numeric(as.character(data[index_former,i]))
      index_latter<-which(rownames(data) == pairinfor$Disease[j])
      latter[j]=as.numeric(as.character(data[index_latter,i]))
    }
    Middata <- data.frame(Ctrl = former, CRC = latter)
    test<-wilcox.test(Middata$Ctrl,Middata$CRC,paired = TRUE)
    AllNW <- (dim(pairinfor)[1]*(dim(pairinfor)[1]+1))/2
    OriginalStat = AllNW - test$statistic

    Ctrlmean <- mean(Middata$Ctrl)
    CRCmean <- mean(Middata$CRC)
    MeanData <- rbind(MeanData,c(as.character(colnames(data)[i]),Ctrlmean,CRCmean))

    cat(as.character(colnames(data)[i]))
    cat("\n")
    #print(AllNW)

    ### calculate shuffle W stat
    #ShuffleStat <- rep(NA,ShuffleTime)

    ShuffleStat <- lapply(ShufflePair, function(x){
      former<-rep(NA,dim(x)[1])
      latter<-rep(NA,dim(x)[1])
      for (j in 1:dim(x)[1]) {
        index_former<-which(rownames(data) == x$Ctl[j])
        former[j]=as.numeric(as.character(data[index_former,i]))
        index_latter<-which(rownames(data) == x$Disease[j])
        latter[j]=as.numeric(as.character(data[index_latter,i]))
      }

      Middata.mid <- data.frame(Ctrl = former, CRC = latter)
      test1<-wilcox.test(Middata.mid$Ctrl,Middata.mid$CRC,paired = TRUE)
      #ShuffleStat[Iter] = AllNW - test1$statistic
      return(AllNW - test1$statistic)
    })
    ShuffleStat <- unlist(ShuffleStat)
    FinalMatrix <- rbind(FinalMatrix,c(as.character(colnames(data)[i]),OriginalStat,ShuffleStat))
  }
  colnames(FinalMatrix) <- c("Feature","OriginalSata.W",paste("ShuffleStat.W.",1:ShuffleTime,sep = ""))

  if (!is.null(ShuffleWstat)) {
    write.csv(pairinfor,paste(SavePath,"/",ShuffleWstat,".csv",sep = ''),row.names = F)
  }

  Mid.Matrix <- matrix(ncol = 6,nrow = 0)
  for (I.index in 1:dim(FinalMatrix)[1]) {
    Increasing.Rank.Min <- rank(as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "min")[1]
    Decreasing.Rank.Min <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "min")[1]
    Increasing.Rank.Max <- rank(as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "max")[1]
    Decreasing.Rank.Max <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "max")[1]
    Increasing.Rank.Average <- rank(as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "average")[1]
    Decreasing.Rank.Average <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "average")[1]
    Mid.Matrix <- rbind(Mid.Matrix,c(Increasing.Rank.Min,Decreasing.Rank.Min,Increasing.Rank.Max,Decreasing.Rank.Max,Increasing.Rank.Average,Decreasing.Rank.Average))
  }
  Mid.Matrix <- Mid.Matrix %>% data.frame() %>% remove_rownames() %>%
    dplyr::rename(Increasing.Rank.Min=1,Decreasing.Rank.Min=2,Increasing.Rank.Max=3,Decreasing.Rank.Max=4,Increasing.Rank.Average=5,Decreasing.Rank.Average=6)

  Mid.Matrix$Decre.maxRank.P <- Mid.Matrix$Decreasing.Rank.Max/(ShuffleTime+1)
  Mid.Matrix$Decre.aveRank.P <- Mid.Matrix$Decreasing.Rank.Average/(ShuffleTime+1)
  Mid.Matrix$Decre.minRank.P <- Mid.Matrix$Decreasing.Rank.Min/(ShuffleTime+1)
  Mid.Matrix$Incre.maxRank.P <- Mid.Matrix$Increasing.Rank.Max/(ShuffleTime+1)
  Mid.Matrix$Incre.aveRank.P <- Mid.Matrix$Increasing.Rank.Average/(ShuffleTime+1)
  Mid.Matrix$Incre.minRank.P <- Mid.Matrix$Increasing.Rank.Min/(ShuffleTime+1)
  Mid.Matrix$Species <- FinalMatrix[,1] #Feature

  Mid.Matrix <- Mid.Matrix %>% data.frame() %>%
    dplyr::arrange(Decre.minRank.P) %>%
    dplyr::mutate(Decre.minRank.P.FDR = p.adjust(.$Decre.minRank.P,method = "BH",n=length(.$Decre.minRank.P))) %>%
    dplyr::arrange(Decre.maxRank.P) %>%
    dplyr::mutate(Decre.maxRank.P.FDR = p.adjust(.$Decre.maxRank.P,method = "BH",n=length(.$Decre.maxRank.P))) %>%
    dplyr::arrange(Decre.aveRank.P) %>%
    dplyr::mutate(Decre.aveRank.P.FDR = p.adjust(.$Decre.aveRank.P,method = "BH",n=length(.$Decre.aveRank.P))) %>%
    dplyr::arrange(Incre.aveRank.P) %>%
    dplyr::mutate(Incre.aveRank.P.FDR = p.adjust(.$Incre.aveRank.P,method = "BH",n=length(.$Incre.aveRank.P)))

  MeanData <- MeanData %>% data.frame() %>% dplyr::rename(Species=1,Ctlmean=2,Dismean=3)
  Mid.Matrix <- merge(Mid.Matrix,MeanData,by="Species") %>% data.frame() %>% mutate(Enrieched = if_else(Decre.aveRank.P <= PvalueCutoff,"Disease",if_else(Incre.aveRank.P <= PvalueCutoff,"Ctrl","N.S.")))

  cat("All done\n")
  return(Mid.Matrix)
}# END - function: Pair_Find

