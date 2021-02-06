pair_find<-function(data=data,RAN_num=15,RAP_num=30,k="euclidean",rawoff=0.05, groupnum, filename){
  suppressMessages(library(tidyverse))
  suppressMessages(library(fdrtool))
  suppressMessages(library(qvalue))
  
  RAN<-data[1:RAN_num,]
  RAP<-data[(RAN_num+1):(RAN_num+RAP_num),]
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
  
  results_sample_pair<-matrix(nrow=0,ncol=(2*RAP_knn_num+2*RAP_knn_num+1))
  
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
        
        new_RANP<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAP)
        new_RANP_dis<-dist(new_RANP,method = k)
        new_RANP_dis<-as.matrix(new_RANP_dis)
        wa_RANP_dis<-sort(new_RANP_dis[,1])[2:(RAP_knn_num+1)]
        knn_RANP_mean<-mean(as.numeric(wa_RANP_dis))
        knn_RANP_sd<-sd(as.numeric(wa_RANP_dis))
        
        name_NN<-c(rep(0,RAN_knn_num))
        name_NP<-c(rep(0,RAP_knn_num))
        
        if(knn_RANP_mean <= (knn_RANN_mean+knn_RANP_sd)){
          for (m in 1:RAN_knn_num){
            RANN_sample<-rownames(new_RANN_dis)[which(new_RANN_dis[,1] == wa_RANN_dis[m])]
            name_NN[m]<-RANN_sample
          }
          for (h in 1:RAP_knn_num){
            RANP_sample<-rownames(new_RANP_dis)[which(new_RANP_dis[,1] == wa_RANP_dis[h])]
            name_NP[h]<-RANP_sample
          }
          new.row<-c(results_whole[j,1],rownames(data[which(rownames(data) == results_whole[j,1]),]), name_NN, wa_RANN_dis, name_NP, wa_RANP_dis)
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
          new.row<-c(results_whole[j,1],rownames(data[which(rownames(data) == results_whole[j,1]),]), name_PN, wa_RAPN_dis, name_PP, wa_RAPP_dis)
          #print(new.row)
          results_sample_pair<-rbind(results_sample_pair,new.row)
        }
      }
    }
  }
  res<-results_sample_pair[,1:(2*RAP_knn_num+2*RAP_knn_num+1)]
  res<-res[,c(2:(2+RAN_knn_num),(2+1+2*RAN_knn_num):(2+2*RAN_knn_num+RAP_knn_num))] %>%
    data.frame() %>% remove_rownames()
  #return(res)
  #sum(str_detect(res$X1,groupnum$Var1[1]%>% as.vector()))
  
  pairinfor = matrix(ncol = 2,nrow = 0)
  for (i in 1:dim(res)[1]) {
    if (i <= sum(str_detect(res$X1,groupnum$Var1[1]%>% as.vector()))) {
      for (j in (2+RAN_knn_num):(1+RAN_knn_num+RAP_knn_num)) {
        pairinfor<-rbind(pairinfor,c(res$X1[i]%>% as.vector(),res[i,j]%>% as.vector()))
      }
    }else{
      for (j in 2:(1+RAN_knn_num)) {
        pairinfor<-rbind(pairinfor,c(res[i,j]%>% as.vector(),res$X1[i]%>% as.vector()))
      }
    }
  }
  
  pairinfor <- pairinfor %>% data.frame() %>% rename(Ctl=1,Disease=2) %>% arrange(Ctl,Disease) %>% unique()
  write.csv(pairinfor,file = filename,row.names = F)
  cat(paste("the redundant pair number is ",dim(pairinfor)[1],"\n",sep = ''))
  
  res2<-matrix(ncol = 7,nrow = 0)
  for (i in 1:dim(data)[2]) {
    former<-rep(NA,dim(pairinfor)[1])
    latter<-rep(NA,dim(pairinfor)[1])
    for (j in 1:dim(pairinfor)[1]) {
      index_former<-which(rownames(data) == pairinfor$Ctl[j])
      former[j]=as.numeric(as.character(data[index_former,i]))
      index_latter<-which(rownames(data) == pairinfor$Disease[j])
      latter[j]=as.numeric(as.character(data[index_latter,i]))
    }
    diff<-former-latter
    if (length(diff[diff>0]) >length(diff[diff<0])) {
      enrich="Ctl"
    }else if(length(diff[diff>0]) <length(diff[diff<0])){
      enrich="Disease"
    }else{
      enrich="nodiff"
    }
    test<-wilcox.test(former,latter,paired = TRUE)
    #cat(test$p.value)
    if (test$p.value >= rawoff | is.na(test$p.value)) {
      enrich = "nodiff"
    }
    res2<-rbind(res2,c(as.character(colnames(data)[i]),test$p.value,mean(former),mean(latter),median(former),median(latter),enrich))
  }
  res2 <- res2 %>% data.frame() %>% rename(Species=1,pvalue=2,Ctlmean=3,Dismean=4,Ctlmedian=5,Dismedian=6,Enrichment=7) %>%
    mutate(pvalue = as.numeric(as.character(pvalue))) %>%
    arrange(pvalue) %>% na.omit() %>%mutate(pvalue = as.numeric(as.character(pvalue)))
  res2<-na.omit(res2)
  #res2$qvalue = qvalue(p = res2$pvalue)$qvalues
  res2$adj.fdr.p<-p.adjust(res2$pvalue,method="fdr",n=length(res2$pvalue))
  #res2$fdr.qval=fdrtool(res2$pvalue,statistic="pvalue")$qval
  #res2$adj.bonferroni.p <- p.adjust(res2$pvalue,method="bonferroni",n=length(res2$pvalue))
  return(res2)
}
