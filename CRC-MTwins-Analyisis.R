setwd("D:/CRC-Pair")
#library(curatedMetagenomicData)
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggstar)
library(fdrtool)
library(qvalue)
library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(ROCR)
library("caTools")
library(sampling)
library(scales)
library(ggpmisc)
library(ggrepel)
library(UpSetR)
library(gridExtra)
library(gridGraphics)
library(pheatmap)
library(ggthemes)
library(ggsci)
library(ComplexHeatmap)
library(ggrepel)
library(randomcoloR)
##########################--------- NEW 8 STUDY ---------########################
####### Analysis Function ########
pair_find<-function(data=data,RAN_num=15,RAP_num=30,k="euclidean",rawoff=0.05){
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
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
####### Merge Data with Speices Abundance and Metadata ########

#### @0 Download Data ####
Study <- c("FengQ_2015","HanniganGD_2017","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014")
suffix<-c(".metaphlan_bugs_list.stool",".pathabundance_relab.stool") #gene_family Data
#metadata <- c("subjectID","study_condition","disease","disease_subtype","country","number_bases","body_site","age","age_category","gender","BMI","")

for (i in Study) {
  for (j in suffix) {
    datanew <- eval(parse(text = paste(i,j,"()",sep = "")))
    datanew@phenoData@data %>% write.csv(file = paste("Data/",i,".metadata.csv",sep = ""),quote = F,row.names = TRUE,col.names = TRUE)
    datanew@assayData$exprs %>% write.table(file = paste("Data/",i,j,".txt",sep = ""),quote = F,row.names = TRUE,col.names = TRUE,sep = '\t')
  }
}

#### @1 Merge MetaData ####
#list.files("Old/",pattern = "[.txt]")
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

MetaData <- data.frame()
for (name in Study) {
  metadata<-read.csv(paste("Data/",name,".metadata.csv",sep = ''),row.names = 1,check.names = F)
  MetaData<-metadata %>%  rownames_to_column() %>% filter(study_condition == "control" | study_condition == "CRC") %>% 
    mutate(study_condition = factor(study_condition,levels = c("control","CRC"))) %>%
    arrange(study_condition) %>% select(rowname,study_condition) %>% mutate(Study=name) %>% 
    rbind(MetaData)
}

MetaData[(MetaData$Study != "PRJEB27928" & !str_detect(MetaData$rowname,"\\.21\\.0")),] %>% dim()
write.csv(MetaData,"Species-8Study-20201010/EightStudies-MetaData.csv",row.names = F)
write.csv(MetaData,"Pathway-8Study-20201010/EightStudies-MetaData.csv",row.names = F)

#### @2 Merge Species Abundance ####
SpeciesData <- read_excel("RawListSpecies.xlsx")
SpeciesData <- SpeciesData[str_detect(SpeciesData$Species,"k__Bacteria") & str_detect(SpeciesData$Species,"s__") & !str_detect(SpeciesData$Species,"t__") ,] %>% data.frame()
SpeciesData$Species <- SpeciesData$Species %>% as.character() %>% str_replace_all(".*s__","")

for (name in Study) {
  print(name)
  Speciesdata <- read.table(paste("Data/",name,".metaphlan_bugs_list.stool.txt",sep = ''),row.names = 1,check.names = F,sep = '\t',header = T)
  print(dim(Speciesdata))
  Speciesdata <- Speciesdata[str_detect(rownames(Speciesdata),"k__Bacteria") & str_detect(rownames(Speciesdata),"s__") & !str_detect(rownames(Speciesdata),"t__") ,] %>% data.frame()
  rownames(Speciesdata) = rownames(Speciesdata) %>% as.character() %>% str_replace_all(".*s__","")
  Speciesdata <- Speciesdata %>% rownames_to_column()
  SpeciesData <-  merge(SpeciesData,Speciesdata,by.x = "Species",by.y = "rowname",all = T)
  print(dim(SpeciesData))
}

SpeciesData[is.na(SpeciesData)] = 0

write.table(SpeciesData,"Species-8Study-20201010/EightStudies-SpeciesAbundance.txt",row.names = F,sep = '\t',quote = F)

#### @3 Merge Pathway Abundance ####
PathwayData <- read_excel("RawListPathway.xlsx")
PathwayData$Pathway <- PathwayData$Pathway %>% str_replace_all("-",".") #%>% str_replace_all(" ",".") %>% str_replace_all(":",".") %>% str_replace_all("\\(",".") %>% str_replace_all("\\)",".")
for (name in Study) {
  print(name)
  if (name == "PRJEB27928") {
    Pathwaydata <- read.csv(paste("Data/",name,".pathabundance_relab.csv",sep = ''),check.names = F,row.names = 1)
  }else{
    Pathwaydata <- read.table(paste("Data/",name,".pathabundance_relab.stool.txt",sep = ''),check.names = F,sep = '\t',header = T,row.names = 1)
  }
  
  Pathwaydata <- Pathwaydata[!str_detect(rownames(Pathwaydata),"\\|"),] %>% data.frame()
  Pathwaydata <- Pathwaydata %>% data.frame() %>% rownames_to_column()
  Pathwaydata$rowname <- Pathwaydata$rowname %>% str_replace_all("-",".")
  PathwayData <- merge(PathwayData,Pathwaydata,by.x = "Pathway",by.y = "rowname",all = T)
}

PathwayData[is.na(PathwayData)] = 0

write.table(PathwayData,"Pathway-8Study-20201010/EightStudies-PathwayAbundance.txt",row.names = F,sep = '\t')

################# Species Differential Analysis ###############
SpeciesData <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance.txt",header = T,row.names = 1,sep = '\t')
SpeciesData <- SpeciesData %>% t() %>% data.frame() %>% rownames_to_column()

MetaData <- read.csv("Species-8Study-20201010/EightStudies-MetaData.csv")
MetaData$rowname <- MetaData$rowname %>% str_replace_all("-",".")
SpeciesData2 <- merge(MetaData,SpeciesData,by = "rowname",all.x = T) %>% data.frame() %>% column_to_rownames("rowname")
SpeciesData2 <- as.matrix(SpeciesData2)
SpeciesData2[is.na(SpeciesData2)] = 0
write.table(SpeciesData2,file = "Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",quote = F,sep = "\t")
#### All Paired + Wilcoxon ####
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)

groupnum <- table(SpeciesData2$study_condition) %>% data.frame()
middata1 <- SpeciesData2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
filename=paste("Species-8Study-20201010/All-8Study-Contained.Speices.pair.csv",sep = '')
res<-pair_find(data=middata1,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
#ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
write.csv(res,file = paste("Species-8Study-20201010/All-8Study-Contained-Species-pair-wilcoxonsign-res.csv",sep = ''),row.names = F)


rownames(middata1)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
wilcox_res<-matrix(ncol = 7,nrow = 0)
## wilcox test differential species
for (j in 1:dim(middata1)[2]) {
  test<-wilcox.test(middata1[,j][1:groupnum$Freq[1]],middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]])
  wilcox_res <- rbind(wilcox_res,c(colnames(middata1)[j],test$p.value,mean(middata1[,j][1:groupnum$Freq[1]]),median(middata1[,j][1:groupnum$Freq[1]]),mean(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),median(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),j))
}
wilcox_res <- as.data.frame(wilcox_res)
colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
wilcox_res <- na.omit(wilcox_res)
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
  mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()

wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
write.csv(wilcox_res,file = paste("Species-8Study-20201010/All-8Study-Contained-Species-wilcoxon-test.csv",sep = ''),row.names = F)

#### Every Study Pair + Wilcoxon ####
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (name in Study) {
  middata2 <- SpeciesData2 %>% filter(Study == name) %>% arrange(study_condition)
  groupnum <- table(middata2$study_condition) %>% data.frame()
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  
  rownames(middata) <- c(paste(groupnum$Var1[1]%>% as.vector(),1:groupnum$Freq[1]),paste(groupnum$Var1[2]%>% as.vector(),1:groupnum$Freq[2])) %>% str_remove_all(" ")
  write.csv(data.frame(NewName=rownames(middata),subjectid=rownames(middata2)),file = paste(name,"-metaphlan2.SubjectID-Newname.csv",sep = ''),row.names = F)
  
  
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential species
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  write.csv(wilcox_res,file = paste(name,"-Metaphlan2-wilcoxonTest.csv",sep = ''),row.names = F)
  assign(paste(name,".wilcoxon",sep = ''),wilcox_res)
  
  # pair analysi wilcoxon pair test
  filename=paste(name,".metaphaln2.pair.csv",sep = '')
  res<-pair_find(data=middata,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
  assign(paste(name,".wilcoxonsign",sep = ''),res)
  write.csv(res,file = paste(name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),row.names = F)
}


######## Study Transfer Study Top50 ########
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
## perform RF model
## Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()

for (name in Study) {
  ## Process Data
  middata <- get(paste(name,".wilcoxonsign",sep = ''))
  middataTop50 <- middata %>% arrange(adj.fdr.p) %>% top_n(-50,adj.fdr.p)
  write.csv(middataTop50,file = paste(name,"-PairTop50.Species.csv",sep = ''),row.names = F)
    
  SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop50$Species)) 
  ## train test data
  TrainData <- SpeciesSelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
  TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
  ## model for self
  set.seed(123)
  split = sample.split(TrainData$Label,SplitRatio = .7)
  train_data = subset(TrainData,split == TRUE)
  test_data  = subset(TrainData,split == FALSE)
    
  TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
  control <- trainControl(method="repeatedcv",number=3,repeats=5)
    
  train_data$Label <- factor(train_data$Label,levels = c(0,1))
  test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
  fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
  rf.pred <- predict(fit.rf, TrainData)
  cm<-confusionMatrix(rf.pred,TrainData$Label)
  confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,BaseModel = name) %>%rbind(confusionMatrixdata)
    
  predob = predict(fit.rf,TrainData,type = "prob")
  pred<-prediction(predob[,2],TrainData$Label)
  perf<-performance(pred,'tpr','fpr')
    
  #extrac plot ROC data
  modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = name,BaseModel = name) %>%rbind(modelROC)
    
  auc<-performance(pred,"auc")
  auc<-unlist(slot(auc,"y.values"))
  model.AUC <- data.frame(Predict=name,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
    
  ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
    mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name,BaseModel = name) %>% rbind(ModelImportance)
    
  ## model for Other All Studies
  for (name2 in Study) {
    if (name != name2) {
      cat(paste(name,"\t",name2,"\n"))
      
      test.data <- SpeciesSelectData %>% filter(Study == name2) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
        
      control <- trainControl(method="repeatedcv",number=3,repeats=5)
      
      test.data$Label <- factor(test.data$Label,levels = c(0,1))
        
      rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
        
      rf.pred <- predict(rf.fit, test.data)
      cm<-confusionMatrix(rf.pred,test.data$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name2,BaseModel = name) %>%rbind(confusionMatrixdata)
        
      predob = predict(rf.fit,test.data,type = "prob")
      pred<-prediction(predob[,2],test.data$Label)
      perf<-performance(pred,'tpr','fpr')
        
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name2,BaseModel = name) %>%
        rbind(modelROC)
        
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(Predict=name2,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
        
      ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name2,BaseModel = name) %>% rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("Study2Study.Top50.",name,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("Study2Study.Top50.",name,".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("Study2Study.Top50.",name,".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = "Study2Study.Top50.RF.model.Importance.csv",row.names = F)

######## Study Transfer Study Top30 ########
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
## perform RF model
## Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()

for (name in Study) {
  ## Process Data
  middata <- get(paste(name,".wilcoxonsign",sep = ''))
  middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-30,adj.fdr.p)
  #write.csv(middataTop50,file = paste(name,"-PairTop50.Species.csv",sep = ''),row.names = F)
  
  SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
  ## train test data
  TrainData <- SpeciesSelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
    select(-study_condition,-Study)
  TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
  
  ## model for self
  set.seed(123)
  split = sample.split(TrainData$Label,SplitRatio = .7)
  train_data = subset(TrainData,split == TRUE)
  test_data  = subset(TrainData,split == FALSE)
  
  TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
  control <- trainControl(method="repeatedcv",number=3,repeats=5)
  
  train_data$Label <- factor(train_data$Label,levels = c(0,1))
  test_data$Label <- factor(test_data$Label,levels = c(0,1))
  
  fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
  
  rf.pred <- predict(fit.rf, TrainData)
  cm<-confusionMatrix(rf.pred,TrainData$Label)
  confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,BaseModel = name) %>%rbind(confusionMatrixdata)
  
  predob = predict(fit.rf,TrainData,type = "prob")
  pred<-prediction(predob[,2],TrainData$Label)
  perf<-performance(pred,'tpr','fpr')
  
  #extrac plot ROC data
  modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = name,BaseModel = name) %>%rbind(modelROC)
  
  auc<-performance(pred,"auc")
  auc<-unlist(slot(auc,"y.values"))
  model.AUC <- data.frame(Predict=name,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
  
  ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
    mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name,BaseModel = name) %>% rbind(ModelImportance)
  
  ## model for Other All Studies
  for (name2 in Study) {
    if (name != name2) {
      cat(paste(name,"\t",name2,"\n"))
      
      test.data <- SpeciesSelectData %>% filter(Study == name2) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      control <- trainControl(method="repeatedcv",number=3,repeats=5)
      
      test.data$Label <- factor(test.data$Label,levels = c(0,1))
      
      rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(rf.fit, test.data)
      cm<-confusionMatrix(rf.pred,test.data$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name2,BaseModel = name) %>%rbind(confusionMatrixdata)
      
      predob = predict(rf.fit,test.data,type = "prob")
      pred<-prediction(predob[,2],test.data$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name2,BaseModel = name) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(Predict=name2,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
      
      ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name2,BaseModel = name) %>% rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("Study2Study.Top30.",name,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("Study2Study.Top30.",name,".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("Study2Study.Top30.",name,".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = "Study2Study.Top30.RF.model.Importance.csv",row.names = F)


######## Study Transfer Study Top10,20,30,40,50 ########
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
## perform RF model
## Store Data

for (kcount in seq(10,50,10)) {
  ModelImportance <- data.frame()
  confusionMatrixdata <- data.frame()
  modelROC <- data.frame()
  model.AUC <- data.frame()
  
  for (name in Study) {
    ## Process Data
    #middata <- get(paste(name,".wilcoxonsign",sep = ''))
    middata <- read.csv(paste(name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''))
    middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-kcount,adj.fdr.p)
    #write.csv(middataTop50,file = paste(name,"-PairTop50.Species.csv",sep = ''),row.names = F)
    
    SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
    ## train test data
    TrainData <- SpeciesSelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
    ## model for self
    set.seed(123)
    split = sample.split(TrainData$Label,SplitRatio = .7)
    train_data = subset(TrainData,split == TRUE)
    test_data  = subset(TrainData,split == FALSE)
    
    TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
    control <- trainControl(method="repeatedcv",number=2,repeats=5)
    
    train_data$Label <- factor(train_data$Label,levels = c(0,1))
    test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
    fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(fit.rf, TrainData)
    cm<-confusionMatrix(rf.pred,TrainData$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,BaseModel = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(fit.rf,TrainData,type = "prob")
    pred<-prediction(predob[,2],TrainData$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = name,BaseModel = name) %>%rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict=name,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
    
    ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name,BaseModel = name) %>% rbind(ModelImportance)
    
    rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)

    ## model for Other All Studies
    for (name2 in Study) {
      if (name != name2) {
        cat(paste(name,"\t",name2,"\n"))
        
        test.data <- SpeciesSelectData %>% filter(Study == name2) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
          select(-study_condition,-Study)
        
        control <- trainControl(method="repeatedcv",number=3,repeats=5)
        
        test.data$Label <- factor(test.data$Label,levels = c(0,1))
        
        rf.pred <- predict(rf.fit, test.data)
        cm<-confusionMatrix(rf.pred,test.data$Label)
        confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name2,BaseModel = name) %>%rbind(confusionMatrixdata)
        
        predob = predict(rf.fit,test.data,type = "prob")
        pred<-prediction(predob[,2],test.data$Label)
        perf<-performance(pred,'tpr','fpr')
        
        #extrac plot ROC data
        modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name2,BaseModel = name) %>%
          rbind(modelROC)
        
        auc<-performance(pred,"auc")
        auc<-unlist(slot(auc,"y.values"))
        model.AUC <- data.frame(Predict=name2,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
        
        ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
          mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name2,BaseModel = name) %>% rbind(ModelImportance)
      }
    }
  }
  
  write.csv(confusionMatrixdata,file = paste("Study2Study.Top",kcount,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
  write.csv(model.AUC,file = paste("Study2Study.Top",kcount,".RF.model.AUC.csv",sep = ''),row.names = F)
  write.csv(modelROC,file = paste("Study2Study.Top",kcount,".RF.model.ROC.csv",sep = ''),row.names = F)
  write.csv(ModelImportance,file = paste("Study2Study.Top",kcount,".RF.model.Importance.csv",sep = ''),row.names = F)
}

######## Intra/Within Study CrossValidation ########
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()
for (kcount in seq(10,50,10)) {
  for (name in Study) {
    middata <- read.csv(paste(name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''))
    middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-kcount,adj.fdr.p)
    
    SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
    ## Inter Study CV5 
    Traindata <- SpeciesSelectData %>% filter(Study == name) %>% select(c(Study,study_condition,middataTop30$Species)) %>% 
      arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1))
    Traindata$Label <- factor(Traindata$Label,levels = c(0,1))
    
    for (repeatn in 1:10) {
      for (k in 1:5) {
        Folds <- createFolds(y = Traindata$study_condition,k = 5)
        
        TrainData <- Traindata[-Folds[[k]],] %>% data.frame() %>% 
          select(-study_condition,-Study)
        
        TestData <- Traindata[Folds[[k]],] %>% data.frame() %>% 
          select(-study_condition,-Study)
        
        TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
        TestData$Label <- factor(TestData$Label,levels = c(0,1))
        # cross training model
        control <- trainControl(method="repeatedcv",repeats=5)
        
        rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
        
        rf.pred <- predict(rf.fit, TestData)
        cm<-confusionMatrix(rf.pred,TestData$Label)
        confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>%rbind(confusionMatrixdata)
        
        predob = predict(rf.fit,TestData,type = "prob")
        pred<-prediction(predob[,2],TestData$Label)
        perf<-performance(pred,'tpr','fpr')
        
        #extrac plot ROC data
        modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>%
          rbind(modelROC)
        
        auc<-performance(pred,"auc")
        auc<-unlist(slot(auc,"y.values"))
        model.AUC <- data.frame(AUC=auc,Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>% rbind(model.AUC)
        
        ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
          mutate(Rank = floor(rank(-MeanDecreaseGini)),Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>% rbind(ModelImportance)
      }
    }
  }
}

write.csv(confusionMatrixdata,file = paste("IntraStudy.CV5Top10-50",".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("IntraStudy.CV5Top10-50",".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("IntraStudy.CV5Top10-50",".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("IntraStudy.CV5Top10-50",".RF.model.Importance.csv",sep = ''),row.names = F)




######## LODO Species #######
#### Differential Species Analysis ####
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

for (name in Study) {
  middata2 <- SpeciesData2 %>% filter(Study != name) %>% arrange(study_condition)
  groupnum <- table(middata2$study_condition) %>% data.frame()
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  
  rownames(middata) <- c(paste(groupnum$Var1[1]%>% as.vector(),1:groupnum$Freq[1]),paste(groupnum$Var1[2]%>% as.vector(),1:groupnum$Freq[2])) %>% str_remove_all(" ")
  write.csv(data.frame(NewName=rownames(middata),subjectid=rownames(middata2)),file = paste("LODO-Exclude.",name,"-metaphlan2.SubjectID-Newname.csv",sep = ''),row.names = F)
  
  
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential species
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  write.csv(wilcox_res,file = paste("LODO-Exclude.",name,"-Metaphlan2-wilcoxonTest.csv",sep = ''),row.names = F)
  assign(paste(name,".LODO.wilcoxon",sep = ''),wilcox_res)
  
  # pair analysi wilcoxon pair test
  filename=paste("LODO-Exclude.",name,".metaphaln2.pair.csv",sep = '')
  res<-pair_find(data=middata,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
  assign(paste(name,".LODO.wilcoxonsign",sep = ''),res)
  write.csv(res,file = paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),row.names = F)
}

#### LODO RF model Top10-50 ####
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (kcount in seq(10,50,10)) {
  ModelImportance <- data.frame()
  confusionMatrixdata <- data.frame()
  modelROC <- data.frame()
  model.AUC <- data.frame()
  
  for (name in Study) {
    ## Process Data
    #middata <- get(paste(name,".LODO.wilcoxonsign",sep = ''))
    middata <- read.csv(paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''))
    middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-kcount,adj.fdr.p)
    #write.csv(middataTop50,file = paste(name,"-PairTop50.Species.csv",sep = ''),row.names = F)
    
    SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
    ## train test data
    TrainData <- SpeciesSelectData %>% filter(Study != name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
    ## model for self
    set.seed(123)
    split = sample.split(TrainData$Label,SplitRatio = .7)
    train_data = subset(TrainData,split == TRUE)
    test_data  = subset(TrainData,split == FALSE)
    
    TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
    control <- trainControl(method="repeatedcv",number=3,repeats=5)
    
    train_data$Label <- factor(train_data$Label,levels = c(0,1))
    test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
    fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(fit.rf, test_data)
    cm<-confusionMatrix(rf.pred,test_data$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(fit.rf,test_data,type = "prob")
    pred<-prediction(predob[,2],test_data$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict="Self",AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = "Self",ModelExcludeStudy = name) %>% rbind(ModelImportance)
    
    ## model for exclulded study
    TestData <- SpeciesSelectData %>% filter(Study != name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TestData$Label <- factor(TestData$Label,levels = c(0,1))
    
    control <- trainControl(method="repeatedcv",number=5,repeats=5)
    rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(rf.fit, TestData)
    cm<-confusionMatrix(rf.pred,TestData$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(rf.fit,TestData,type = "prob")
    pred<-prediction(predob[,2],TestData$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name,ModelExcludeStudy = name) %>%
      rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict=name,AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = name,ModelExcludeStudy = name) %>% rbind(ModelImportance)
  }
  
  write.csv(confusionMatrixdata,file = paste("LODO.Top",kcount,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
  write.csv(model.AUC,file = paste("LODO.Top",kcount,".RF.model.AUC.csv",sep = ''),row.names = F)
  write.csv(modelROC,file = paste("LODO.Top",kcount,".RF.model.ROC.csv",sep = ''),row.names = F)
  write.csv(ModelImportance,file = paste("LODO.Top",kcount,".RF.model.Importance.csv",sep = ''),row.names = F)
}

######## CrossValidatation RF Model for All Study Pair Top10-50######## 
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

#Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()
for (kcount in seq(10,50,10)){
  middata <- read.csv(paste("Species-8Study-20201010/All-8Study-Contained-Species-pair-wilcoxonsign-res.csv",sep = ''))
  middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-kcount,adj.fdr.p)
  for (i in 1:20) {
    print(i)
    for (k in 1:10) {
      Folds <- createFolds(y = SpeciesData2$study_condition,k = 10)
      #sample index for train and test data
      #Index <- sample(1:10,5,replace = F)
      #two for two
      TrainData <- SpeciesData2[-Folds[[k]],] %>% data.frame() %>% select(c(Study,study_condition,middataTop30$Species)) %>% 
        mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      TestData <- SpeciesData2[Folds[[k]],] %>% data.frame() %>% select(c(Study,study_condition,middataTop30$Species)) %>% 
        mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
      TestData$Label <- factor(TestData$Label,levels = c(0,1))
      # cross training model
      control <- trainControl(method="repeatedcv",number=3,repeats=5)
      
      rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(rf.fit, TestData)
      cm<-confusionMatrix(rf.pred,TestData$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(RepeatTimes = i,Fold=k,FatureCount=kcount) %>%rbind(confusionMatrixdata)
      
      predob = predict(rf.fit,TestData,type = "prob")
      pred<-prediction(predob[,2],TestData$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(RepeatTimes = i,Fold=k,FatureCount=kcount) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(AUC=auc,RepeatTimes = i,Fold=k,FatureCount=kcount) %>% rbind(model.AUC)
      
      ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini)),RepeatTimes = i,Fold=k,FatureCount=kcount) %>% rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("AllStudy.CV10Top10-50",".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("AllStudy.CV10Top10-50",".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("AllStudy.CV10Top10-50",".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("AllStudy.CV10Top10-50",".RF.model.Importance.csv",sep = ''),row.names = F)

######## Sampling Repeat 10 Times #########
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
for (k in 1:10) {
  for (i in SampleNumber) {
    SamplingIndex <- strata(SpeciesData2,stratanames = 'study_condition',size=c(i,i),method = 'srswor')
    SpeciesSelectData <- getdata(SpeciesData2,SamplingIndex) %>% arrange(study_condition)
    write.csv(SpeciesSelectData,file = paste("Sampling-",i,".",k,"-Metaphlan2.Species-Abundance.csv",sep = ""))
    ##
    groupnum <- table(SpeciesSelectData$study_condition) %>% data.frame()
    SpeciesSelectData <- SpeciesSelectData %>% arrange(study_condition) %>% select(-study_condition,-Study,-ID_unit, -Prob,-Stratum)
    rownames(SpeciesSelectData)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
    
    wilcox_res<-matrix(ncol = 7,nrow = 0)
    ## wilcox test differential species
    for (j in 1:dim(SpeciesSelectData)[2]) {
      test<-wilcox.test(SpeciesSelectData[,j][1:groupnum$Freq[1]],SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]])
      wilcox_res <- rbind(wilcox_res,c(colnames(SpeciesSelectData)[j],test$p.value,mean(SpeciesSelectData[,j][1:groupnum$Freq[1]]),median(SpeciesSelectData[,j][1:groupnum$Freq[1]]),mean(SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]]),median(SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]]),j))
    }
    wilcox_res <- as.data.frame(wilcox_res)
    colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
    wilcox_res <- na.omit(wilcox_res)
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
      mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
    wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
    
    wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
    wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
    write.csv2(wilcox_res,file = paste("Sampling-",i,".",k,"-Metaphlan2.Species-wilcoxonTest.csv",sep = ''),row.names = F)
    
    filename=paste("Sampling-",i,".",k,".Speices.pair.csv",sep = '')
    res2<-pair_find(data=SpeciesSelectData,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
    #ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
    write.csv(res2,file = paste("Sampling-",i,".",k,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''),row.names = F)
  }
  
}

######## Sampling RF to Predict Extra All with Repeat 10 Times #########
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
SpeciesData2 <- read.table("Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
# Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()

for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      SampleData <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-Abundance.csv",sep = ''),row.names = 1)
      
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      
      TrainData <- SpeciesData2[rownames(SpeciesData2) %in% rownames(SampleData),] %>% data.frame() %>% 
        mutate(Label = if_else(study_condition == "control",0,1)) %>% select(c(mid2$Species,Label))
      TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
      
      
      TestData <- SpeciesData2[!rownames(SpeciesData2) %in% rownames(SampleData),] %>% data.frame()%>% 
        mutate(Label = if_else(study_condition == "control",0,1)) %>% select(c(mid2$Species,Label))
      TestData$Label <- factor(TestData$Label,levels = c(0,1))
      
      control <- trainControl(method="repeatedcv",repeats=5)
      
      fit.rf <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(fit.rf, TestData)
      cm<-confusionMatrix(rf.pred,TestData$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% 
        mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>% 
        rbind(confusionMatrixdata)
      
      predob = predict(fit.rf,TestData,type = "prob")
      pred<-prediction(predob[,2],TestData$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% 
        mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(AUC=auc) %>% mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>% rbind(model.AUC)
      
      ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini))) %>% mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>%
        rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("Sampling.RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("Sampling.RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("Sampling.RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("Sampling.RF.model.Importance.csv",sep = ''),row.names = F)


######## Draw Picture for Species Analysis #######
#### Diff Species in All Pvalue Analysis like Meta ####
AllDiffData <- read.csv("Species-8Study-20201010/All-8Study-Contained-Species-pair-wilcoxonsign-res.csv",stringsAsFactors = F)

for (feature in seq(10,50,5)) {
  AllDiffSpecies <- AllDiffData %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% select(Species,Enrichment,adj.fdr.p) %>% mutate(Study="All")
  OrderData <- AllDiffData %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% arrange(Enrichment)
  Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
  for (name in Study) {
    mid <- read.csv(paste(name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    AllDiffSpecies<-mid %>% arrange(adj.fdr.p) %>% select(c("Species","Enrichment","adj.fdr.p")) %>% mutate(Study=name) %>% 
      filter(Species %in% OrderData$Species) %>% rbind(AllDiffSpecies)
  }
  AllDiffSpecies <- AllDiffSpecies %>% mutate(if_else(adj.fdr.p == 0,1e-350,adj.fdr.p)) %>% 
    mutate(adj.fdr.p2 = -log10(adj.fdr.p)) %>%
    mutate(Direction = if_else(Enrichment == "Ctl",-adj.fdr.p2,adj.fdr.p2))
  
  AllDiffSpecies$Species <- factor(AllDiffSpecies$Species,levels = rev(OrderData$Species))
  
  p<-ggplot(AllDiffSpecies,aes(x=Species, y=Direction, shape = Study)) + 
    geom_point(aes(color = Enrichment),size = 2.3) + theme_bw() +
    scale_color_manual(name="Enrich", labels = c("Ctl", "CRC","NoDiff"), values = c("Ctl"="blue3", "Disease"="red2","nodiff" = "grey60")) + 
    labs(x="",y="-log10(FDR)") + #scale_y_continuous(limits = c(-200,320))+
    scale_shape_manual(values = c(17,22,3,9,16,8,1,5,20))+
    geom_hline(yintercept=c(-log10(0.01),log10(0.01)),color="purple",linetype = 2)+
    #geom_vline(xintercept=3.5,color="purple",linetype = 2)+
    theme(axis.text.y = element_text(face = "bold.italic"))+
    coord_flip()
  ggsave(p,filename = paste("DiffSpecies-Study.Top",feature,".pdf"),height = 8,width = 7)
  
}
#### Diff TopFeature Pvalue Distribution ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
AllDiffData <- read.csv("Species-8Study-20201010/All-8Study-Contained-Species-pair-wilcoxonsign-res.csv",stringsAsFactors = F)
AllPair <- read.csv("Species-8Study-20201010/All-8Study-Contained.Speices.pair.csv",stringsAsFactors = F)
for (feature in seq(10,50,5)) {
  AllDiffSpecies <- AllDiffData %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% select(Species,Enrichment,adj.fdr.p) %>% 
    mutate(Study="All",PairN=dim(AllPair)[1])
  for (name in Study) {
    mid <- read.csv(paste(name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    PairNumber <- read.csv(paste(name,".metaphaln2.pair.csv",sep = ''),stringsAsFactors = F)
    AllDiffSpecies<-mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>%
      select(c("Species","Enrichment","adj.fdr.p")) %>% mutate(Study=name,PairN=dim(PairNumber)[1]) %>% 
      rbind(AllDiffSpecies)
  }

  p<-ggplot(AllDiffSpecies,aes(x=adj.fdr.p,fill=Study))+geom_histogram(bins = 40)+
    geom_vline(xintercept = c(0.01))+
    facet_wrap(~Study)+theme_bw()+theme(legend.position = "none")+labs(x="FDR",y=paste("Top ",feature," Species Count",sep = ''))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
    theme(axis.text=element_text(size=8),axis.title = element_text(size=12))
  
  ggsave(p,filename = paste("DiffSpecies-Study.Top",feature,".PvalueDistribution.pdf"),height = 4,width = 4)
}

#### Diff Samples Pair Number ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
AllDiffData <- read.csv("Species-8Study-20201010/EightStudies-MetaData.csv",stringsAsFactors = F)
AllPair <- read.csv("Species-8Study-20201010/All-8Study-Contained.Speices.pair.csv",stringsAsFactors = F)
SamplePair <- data.frame(SampleN = dim(AllDiffData)[1],PairN=dim(AllPair)[1])
for (name in Study) {
  SampleData <- read.csv(paste(name,"-metaphlan2.SubjectID-Newname.csv",sep = ''),stringsAsFactors = F)
  SampleP <- read.csv(paste(name,".metaphaln2.pair.csv",sep = ''),stringsAsFactors = F)
  SamplePair <- data.frame(SampleN = dim(SampleData)[1],PairN=dim(SampleP)[1]) %>% rbind(SamplePair)
}

formula <- y ~ x 
p<-ggplot(SamplePair,aes(x=SampleN,y=PairN)) + 
  geom_point() + theme_bw() + #geom_smooth() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #stat_fit_deviations(method = "lm", formula = formula, colour = "red")+ ##library(gginnards)
  geom_smooth(method = "lm", formula = formula) + 
  stat_poly_eq(aes(label =  paste(stat(eq.label), stat(adj.rr.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE,size = 5)+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))

ggsave(p,filename = paste("DiffSpecies-Study.SamplePairNumber.pdf"),height = 4,width = 4)

#### Diff Venn Plot with Meta-Analysis Results ####
#20 Species
Meta <- read.csv("Diff-Species-MetaAnalysis.csv")
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
#note: You need manual assignment value to feature, otherwise it will generate the brokedown pdf. I donnot why.
for (feature in seq(10,50,10)) { #
  # feature=50
  Input_list=list()
  Input_list[["Meta"]] = Meta$Species
  for (name in Study) {
    mid <- read.csv(paste(name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    topSpecies <- mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p)
    Input_list[[name]]<-topSpecies$Species
  }
  
  pdf(file=paste("Diff.Top.",feature,".InstersectMeta.pdf",sep = ''),height = 4,width = 6)
  upset(fromList(Input_list), sets = c(Study,"Meta"), mb.ratio = c(0.5, 0.5), order.by = "freq", 
        nsets = 9, number.angles = 0, point.size = 2, line.size = 1, mainbar.y.label = "Species Count",
        sets.x.label = "Species Count", text.scale = c(2, 2, 2, 2,2, 1.2))
  dev.off()
}

#### Diff RF All CrossValidation 10 Folds ####
# AUC 
AUCdata <- read.csv("AllStudy.CV10Top10-50.RF.model.AUC.csv")
AUCdata <- AUCdata %>% group_by(FatureCount,RepeatTimes) %>% summarise(AUCmean=mean(AUC))
AUCdata$FatureCount <- factor(AUCdata$FatureCount,levels = seq(10,50,10))
p<-ggboxplot(AUCdata,x="FatureCount",y="AUCmean",color = "FatureCount",add = "dotplot",shape="FatureCount")+
  labs(x="Feature Counts",y="Mean AUC CV 10 Folds")+
  theme(legend.position = "none")
ggsave(p,filename = "AllStudy.CV10.AUC.pdf",height = 2,width = 2)

# ROC
ROCdata <- read.csv("AllStudy.CV10Top10-50.RF.model.ROC.csv")
for (repeatn in 1:20) {
  for (featuren in seq(10,50,10)) {
    ROCdata2 <- ROCdata %>% filter(FatureCount == featuren & RepeatTimes == repeatn)
    ROCdata2$Fold <- factor(ROCdata2$Fold,levels = 1:10)
    p<-ggplot(ROCdata2) + geom_path(aes(x=FPR,y=TPR,color=Fold))+
      labs(x = "False positive rate", y = "Ture positive rate") +
      #theme(plot.title = element_text(face = 'bold',size=15))+
      theme_bw() + 
      theme(legend.position = "top")
    ggsave(p,filename = paste("AllStudy.CV10.ROC.",repeatn,".",featuren,".pdf",sep = ''),height = 2,width = 2)
  }
}

# Importance
for (repeatn in 1:20) {
  for (featuren in seq(10,50,10)){
    Importancedata <- read.csv("AllStudy.CV10Top10-50.RF.model.Importance.csv")
    OrderData <- read.csv("Species-8Study-20201010/All-8Study-Contained-Species-pair-wilcoxonsign-res.csv")
    OrderData <- OrderData %>% arrange(adj.fdr.p) %>% top_n(-featuren,adj.fdr.p) %>% arrange(Enrichment)
    Importancedata2 <- Importancedata %>% filter(FatureCount == featuren) %>% filter(RepeatTimes == repeatn)
    Importancedata2 <- Importancedata2 %>% select(Rank,rowname,Fold) %>% spread(Fold,Rank)
    Importancedata2$rowname <- factor(Importancedata2$rowname,levels = OrderData$Species)
    Importancedata2 <- Importancedata2[order(Importancedata2$rowname),] %>% data.frame() %>% remove_rownames() %>% column_to_rownames("rowname")
    
    Number <- table(OrderData$Enrichment) %>% data.frame()
    
    annotation_row = data.frame(EnrichGroup = factor(rep(c("Ctl","CRC"),c(Number$Freq[1],Number$Freq[2]))))
    rownames(annotation_row) <- rownames(Importancedata2)
    anno_color=list(EnrichGroup = c(Ctl="#0066CC",CRC="#CC0000"))
    
    pdf(paste("AllStudy.CV10.ImportanceRank.",repeatn,".",featuren,".pdf",sep = ''),height = 10,width = 10)
    pheatmap(Importancedata2, 
             display_numbers = round(Importancedata2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
             annotation_row = annotation_row, 
             annotation_colors = anno_color,
             color = colorRampPalette(c("#FC4E07", "white","#00AFBB" ))(50),cluster_rows = F,cluster_cols = F,
             border_color = "black",
             cellwidth = 20,
             cellheight = 12,
             #fontsize_row = 5,
             fontsize_col = 12)
    dev.off()
    
  }
}

#### LODO Diff Species Intersect Meta Study  #### 
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
Meta <- read.csv("Diff-Species-MetaAnalysis.csv")

# note : maunal give feature value : 
feature=50
Input_list=list()
Input_list[["Meta"]] = Meta$Species
for (name in Study) {
  mid <- read.csv(paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
  topSpecies <- mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p)
  Input_list[[name]]<-topSpecies$Species
}
pdf(file=paste("LODO.Diff.Top.",feature,".InstersectMeta.pdf",sep = ''),height = 4,width = 6)
upset(fromList(Input_list), sets = c(Study,"Meta"), mb.ratio = c(0.5, 0.5), order.by = "freq", 
      nsets = 9, number.angles = 0, point.size = 2, line.size = 1, mainbar.y.label = "Species Count",
      sets.x.label = "Species Count", text.scale = c(2, 2, 2, 2,2, 1.2))
dev.off()

#### LODO Top Feature P value Distribution ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  AllTopSpecies <- data.frame()
  for (name in Study) {
    mid <- read.csv(paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    topSpecies <- mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% select(Species,adj.fdr.p,Enrichment) %>%
      mutate(Study = name)
    
    AllTopSpecies <- rbind(AllTopSpecies,topSpecies)
  }
  # plot picture
  p<-ggplot(AllTopSpecies,aes(x=adj.fdr.p,fill=Study))+geom_histogram(bins = 40)+
    geom_vline(xintercept = c(0.01))+
    facet_wrap(~Study)+theme_bw()+theme(legend.position = "none")+labs(x="FDR",y=paste("Top ",feature," Species Count",sep = ''))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
    theme(axis.text=element_text(size=6),axis.title = element_text(size=12))
  
  ggsave(p,filename = paste("LODO-DiffSpecies-Study.Top",feature,".PvalueDistribution.pdf"),height = 4,width = 4)
}

Pvalue0Species <- AllTopSpecies %>% filter(adj.fdr.p == 0)
Pvalue0Data <- AllTopSpecies %>% subset(Species %in% Pvalue0Species$Species)
# Fusobacterium_nucleatum, Gemella_morbillorum,Parvimonas_unclassified,Parvimonas_micra,Peptostreptococcus_stomatis,Solobacterium_moorei,Porphyromonas_asaccharolytica
Pvalue0Data2 <- Pvalue0Data %>% select(-Enrichment) %>% spread(Study,adj.fdr.p) %>% remove_rownames() %>% column_to_rownames("Species")
pdf("LODO.Diff.Pvalue0Species.pdf",width = 8,height = 6)
pheatmap(Pvalue0Data2, 
         #display_numbers = matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
         #annotation_row = annotation_row, 
         #annotation_colors = anno_color,
         display_numbers = T,number_format = "%.1e",
         color = colorRampPalette(c("#00AFBB", "white", "#FC4E07"))(50),cluster_rows = F,cluster_cols = F,
         border_color = "black",
         cellwidth = 40,
         cellheight = 20,
         #fontsize_row = 5,
         fontsize_col = 15,
         legend = F,gaps_col = c(1,2,3,4,5,6))
dev.off()

#### LODO Sample Pair Number ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
SamplePair <- data.frame()
for (name in Study) {
  SampleData <- read.csv(paste("LODO-Exclude.",name,"-metaphlan2.SubjectID-Newname.csv",sep = ''),stringsAsFactors = F)
  SampleP <- read.csv(paste("LODO-Exclude.",name,".metaphaln2.pair.csv",sep = ''),stringsAsFactors = F)
  SamplePair <- data.frame(SampleN = dim(SampleData)[1],PairN=dim(SampleP)[1]) %>% rbind(SamplePair)
}

formula <- y ~ x 
p<-ggplot(SamplePair,aes(x=SampleN,y=PairN)) + 
  geom_point() + theme_bw() + #geom_smooth() + 
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #stat_fit_deviations(method = "lm", formula = formula, colour = "red")+ ##library(gginnards)
  geom_smooth(method = "lm", formula = formula) + 
  stat_poly_eq(aes(label =  paste(stat(eq.label), stat(adj.rr.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE,size = 5)+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))

ggsave(p,filename = paste("LODO.DiffSpecies-Study.SamplePairNumber.pdf"),height = 4,width = 4)


#### LODO RF Plot ####
#AUC
AllAUC <- data.frame()
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.AUC.csv",sep = ''))
  AllAUC <- mid %>% mutate(FeatureCount = feature) %>% rbind(AllAUC)
}
AllAUC$FeatureCount <- factor(AllAUC$FeatureCount,levels = seq(10,50,10))
AllAUC$"Condition" = if_else(AllAUC$Predict == "Self","Self","ExcludeStudy")
p<-ggplot(AllAUC%>%filter(Condition == "Self"),aes(FeatureCount,AUC))+
  geom_point(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),size=3)+
  geom_line(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  xlab("No. of features used")+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))

ggsave(p,filename = "LODO-ExcludeStudy.Self.AUC.pdf",height = 4,width = 6)

p<-ggplot(AllAUC%>%filter(Condition == "ExcludeStudy"),aes(FeatureCount,AUC))+
  geom_point(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),size=3)+
  geom_line(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  xlab("No. of features used")+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))+ylim(0.5,1)

ggsave(p,filename = "LODO-ExcludeStudy.ExcludeStudy.AUC2.pdf",height = 4,width = 4)

# ROC
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.ROC.csv",sep = ''))
  p<-ggplot(mid%>%filter(Predict == "Self")) + geom_path(aes(x=FPR,y=TPR,color=ModelExcludeStudy))+
    labs(x = "False positive rate", y = "Ture positive rate") +
    #theme(plot.title = element_text(face = 'bold',size=15))+
    theme_few() + 
    theme(legend.position = "top",legend.title = element_blank())+
    theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
  ggsave(p,filename = paste("LODO-ExcludeStudy.ROC.Self.",feature,".pdf",sep = ''),height = 3,width = 3)
  
  p<-ggplot(mid%>%filter(Predict != "Self")) + geom_path(aes(x=FPR,y=TPR,color=ModelExcludeStudy))+
    labs(x = "False positive rate", y = "Ture positive rate") +
    #theme(plot.title = element_text(face = 'bold',size=15))+
    theme_few() + 
    theme(legend.position = "top",legend.title = element_blank())+
    theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
  ggsave(p,filename = paste("LODO-ExcludeStudy.ROC.ExcludedStudy.",feature,".pdf",sep = ''),height = 3,width = 3)
  
}

# Importance => row name: basemodel; column name: predict model
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.Importance.csv",sep = ''))
  mid2 <- mid %>% filter(Predict == "Self") %>% select(-MeanDecreaseGini,-Predict) %>% spread(ModelExcludeStudy,Rank) %>%
    remove_rownames() %>% column_to_rownames("rowname")
  
  pdf(paste("LODO-ExcludeStudy.Importance.Self.",feature,".pdf",sep = ''),height = 12,width = 10)
  pheatmap(mid2, 
           display_numbers = round(mid2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
           #annotation_row = annotation_row, 
           #annotation_colors = anno_color,
           color = colorRampPalette(c("#FC4E07", "white","#00AFBB" ))(50),cluster_rows = F,cluster_cols = F,
           border_color = "black",
           cellwidth = 30,
           cellheight = 6.5,
           fontsize_row = 6.5,
           fontsize_col = 12,
           legend = T,angle_col = 45)
  dev.off()
  
  mid2 <- mid %>% filter(Predict != "Self") %>% select(-MeanDecreaseGini,-Predict) %>% spread(ModelExcludeStudy,Rank) %>%
    remove_rownames() %>% column_to_rownames("rowname")
  
  pdf(paste("LODO-ExcludeStudy.Importance.ExcludeStudy.",feature,".pdf",sep = ''),height = 12,width = 10)
  pheatmap(mid2, 
           display_numbers = round(mid2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
           #annotation_row = annotation_row, 
           #annotation_colors = anno_color,
           color = colorRampPalette(c("#FC4E07", "white","#00AFBB" ))(50),cluster_rows = F,cluster_cols = F,
           border_color = "black",
           cellwidth = 30,
           cellheight = 6.5,
           fontsize_row = 6.5,
           fontsize_col = 12,
           legend = T,angle_col = 45)
  dev.off()
}

#### Study to Study RF Plot ####
#AUC
#revise the picture need to change the work directory
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.AUC.csv",sep = ''))
  mid2 <- mid %>% spread(Predict,AUC) %>% remove_rownames() %>% column_to_rownames("BaseModel")
  
  pdf(paste("Study2Study.Top",feature,".RF.model.AUC",".pdf",sep = ''),height = 4,width = 6)
  pheatmap(mid2, 
           display_numbers = round(mid2,2), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
           #annotation_row = annotation_row, 
           #annotation_colors = anno_color,
           color = colorRampPalette(c("#00AFBB", "Yellow2"))(50),
           cluster_rows = F,cluster_cols = F,
           border_color = "black",
           cellwidth = 35,
           cellheight = 20,
           #fontsize_row = 5,
           fontsize_col = 12)
  dev.off()
}

#ROC
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.ROC.csv",sep = ''))
  for (name in Study) {
    mid2 <- mid %>% filter(BaseModel == name) %>% select(-BaseModel)
    p<-ggplot(mid2) + geom_path(aes(x=FPR,y=TPR,color=Predict))+
      labs(x = "False positive rate", y = "Ture positive rate") +
      #theme(plot.title = element_text(face = 'bold',size=15))+
      theme_few() + 
      theme(legend.position = "top",legend.title = element_blank())+
      theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
    ggsave(p,filename = paste("Study2Study.ROC.Top.",feature,".","BaseModel.",name,".pdf",sep = ''),height = 3,width = 3)
  }
}

## Importance
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.Importance.csv",sep = ''))
  mid2 <- mid %>% filter(StudyModel != BaseModel) %>% select(rowname,Rank,StudyModel,BaseModel) %>% filter(StudyModel == "PRJEB27928")
  mid2 <- mid %>% filter(StudyModel != BaseModel) %>% select(rowname,Rank,StudyModel,BaseModel) %>% filter(StudyModel == "PRJDB4176") %>%
    filter(BaseModel == "PRJEB27928") %>% rbind(mid2)
  mid3 <- mid2 %>% select(-StudyModel) %>% spread(BaseModel,Rank) %>% remove_rownames() %>% column_to_rownames("rowname")
  write.csv(mid3,paste("Study2Study.Top.",feature,".Integration.csv"))
}

# Every Study
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.Importance.csv",sep = ''))
  for (name in Study) {
    mid2 <- mid %>% filter(BaseModel == name) %>% select(-MeanDecreaseGini,-BaseModel) %>% spread(StudyModel,Rank) %>%
      remove_rownames() %>% column_to_rownames("rowname")
    
    pdf(paste("Study2Study.Top.",feature,".RF.model.Importance.BaseModel.",name,".pdf",sep = ''),height = 10,width = 8)
    pheatmap(mid2, 
             display_numbers = round(mid2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
             #annotation_row = annotation_row, 
             #annotation_colors = anno_color,
             color = colorRampPalette(c("red3","Yellow2"))(50),
             cluster_rows = F,cluster_cols = F,
             border_color = "black",
             cellwidth = 20,
             cellheight = 12,
             #fontsize_row = 5,
             fontsize_col = 10)
    dev.off()
  }
}



#### Sampling Intersect Meta ####
Meta <- read.csv("Diff-Species-MetaAnalysis.csv")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
DataList=list()
for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      Mid[i] <- length(intersect(mid2$Species,Meta$Species))
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  Data <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="IntersectSpecies",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X"))
  
  p<-ggboxplot(Data,x="Sampling",y="IntersectSpecies",color = "Sampling",add = "jitter")+
    labs(x="Sampling",y="Intersect Species Count")+
    theme_few() + theme(legend.position = "none")
  
  ggsave(p,filename = paste("Sampling.Top.",feature,".IntersectMeta.pdf",sep = ''),width = 3,height = 3)
  
}

#### Sampling Sample Pair ####
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
SamplePair <- data.frame()
for (SampleN in SampleNumber) {
  for (i in 1:10) {
    SampleData <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-Abundance.csv",sep = ''))
    PairData <- read.csv(paste("Sampling-",SampleN,".",i,".Speices.pair.csv",sep = ''))
    
    SamplePair<-data.frame(SampleN = dim(SampleData)[1],PairN = dim(PairData)[1]) %>% mutate(RepeatTimes=i,Sampling=SampleN*2) %>% rbind(SamplePair)
  }
}

formula <- y ~ x 
p<-ggplot(SamplePair,aes(x=SampleN,y=PairN)) + 
  geom_point() + theme_bw() + #geom_smooth() + 
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #stat_fit_deviations(method = "lm", formula = formula, colour = "red")+ ##library(gginnards)
  geom_smooth(method = "lm", formula = formula) + 
  stat_poly_eq(aes(label =  paste(stat(eq.label), stat(adj.rr.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE,size = 5)+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))

ggsave(p,filename = paste("Sampling.SamplePairNumber.pdf"),height = 4,width = 4)

#### Sampling Top 10-50 Feature Max Pvlue Distribution ####
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
AllData<-data.frame()
for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      Mid[i] <- max(mid2$adj.fdr.p)
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  AllData <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="MaxPvalue",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X")) %>% mutate(FeatureCount=feature) %>% rbind(AllData)
}
AllData$Sampling <- factor(AllData$Sampling,levels = SampleNumber*2)
AllData$FeatureCount <- factor(AllData$FeatureCount,levels = seq(10,50,10))
p<-ggplot(AllData,aes(x=Sampling,y=MaxPvalue,color=FeatureCount))+geom_boxplot()+
  geom_hline(yintercept = c(0.01))+
  facet_wrap(~FeatureCount)+theme_few()+theme(legend.position = "none")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme(axis.text=element_text(size=6),axis.title = element_text(size=12),axis.text.x = element_text(angle = 45))

ggsave(p,filename = paste("Sampling.Study.TopFeature.MaxPvalueDistribution.pdf"),height = 4,width = 4)

#### Sampling RF ####
# AUC
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
AUCdata <- read.csv("Sampling.RF.model.AUC.csv")

AUCdata$Sampling <- factor(AUCdata$Sampling,levels = SampleNumber*2)
AUCdata$FeatureCount <- factor(AUCdata$FeatureCount,levels = seq(10,50,10))
p<-ggplot(data=AUCdata,aes(Sampling,AUC,color=Sampling)) + geom_boxplot()+geom_jitter()+
  facet_wrap(~FeatureCount) + theme_few() + 
  theme(legend.position = "none")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15),axis.text.x = element_text(angle = 45))

ggsave(p,filename = "Sampling.RF.model.AUC.pdf",height = 6,width = 8)


################# Pathway Differential Analysis ###############
PathwayData <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance.txt",header = T,row.names = 1,sep = '\t')
PathwayData <- PathwayData %>% t() %>% data.frame() %>% rownames_to_column()

MetaData <- read.csv("Pathway-8Study-20201010/EightStudies-MetaData.csv")
MetaData$rowname <- MetaData$rowname %>% str_replace_all("-",".")
PathwayData2 <- merge(MetaData,PathwayData,by = "rowname",all.x = T) %>% data.frame() %>% column_to_rownames("rowname")
PathwayData2 <- as.matrix(PathwayData2)
PathwayData2[is.na(PathwayData2)] = 0
write.table(PathwayData2,file = "Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",quote = F,sep = "\t")
pair_find<-function(data=data,RAN_num=15,RAP_num=30,k="euclidean",rawoff=0.05){
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
  res2 <- res2 %>% data.frame() %>% rename(Pathway=1,pvalue=2,Ctlmean=3,Dismean=4,Ctlmedian=5,Dismedian=6,Enrichment=7) %>%
    mutate(pvalue = as.numeric(as.character(pvalue))) %>%
    arrange(pvalue) %>% na.omit() %>%mutate(pvalue = as.numeric(as.character(pvalue)))
  res2<-na.omit(res2)
  #res2$qvalue = qvalue(p = res2$pvalue)$qvalues
  res2$adj.fdr.p<-p.adjust(res2$pvalue,method="fdr",n=length(res2$pvalue))
  #res2$fdr.qval=fdrtool(res2$pvalue,statistic="pvalue")$qval
  #res2$adj.bonferroni.p <- p.adjust(res2$pvalue,method="bonferroni",n=length(res2$pvalue))
  return(res2)
}

#### All Paired + Wilcoxon ####
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
#PathwayData2 <- PathwayData2 %>% arrange(study_condition)
NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

PathwayData2$study_condition <- factor(PathwayData2$study_condition,levels = c("control","CRC"))
PathwayData2 <- PathwayData2[order(PathwayData2$study_condition),] %>% data.frame()

groupnum <- table(PathwayData2$study_condition) %>% data.frame()
middata1 <- PathwayData2[,-c(1:2)] %>% data.frame() #%>% select(-study_condition,-Study)
filename=paste("Pathway-8Study-20201010/All-8Study-Contained.Pathway.pair.csv",sep = '')
res<-pair_find(data=middata1,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
#ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
res<- merge(res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
write.csv(res,file = paste("Pathway-8Study-20201010/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv",sep = ''),row.names = F)


rownames(middata1)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
wilcox_res<-matrix(ncol = 7,nrow = 0)
## wilcox test differential Pathway
for (j in 1:dim(middata1)[2]) {
  test<-wilcox.test(middata1[,j][1:groupnum$Freq[1]],middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]])
  wilcox_res <- rbind(wilcox_res,c(colnames(middata1)[j],test$p.value,mean(middata1[,j][1:groupnum$Freq[1]]),median(middata1[,j][1:groupnum$Freq[1]]),mean(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),median(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),j))
}
wilcox_res <- as.data.frame(wilcox_res)
colnames(wilcox_res)=c("Pathway","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
wilcox_res <- na.omit(wilcox_res)
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
  mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()

wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
wilcox_res<- merge(wilcox_res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
write.csv(wilcox_res,file = paste("Pathway-8Study-20201010/All-8Study-Contained-Pathway-wilcoxon-test.csv",sep = ''),row.names = F)

#### Every Study Pair + Wilcoxon ####
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
write.csv(NameData,file = "Pathway-8Study-20201010/Pathway.Rename.csv",row.names = F)

colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

PathwayData2$study_condition <- factor(PathwayData2$study_condition,levels = c("control","CRC"))
PathwayData2 <- PathwayData2[order(PathwayData2$study_condition),] %>% data.frame()

Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (name in Study) {
  middata2 <- PathwayData2 %>% filter(Study == name) %>% arrange(study_condition)
  groupnum <- table(middata2$study_condition) %>% data.frame()
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  
  rownames(middata) <- c(paste(groupnum$Var1[1]%>% as.vector(),1:groupnum$Freq[1]),paste(groupnum$Var1[2]%>% as.vector(),1:groupnum$Freq[2])) %>% str_remove_all(" ")
  
  write.csv(data.frame(NewName=rownames(middata),subjectid=rownames(middata2)),file = paste(name,"-Humann2.SubjectID-Newname.csv",sep = ''),row.names = F)
  
  
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential Pathway
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Pathway","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  # merge final Pathway Data
  wilcox_res<- merge(wilcox_res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
  
  write.csv(wilcox_res,file = paste(name,"-Humann2-wilcoxonTest.csv",sep = ''),row.names = F)
  assign(paste(name,".wilcoxon",sep = ''),wilcox_res)
  
  # pair analysi wilcoxon pair test
  filename=paste(name,".metaphaln2.pair.csv",sep = '')
  res<-pair_find(data=middata,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
  assign(paste(name,".wilcoxonsign",sep = ''),res)
  # merge final Pathway Data
  res<- merge(res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
  write.csv(res,file = paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''),row.names = F)
}


######## Study Transfer Study Top50 => neglect ########
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')
## perform RF model
## Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()

for (name in Study) {
  ## Process Data
  #middata <- get(paste(name,".wilcoxonsign",sep = ''))
  middata <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''))
  middataTop50 <- middata %>% arrange(adj.fdr.p) %>% top_n(-50,adj.fdr.p)
  write.csv(middataTop50,file = paste(name,"-PairTop50.Pathway.csv",sep = ''),row.names = F)
  
  PathwaySelectData <- PathwayData2 %>% select(c(Study,study_condition,middataTop50$Pathway)) 
  ## train test data
  TrainData <- PathwaySelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
    select(-study_condition,-Study)
  TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
  
  ## model for self
  set.seed(123)
  split = sample.split(TrainData$Label,SplitRatio = .7)
  train_data = subset(TrainData,split == TRUE)
  test_data  = subset(TrainData,split == FALSE)
  
  TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
  control <- trainControl(method="repeatedcv",number=3,repeats=5)
  
  train_data$Label <- factor(train_data$Label,levels = c(0,1))
  test_data$Label <- factor(test_data$Label,levels = c(0,1))
  
  fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
  
  rf.pred <- predict(fit.rf, TrainData)
  cm<-confusionMatrix(rf.pred,TrainData$Label)
  confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,BaseModel = name) %>%rbind(confusionMatrixdata)
  
  predob = predict(fit.rf,TrainData,type = "prob")
  pred<-prediction(predob[,2],TrainData$Label)
  perf<-performance(pred,'tpr','fpr')
  
  #extrac plot ROC data
  modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = name,BaseModel = name) %>%rbind(modelROC)
  
  auc<-performance(pred,"auc")
  auc<-unlist(slot(auc,"y.values"))
  model.AUC <- data.frame(Predict=name,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
  
  ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
    mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name,BaseModel = name) %>% rbind(ModelImportance)
  
  ## model for Other All Studies
  for (name2 in Study) {
    if (name != name2) {
      cat(paste(name,"\t",name2,"\n"))
      
      test.data <- PathwaySelectData %>% filter(Study == name2) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      control <- trainControl(method="repeatedcv",number=3,repeats=5)
      
      test.data$Label <- factor(test.data$Label,levels = c(0,1))
      
      rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(rf.fit, test.data)
      cm<-confusionMatrix(rf.pred,test.data$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name2,BaseModel = name) %>%rbind(confusionMatrixdata)
      
      predob = predict(rf.fit,test.data,type = "prob")
      pred<-prediction(predob[,2],test.data$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name2,BaseModel = name) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(Predict=name2,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
      
      ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name2,BaseModel = name) %>% rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("Study2Study.Top50.",name,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("Study2Study.Top50.",name,".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("Study2Study.Top50.",name,".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = "Study2Study.Top50.RF.model.Importance.csv",row.names = F)

######## Study Transfer Study Top30 => neglect ########
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
## perform RF model
## Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()

for (name in Study) {
  ## Process Data
  middata <- get(paste(name,".wilcoxonsign",sep = ''))
  middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-30,adj.fdr.p)
  #write.csv(middataTop50,file = paste(name,"-PairTop50.Pathway.csv",sep = ''),row.names = F)
  
  PathwaySelectData <- PathwayData2 %>% select(c(Study,study_condition,middataTop30$Pathway)) 
  ## train test data
  TrainData <- PathwaySelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
    select(-study_condition,-Study)
  TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
  
  ## model for self
  set.seed(123)
  split = sample.split(TrainData$Label,SplitRatio = .7)
  train_data = subset(TrainData,split == TRUE)
  test_data  = subset(TrainData,split == FALSE)
  
  TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
  control <- trainControl(method="repeatedcv",number=3,repeats=5)
  
  train_data$Label <- factor(train_data$Label,levels = c(0,1))
  test_data$Label <- factor(test_data$Label,levels = c(0,1))
  
  fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
  
  rf.pred <- predict(fit.rf, TrainData)
  cm<-confusionMatrix(rf.pred,TrainData$Label)
  confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,BaseModel = name) %>%rbind(confusionMatrixdata)
  
  predob = predict(fit.rf,TrainData,type = "prob")
  pred<-prediction(predob[,2],TrainData$Label)
  perf<-performance(pred,'tpr','fpr')
  
  #extrac plot ROC data
  modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = name,BaseModel = name) %>%rbind(modelROC)
  
  auc<-performance(pred,"auc")
  auc<-unlist(slot(auc,"y.values"))
  model.AUC <- data.frame(Predict=name,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
  
  ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
    mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name,BaseModel = name) %>% rbind(ModelImportance)
  
  ## model for Other All Studies
  for (name2 in Study) {
    if (name != name2) {
      cat(paste(name,"\t",name2,"\n"))
      
      test.data <- PathwaySelectData %>% filter(Study == name2) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      control <- trainControl(method="repeatedcv",number=3,repeats=5)
      
      test.data$Label <- factor(test.data$Label,levels = c(0,1))
      
      rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(rf.fit, test.data)
      cm<-confusionMatrix(rf.pred,test.data$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name2,BaseModel = name) %>%rbind(confusionMatrixdata)
      
      predob = predict(rf.fit,test.data,type = "prob")
      pred<-prediction(predob[,2],test.data$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name2,BaseModel = name) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(Predict=name2,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
      
      ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name2,BaseModel = name) %>% rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("Study2Study.Top30.",name,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("Study2Study.Top30.",name,".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("Study2Study.Top30.",name,".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = "Study2Study.Top30.RF.model.Importance.csv",row.names = F)


######## Study Transfer Study RF model Top10,20,30,40,50 ########
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

## perform RF model
## Store Data
setwd("D:/CRC-Pair/Pathway-8Study-20201010/Res")
for (kcount in seq(10,50,10)) {
  ModelImportance <- data.frame()
  confusionMatrixdata <- data.frame()
  modelROC <- data.frame()
  model.AUC <- data.frame()
  
  for (name in Study) {
    ## Process Data
    #middata <- get(paste(name,".wilcoxonsign",sep = ''))
    middata <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''))
    middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-kcount,adj.fdr.p)
    #write.csv(middataTop50,file = paste(name,"-PairTop50.Pathway.csv",sep = ''),row.names = F)
    
    PathwaySelectData <- PathwayData2 %>% select(c(Study,study_condition,middataTop30$Pathway)) 
    ## train test data
    TrainData <- PathwaySelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
    ## model for self
    split = sample.split(TrainData$Label,SplitRatio = .7)
    train_data = subset(TrainData,split == TRUE)
    test_data  = subset(TrainData,split == FALSE)
    
    TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
    control <- trainControl(method="repeatedcv",number=2,repeats=5)
    
    train_data$Label <- factor(train_data$Label,levels = c(0,1))
    test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
    fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(fit.rf, TrainData)
    cm<-confusionMatrix(rf.pred,TrainData$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,BaseModel = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(fit.rf,TrainData,type = "prob")
    pred<-prediction(predob[,2],TrainData$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = name,BaseModel = name) %>%rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict=name,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
    
    ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name,BaseModel = name) %>% rbind(ModelImportance)
    
    ## model for Other All Studies
    for (name2 in Study) {
      if (name != name2) {
        cat(paste(name,"\t",name2,"\n"))
        
        test.data <- PathwaySelectData %>% filter(Study == name2) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
          select(-study_condition,-Study)
        
        control <- trainControl(method="repeatedcv",repeats=5)
        rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
        
        test.data$Label <- factor(test.data$Label,levels = c(0,1))
        
        rf.pred <- predict(rf.fit, test.data)
        cm<-confusionMatrix(rf.pred,test.data$Label)
        confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name2,BaseModel = name) %>%rbind(confusionMatrixdata)
        
        predob = predict(rf.fit,test.data,type = "prob")
        pred<-prediction(predob[,2],test.data$Label)
        perf<-performance(pred,'tpr','fpr')
        
        #extrac plot ROC data
        modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name2,BaseModel = name) %>%
          rbind(modelROC)
        
        auc<-performance(pred,"auc")
        auc<-unlist(slot(auc,"y.values"))
        model.AUC <- data.frame(Predict=name2,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
        
        ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
          mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name2,BaseModel = name) %>% rbind(ModelImportance)
      }
    }
  }
  
  #colnames(ModelImportance)[1] = "Pathway"
  ModelImportance <- merge(ModelImportance,NameData,by.x = "rowname",by.y = "Rename",all.x =T)

  write.csv(confusionMatrixdata,file = paste("Study22Study.Top",kcount,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
  write.csv(model.AUC,file = paste("Study22Study.Top",kcount,".RF.model.AUC.csv",sep = ''),row.names = F)
  write.csv(modelROC,file = paste("Study22Study.Top",kcount,".RF.model.ROC.csv",sep = ''),row.names = F)
  write.csv(ModelImportance,file = paste("Study22Study.Top",kcount,".RF.model.Importance.csv",sep = ''),row.names = F)
}

######## LODO Pathway #######
#### Differential Pathway Analysis ####
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

for (name in Study) {
  middata2 <- PathwayData2 %>% filter(Study != name) %>% arrange(study_condition)
  groupnum <- table(middata2$study_condition) %>% data.frame()
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  
  rownames(middata) <- c(paste(groupnum$Var1[1]%>% as.vector(),1:groupnum$Freq[1]),paste(groupnum$Var1[2]%>% as.vector(),1:groupnum$Freq[2])) %>% str_remove_all(" ")
  write.csv(data.frame(NewName=rownames(middata),subjectid=rownames(middata2)),file = paste("LODO-Exclude.",name,"-Humann2.SubjectID-Newname.csv",sep = ''),row.names = F)
  
  
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential Pathway
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Pathway","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  wilcox_res <- merge(wilcox_res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
  write.csv(wilcox_res,file = paste("LODO-Exclude.",name,"-Humann2-wilcoxonTest.csv",sep = ''),row.names = F)
  assign(paste(name,".LODO.wilcoxon",sep = ''),wilcox_res)
  
  # pair analysi wilcoxon pair test
  filename=paste("LODO-Exclude.",name,".metaphaln2.pair.csv",sep = '')
  res<-pair_find(data=middata,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
  assign(paste(name,".LODO.wilcoxonsign",sep = ''),res)
  res <- merge(res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
  write.csv(res,file = paste("LODO-Exclude.",name,"-Humann2-PairWilcoxonSign.csv",sep = ''),row.names = F)
}

#### LODO RF model Top10-50 ####
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)

NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (kcount in seq(10,50,10)) {
  ModelImportance <- data.frame()
  confusionMatrixdata <- data.frame()
  modelROC <- data.frame()
  model.AUC <- data.frame()
  
  for (name in Study) {
    ## Process Data
    #middata <- get(paste(name,".LODO.wilcoxonsign",sep = ''))
    middata <- read.csv(paste("LODO-Exclude.",name,"-Humann2-PairWilcoxonSign.csv",sep = ''))
    
    middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-kcount,adj.fdr.p)
    #write.csv(middataTop50,file = paste(name,"-PairTop50.Pathway.csv",sep = ''),row.names = F)
    
    PathwaySelectData <- PathwayData2 %>% select(c(Study,study_condition,middataTop30$Pathway)) 
    ## train test data
    TrainData <- PathwaySelectData %>% filter(Study != name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
    ## model for self
    set.seed(123)
    split = sample.split(TrainData$Label,SplitRatio = .7)
    train_data = subset(TrainData,split == TRUE)
    test_data  = subset(TrainData,split == FALSE)
    
    TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
    control <- trainControl(method="repeatedcv",number=3,repeats=5)
    
    train_data$Label <- factor(train_data$Label,levels = c(0,1))
    test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
    fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(fit.rf, test_data)
    cm<-confusionMatrix(rf.pred,test_data$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(fit.rf,test_data,type = "prob")
    pred<-prediction(predob[,2],test_data$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict="Self",AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = "Self",ModelExcludeStudy = name) %>% rbind(ModelImportance)
    
    ## model for exclulded study
    TestData <- PathwaySelectData %>% filter(Study != name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TestData$Label <- factor(TestData$Label,levels = c(0,1))
    
    control <- trainControl(method="repeatedcv",number=5,repeats=5)
    rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(rf.fit, TestData)
    cm<-confusionMatrix(rf.pred,TestData$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(rf.fit,TestData,type = "prob")
    pred<-prediction(predob[,2],TestData$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name,ModelExcludeStudy = name) %>%
      rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict=name,AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = name,ModelExcludeStudy = name) %>% rbind(ModelImportance)
  }
  
  ModelImportance <- merge(ModelImportance,NameData,by.x = "rowname",by.y = "Rename",all.x =T)
  
  write.csv(confusionMatrixdata,file = paste("LODO.Top",kcount,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
  write.csv(model.AUC,file = paste("LODO.Top",kcount,".RF.model.AUC.csv",sep = ''),row.names = F)
  write.csv(modelROC,file = paste("LODO.Top",kcount,".RF.model.ROC.csv",sep = ''),row.names = F)
  write.csv(ModelImportance,file = paste("LODO.Top",kcount,".RF.model.Importance.csv",sep = ''),row.names = F)
}

######## CrossValidatation RF Model for All Study Pair Top10-50######## 
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

PathwayData2 <- PathwayData2 %>% arrange(study_condition)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

#Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()
for (kcount in seq(10,50,10)){
  middata <- read.csv(paste("Pathway-8Study-20201010/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv",sep = ''))
  middataTop30 <- middata %>% arrange(adj.fdr.p) %>% top_n(-kcount,adj.fdr.p)
  for (i in 1:20) {
    print(i)
    for (k in 1:10) {
      Folds <- createFolds(y = PathwayData2$study_condition,k = 10)
      #sample index for train and test data
      #Index <- sample(1:10,5,replace = F)
      #two for two
      TrainData <- PathwayData2[-Folds[[k]],] %>% data.frame() %>% select(c(Study,study_condition,middataTop30$Pathway)) %>% 
        mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      TestData <- PathwayData2[Folds[[k]],] %>% data.frame() %>% select(c(Study,study_condition,middataTop30$Pathway)) %>% 
        mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
      TestData$Label <- factor(TestData$Label,levels = c(0,1))
      # cross training model
      control <- trainControl(method="repeatedcv",number=3,repeats=5)
      
      rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(rf.fit, TestData)
      cm<-confusionMatrix(rf.pred,TestData$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(RepeatTimes = i,Fold=k,FatureCount=kcount) %>%rbind(confusionMatrixdata)
      
      predob = predict(rf.fit,TestData,type = "prob")
      pred<-prediction(predob[,2],TestData$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(RepeatTimes = i,Fold=k,FatureCount=kcount) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(AUC=auc,RepeatTimes = i,Fold=k,FatureCount=kcount) %>% rbind(model.AUC)
      
      ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini)),RepeatTimes = i,Fold=k,FatureCount=kcount) %>% rbind(ModelImportance)
    }
  }
}
ModelImportance <- merge(ModelImportance,NameData,by.x = "rowname",by.y = "Rename",all.x =T)

write.csv(confusionMatrixdata,file = paste("AllStudy.CV10Top10-50",".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("AllStudy.CV10Top10-50",".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("AllStudy.CV10Top10-50",".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("AllStudy.CV10Top10-50",".RF.model.Importance.csv",sep = ''),row.names = F)

######## Sampling Repeat 10 Times #########
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)

NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

PathwayData2 <- PathwayData2 %>% arrange(study_condition)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
for (k in 1:10) {
  for (i in SampleNumber) {
    SamplingIndex <- strata(PathwayData2,stratanames = 'study_condition',size=c(i,i),method = 'srswor')
    PathwaySelectData <- getdata(PathwayData2,SamplingIndex) %>% arrange(study_condition)
    write.csv(PathwaySelectData,file = paste("Sampling-",i,".",k,"-Humann2.Pathway-Abundance.csv",sep = ""))
    ##
    groupnum <- table(PathwaySelectData$study_condition) %>% data.frame()
    PathwaySelectData <- PathwaySelectData %>% arrange(study_condition) %>% select(-study_condition,-Study,-ID_unit, -Prob,-Stratum)
    rownames(PathwaySelectData)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
    
    wilcox_res<-matrix(ncol = 7,nrow = 0)
    ## wilcox test differential Pathway
    for (j in 1:dim(PathwaySelectData)[2]) {
      test<-wilcox.test(PathwaySelectData[,j][1:groupnum$Freq[1]],PathwaySelectData[,j][(1+groupnum$Freq[1]):dim(PathwaySelectData)[1]])
      wilcox_res <- rbind(wilcox_res,c(colnames(PathwaySelectData)[j],test$p.value,mean(PathwaySelectData[,j][1:groupnum$Freq[1]]),median(PathwaySelectData[,j][1:groupnum$Freq[1]]),mean(PathwaySelectData[,j][(1+groupnum$Freq[1]):dim(PathwaySelectData)[1]]),median(PathwaySelectData[,j][(1+groupnum$Freq[1]):dim(PathwaySelectData)[1]]),j))
    }
    wilcox_res <- as.data.frame(wilcox_res)
    colnames(wilcox_res)=c("Pathway","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
    wilcox_res <- na.omit(wilcox_res)
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
      mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
    wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
    
    wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
    wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
    wilcox_res<- merge(wilcox_res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
    write.csv2(wilcox_res,file = paste("Sampling-",i,".",k,"-Humann2.Pathway-wilcoxonTest.csv",sep = ''),row.names = F)
    
    filename=paste("Sampling-",i,".",k,".Pathway.pair.csv",sep = '')
    res2<-pair_find(data=PathwaySelectData,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
    #ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
    res2<- merge(res2,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
    write.csv(res2,file = paste("Sampling-",i,".",k,"-Humann2.Pathway-PairwilcoxonSign-res.csv",sep = ''),row.names = F)
  }
  
}

######## Sampling RF to Predict Extra All with Repeat 10 Times #########
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)

NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

# Store Data
ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()

for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Humann2.Pathway-PairwilcoxonSign-res.csv",sep = ''))
      SampleData <- read.csv(paste("Sampling-",SampleN,".",i,"-Humann2.Pathway-Abundance.csv",sep = ''),row.names = 1)
      
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      
      TrainData <- PathwayData2[rownames(PathwayData2) %in% rownames(SampleData),] %>% data.frame() %>% 
        mutate(Label = if_else(study_condition == "control",0,1)) %>% select(c(mid2$Pathway,Label))
      TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
      
      
      TestData <- PathwayData2[!rownames(PathwayData2) %in% rownames(SampleData),] %>% data.frame()%>% 
        mutate(Label = if_else(study_condition == "control",0,1)) %>% select(c(mid2$Pathway,Label))
      TestData$Label <- factor(TestData$Label,levels = c(0,1))
      
      control <- trainControl(method="repeatedcv",repeats=5)
      
      fit.rf <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(fit.rf, TestData)
      cm<-confusionMatrix(rf.pred,TestData$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% 
        mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>% 
        rbind(confusionMatrixdata)
      
      predob = predict(fit.rf,TestData,type = "prob")
      pred<-prediction(predob[,2],TestData$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% 
        mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(AUC=auc) %>% mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>% rbind(model.AUC)
      
      ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini))) %>% mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>%
        rbind(ModelImportance)
    }
  }
}

ModelImportance <- merge(ModelImportance,NameData,by.x = "rowname",by.y = "Rename",all.x =T)

write.csv(confusionMatrixdata,file = paste("Sampling.RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("Sampling.RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("Sampling.RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("Sampling.RF.model.Importance.csv",sep = ''),row.names = F)


######## Draw Picture for Pathway Analysis #######
#### Diff Pathway in All Pvalue Analysis like Meta ####
AllDiffData <- read.csv("Pathway-8Study-20201010/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv",stringsAsFactors = F)

for (feature in seq(10,50,5)) {
  AllDiffPathway <- AllDiffData %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% select(Origin,Enrichment,adj.fdr.p) %>% mutate(Study="All")
  OrderData <- AllDiffData %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% arrange(Enrichment)
  Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
  for (name in Study) {
    mid <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    AllDiffPathway<-mid %>% arrange(adj.fdr.p) %>% select(c("Origin","Enrichment","adj.fdr.p")) %>% mutate(Study=name) %>% 
      filter(Origin %in% OrderData$Origin) %>% rbind(AllDiffPathway)
  }
  AllDiffPathway <- AllDiffPathway %>% mutate(if_else(adj.fdr.p == 0,1e-350,adj.fdr.p)) %>% 
    mutate(adj.fdr.p2 = -log10(adj.fdr.p)) %>%
    mutate(Direction = if_else(Enrichment == "Ctl",-adj.fdr.p2,adj.fdr.p2))
  
  AllDiffPathway$Origin <- factor(AllDiffPathway$Origin,levels = rev(OrderData$Origin))
  
  p<-ggplot(AllDiffPathway,aes(x=Origin, y=Direction, shape = Study)) + 
    geom_point(aes(color = Enrichment),size = 2.3) + theme_bw() +
    scale_color_manual(name="Enrich", labels = c("Ctl", "CRC","NoDiff"), values = c("Ctl"="blue3", "Disease"="red2","nodiff" = "grey60")) + 
    labs(x="",y="-log10(FDR)") + #scale_y_continuous(limits = c(-200,320))+
    scale_shape_manual(values = c(17,22,3,9,16,8,1,5,20))+
    geom_hline(yintercept=c(-log10(0.01),log10(0.01)),color="purple",linetype = 2)+
    #geom_vline(xintercept=3.5,color="purple",linetype = 2)+
    theme(axis.text.y = element_text(face = "bold.italic"))+
    coord_flip()
  ggsave(p,filename = paste("DiffPathway-Study.Top",feature,".pdf"),height = 8,width = 7)
  
}
#### Diff TopFeature Pvalue Distribution ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
AllDiffData <- read.csv("Pathway-8Study-20201010/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv",stringsAsFactors = F)
AllPair <- read.csv("Pathway-8Study-20201010/All-8Study-Contained.Pathway.pair.csv",stringsAsFactors = F)
for (feature in seq(10,50,5)) {
  AllDiffPathway <- AllDiffData %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% select(Pathway,Enrichment,adj.fdr.p) %>% 
    mutate(Study="All",PairN=dim(AllPair)[1])
  for (name in Study) {
    mid <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    PairNumber <- read.csv(paste(name,".metaphaln2.pair.csv",sep = ''),stringsAsFactors = F)
    AllDiffPathway<-mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>%
      select(c("Pathway","Enrichment","adj.fdr.p")) %>% mutate(Study=name,PairN=dim(PairNumber)[1]) %>% 
      rbind(AllDiffPathway)
  }
  
  p<-ggplot(AllDiffPathway,aes(x=adj.fdr.p,fill=Study))+geom_histogram(bins = 40)+
    geom_vline(xintercept = c(0.05))+
    facet_wrap(~Study)+theme_bw()+theme(legend.position = "none")+labs(x="FDR",y=paste("Top ",feature," Pathway Count",sep = ''))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
    theme(axis.text=element_text(size=8),axis.title = element_text(size=12))
  
  ggsave(p,filename = paste("DiffPathway-Study.Top",feature,".PvalueDistribution.pdf"),height = 4,width = 4)
}

#### Diff Samples Pair Number ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
AllDiffData <- read.csv("Pathway-8Study-20201010/EightStudies-MetaData.csv",stringsAsFactors = F)
AllPair <- read.csv("Pathway-8Study-20201010/All-8Study-Contained.Pathway.pair.csv",stringsAsFactors = F)
SamplePair <- data.frame(SampleN = dim(AllDiffData)[1],PairN=dim(AllPair)[1])
for (name in Study) {
  SampleData <- read.csv(paste(name,"-Humann2.SubjectID-Newname.csv",sep = ''),stringsAsFactors = F)
  SampleP <- read.csv(paste(name,".metaphaln2.pair.csv",sep = ''),stringsAsFactors = F)
  SamplePair <- data.frame(SampleN = dim(SampleData)[1],PairN=dim(SampleP)[1]) %>% rbind(SamplePair)
}

formula <- y ~ x 
p<-ggplot(SamplePair,aes(x=SampleN,y=PairN)) + 
  geom_point() + theme_bw() + #geom_smooth() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #stat_fit_deviations(method = "lm", formula = formula, colour = "red")+ ##library(gginnards)
  geom_smooth(method = "lm", formula = formula) + 
  stat_poly_eq(aes(label =  paste(stat(eq.label), stat(adj.rr.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE,size = 5)+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))

ggsave(p,filename = paste("DiffPathway-Study.SamplePairNumber.pdf"),height = 4,width = 4)

#### Diff Venn Plot with Meta-Analysis Results ####
#20 Pathway
Meta <- read.csv("Diff-Humann2Pathway-MetaAnalysis.csv")
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
#note: You need manual assignment value to feature, otherwise it will generate the brokedown pdf. I donnot why.
for (feature in seq(10,50,10)) { #
  # feature=10
  Input_list=list()
  Input_list[["Meta"]] = Meta$Pathway
  for (name in Study) {
    mid <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    topPathway <- mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p)
    Input_list[[name]]<-topPathway$Origin
  }
  
  pdf(file=paste("Diff.Top.",feature,".InstersectMeta.pdf",sep = ''),height = 4,width = 6)
  upset(fromList(Input_list), sets = c(Study,"Meta"), mb.ratio = c(0.5, 0.5), order.by = "freq", 
        nsets = 9, number.angles = 0, point.size = 2, line.size = 1, mainbar.y.label = "Pathway Count",
        sets.x.label = "Pathway Count", text.scale = c(2, 2, 2, 2,2, 1.2))
  dev.off()
}

#### Diff RF All CrossValidation 10 Folds ####
# AUC 
AUCdata <- read.csv("AllStudy.CV10Top10-50.RF.model.AUC.csv")
AUCdata <- AUCdata %>% group_by(FatureCount,RepeatTimes) %>% summarise(AUCmean=mean(AUC))
AUCdata$FatureCount <- factor(AUCdata$FatureCount,levels = seq(10,50,10))
p<-ggboxplot(AUCdata,x="FatureCount",y="AUCmean",color = "FatureCount",add = "dotplot",shape="FatureCount")+
  labs(x="Feature Counts",y="Mean AUC CV 10 Folds")+
  theme(legend.position = "none")
ggsave(p,filename = "AllStudy.CV10.AUC.pdf",height = 2,width = 2)

# ROC
ROCdata <- read.csv("AllStudy.CV10Top10-50.RF.model.ROC.csv")
for (repeatn in 1:20) {
  for (featuren in seq(10,50,10)) {
    ROCdata2 <- ROCdata %>% filter(FatureCount == featuren & RepeatTimes == repeatn)
    ROCdata2$Fold <- factor(ROCdata2$Fold,levels = 1:10)
    p<-ggplot(ROCdata2) + geom_path(aes(x=FPR,y=TPR,color=Fold))+
      labs(x = "False positive rate", y = "Ture positive rate") +
      #theme(plot.title = element_text(face = 'bold',size=15))+
      theme_bw() + 
      theme(legend.position = "top")
    ggsave(p,filename = paste("AllStudy.CV10.ROC.",repeatn,".",featuren,".pdf",sep = ''),height = 2,width = 2)
  }
}

# Importance
for (repeatn in 1:20) {
  for (featuren in seq(10,50,10)){
    Importancedata <- read.csv("AllStudy.CV10Top10-50.RF.model.Importance.csv")
    OrderData <- read.csv("Pathway-8Study-20201010/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")
    OrderData <- OrderData %>% arrange(adj.fdr.p) %>% top_n(-featuren,adj.fdr.p) %>% arrange(Enrichment)
    Importancedata2 <- Importancedata %>% filter(FatureCount == featuren) %>% filter(RepeatTimes == repeatn)
    Importancedata2 <- Importancedata2 %>% select(Rank,Origin,Fold) %>% spread(Fold,Rank)
    Importancedata2$Origin <- factor(Importancedata2$Origin,levels = OrderData$Origin)
    Importancedata2 <- Importancedata2[order(Importancedata2$Origin),] %>% data.frame() %>% remove_rownames() %>% column_to_rownames("Origin")
    
    Number <- table(OrderData$Enrichment) %>% data.frame()
    
    annotation_row = data.frame(EnrichGroup = factor(rep(c("Ctl","CRC"),c(Number$Freq[1],Number$Freq[2]))))
    rownames(annotation_row) <- rownames(Importancedata2)
    anno_color=list(EnrichGroup = c(Ctl="#0066CC",CRC="#CC0000"))
    
    pdf(paste("AllStudy.CV10.ImportanceRank.",repeatn,".",featuren,".pdf",sep = ''),height = 10,width = 12)
    pheatmap(Importancedata2, 
             display_numbers = round(Importancedata2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
             annotation_row = annotation_row, 
             annotation_colors = anno_color,
             color = colorRampPalette(c("#FC4E07", "white","#00AFBB" ))(50),cluster_rows = F,cluster_cols = F,
             border_color = "black",
             cellwidth = 20,
             cellheight = 12,
             #fontsize_row = 5,
             fontsize_col = 12)
    dev.off()
    
  }
}

#### LODO Diff Pathway Intersect Meta Study  #### 
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
Meta <- read.csv("Diff-Humann2Pathway-MetaAnalysis.csv")

# note : maunal give feature value : 
feature=10
Input_list=list()
Input_list[["Meta"]] = Meta$Pathway
for (name in Study) {
  mid <- read.csv(paste("LODO-Exclude.",name,"-Humann2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
  topPathway <- mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p)
  Input_list[[name]]<-topPathway$Origin
}
pdf(file=paste("LODO.Diff.Top.",feature,".InstersectMeta.pdf",sep = ''),height = 4,width = 6)
upset(fromList(Input_list), sets = c(Study,"Meta"), mb.ratio = c(0.5, 0.5), order.by = "freq", 
      nsets = 9, number.angles = 0, point.size = 2, line.size = 1, mainbar.y.label = "Pathway Count",
      sets.x.label = "Pathway Count", text.scale = c(2, 2, 2, 2,2, 1.2))
dev.off()

#### LODO Top Feature P value Distribution ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  AllTopPathway <- data.frame()
  for (name in Study) {
    mid <- read.csv(paste("LODO-Exclude.",name,"-Humann2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
    topPathway <- mid %>% arrange(adj.fdr.p) %>% top_n(-feature,adj.fdr.p) %>% select(Origin,adj.fdr.p,Enrichment) %>%
      mutate(Study = name)
    
    AllTopPathway <- rbind(AllTopPathway,topPathway)
  }
  # plot picture
  p<-ggplot(AllTopPathway,aes(x=adj.fdr.p,fill=Study))+geom_histogram(bins = 40)+
    geom_vline(xintercept = c(0.05))+
    facet_wrap(~Study)+theme_bw()+theme(legend.position = "none")+labs(x="FDR",y=paste("Top ",feature," Pathway Count",sep = ''))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
    theme(axis.text=element_text(size=6),axis.title = element_text(size=12))
  
  ggsave(p,filename = paste("LODO-DiffPathway-Study.Top",feature,".PvalueDistribution.pdf"),height = 4,width = 4)
}

#donnot exist Pvalue = 0 Pathway
Pvalue0Pathway <- AllTopPathway %>% filter(adj.fdr.p == 0)
Pvalue0Data <- AllTopPathway %>% subset(Origin %in% Pvalue0Pathway$Origin)
# Fusobacterium_nucleatum, Gemella_morbillorum,Parvimonas_unclassified,Parvimonas_micra,Peptostreptococcus_stomatis,Solobacterium_moorei,Porphyromonas_asaccharolytica
Pvalue0Data2 <- Pvalue0Data %>% select(-Enrichment) %>% spread(Study,adj.fdr.p) %>% remove_rownames() %>% column_to_rownames("Pathway")
pdf("LODO.Diff.Pvalue0Pathway.pdf",width = 8,height = 6)
pheatmap(Pvalue0Data2, 
         #display_numbers = matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
         #annotation_row = annotation_row, 
         #annotation_colors = anno_color,
         display_numbers = T,number_format = "%.1e",
         color = colorRampPalette(c("#00AFBB", "white", "#FC4E07"))(50),cluster_rows = F,cluster_cols = F,
         border_color = "black",
         cellwidth = 40,
         cellheight = 20,
         #fontsize_row = 5,
         fontsize_col = 15,
         legend = F,gaps_col = c(1,2,3,4,5,6))
dev.off()

#### LODO Sample Pair Number ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
SamplePair <- data.frame()
for (name in Study) {
  SampleData <- read.csv(paste("LODO-Exclude.",name,"-Humann2.SubjectID-Newname.csv",sep = ''),stringsAsFactors = F)
  SampleP <- read.csv(paste("LODO-Exclude.",name,".metaphaln2.pair.csv",sep = ''),stringsAsFactors = F)
  SamplePair <- data.frame(SampleN = dim(SampleData)[1],PairN=dim(SampleP)[1]) %>% rbind(SamplePair)
}

formula <- y ~ x 
p<-ggplot(SamplePair,aes(x=SampleN,y=PairN)) + 
  geom_point() + theme_bw() + #geom_smooth() + 
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #stat_fit_deviations(method = "lm", formula = formula, colour = "red")+ ##library(gginnards)
  geom_smooth(method = "lm", formula = formula) + 
  stat_poly_eq(aes(label =  paste(stat(eq.label), stat(adj.rr.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE,size = 5)+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))

ggsave(p,filename = paste("LODO.DiffPathway-Study.SamplePairNumber.pdf"),height = 4,width = 4)


#### LODO RF Plot ####
#AUC
AllAUC <- data.frame()
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.AUC.csv",sep = ''))
  AllAUC <- mid %>% mutate(FeatureCount = feature) %>% rbind(AllAUC)
}
AllAUC$FeatureCount <- factor(AllAUC$FeatureCount,levels = seq(10,50,10))
AllAUC$"Condition" = if_else(AllAUC$Predict == "Self","Self","ExcludeStudy")
p<-ggplot(AllAUC%>%filter(Condition == "Self"),aes(FeatureCount,AUC))+
  geom_point(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),size=3)+
  geom_line(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  xlab("No. of features used")+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))

ggsave(p,filename = "LODO-ExcludeStudy.Self.AUC.pdf",height = 4,width = 6)

p<-ggplot(AllAUC%>%filter(Condition == "ExcludeStudy"),aes(FeatureCount,AUC))+
  geom_point(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),size=3)+
  geom_line(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  xlab("No. of features used")+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))

ggsave(p,filename = "LODO-ExcludeStudy.ExcludeStudy.AUC.pdf",height = 4,width = 6)

# ROC
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.ROC.csv",sep = ''))
  p<-ggplot(mid%>%filter(Predict == "Self")) + geom_path(aes(x=FPR,y=TPR,color=ModelExcludeStudy))+
    labs(x = "False positive rate", y = "Ture positive rate") +
    #theme(plot.title = element_text(face = 'bold',size=15))+
    theme_few() + 
    theme(legend.position = "top",legend.title = element_blank())+
    theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
  ggsave(p,filename = paste("LODO-ExcludeStudy.ROC.Self.",feature,".pdf",sep = ''),height = 3,width = 3)
  
  p<-ggplot(mid%>%filter(Predict != "Self")) + geom_path(aes(x=FPR,y=TPR,color=ModelExcludeStudy))+
    labs(x = "False positive rate", y = "Ture positive rate") +
    #theme(plot.title = element_text(face = 'bold',size=15))+
    theme_few() + 
    theme(legend.position = "top",legend.title = element_blank())+
    theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
  ggsave(p,filename = paste("LODO-ExcludeStudy.ROC.ExcludedStudy.",feature,".pdf",sep = ''),height = 3,width = 3)
  
}

# Importance => row name: basemodel; column name: predict model
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.Importance.csv",sep = ''))
  mid2 <- mid %>% filter(Predict == "Self") %>% select(-MeanDecreaseGini,-Predict) %>% spread(ModelExcludeStudy,Rank) %>%
    remove_rownames() %>% column_to_rownames("rowname")
  
  pdf(paste("LODO-ExcludeStudy.Importance.Self.",feature,".pdf",sep = ''),height = 10,width = 6)
  pheatmap(mid2, 
           display_numbers = round(mid2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
           #annotation_row = annotation_row, 
           #annotation_colors = anno_color,
           color = colorRampPalette(c("#FC4E07", "white","#00AFBB" ))(50),cluster_rows = F,cluster_cols = F,
           border_color = "black",
           cellwidth = 30,
           cellheight = 6.5,
           fontsize_row = 6.5,
           fontsize_col = 12,
           legend = T,angle_col = 45)
  dev.off()
  
  mid2 <- mid %>% filter(Predict != "Self") %>% select(-MeanDecreaseGini,-Predict) %>% spread(ModelExcludeStudy,Rank) %>%
    remove_rownames() %>% column_to_rownames("rowname")
  
  pdf(paste("LODO-ExcludeStudy.Importance.ExcludeStudy.",feature,".pdf",sep = ''),height = 10,width = 6)
  pheatmap(mid2, 
           display_numbers = round(mid2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
           #annotation_row = annotation_row, 
           #annotation_colors = anno_color,
           color = colorRampPalette(c("#FC4E07", "white","#00AFBB" ))(50),cluster_rows = F,cluster_cols = F,
           border_color = "black",
           cellwidth = 30,
           cellheight = 6.5,
           fontsize_row = 6.5,
           fontsize_col = 12,
           legend = T,angle_col = 45)
  dev.off()
}

#### Study to Study RF Plot ####
#AUC
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.AUC.csv",sep = ''))
  mid2 <- mid %>% spread(Predict,AUC) %>% remove_rownames() %>% column_to_rownames("BaseModel")
  
  pdf(paste("Study2Study.Top",feature,".RF.model.AUC",".pdf",sep = ''),height = 4,width = 6)
  pheatmap(mid2, 
           display_numbers = round(mid2,2), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
           #annotation_row = annotation_row, 
           #annotation_colors = anno_color,
           color = colorRampPalette(c("#00AFBB", "Yellow2"))(50),
           cluster_rows = F,cluster_cols = F,
           border_color = "black",
           cellwidth = 35,
           cellheight = 20,
           #fontsize_row = 5,
           fontsize_col = 12)
  dev.off()
}

#ROC
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.ROC.csv",sep = ''))
  for (name in Study) {
    mid2 <- mid %>% filter(BaseModel == name) %>% select(-BaseModel)
    p<-ggplot(mid2) + geom_path(aes(x=FPR,y=TPR,color=Predict))+
      labs(x = "False positive rate", y = "Ture positive rate") +
      #theme(plot.title = element_text(face = 'bold',size=15))+
      theme_few() + 
      theme(legend.position = "top",legend.title = element_blank())+
      theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
    ggsave(p,filename = paste("Study2Study.ROC.Top.",feature,".","BaseModel.",name,".pdf",sep = ''),height = 3,width = 3)
  }
}

## Importance
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.Importance.csv",sep = ''))
  mid2 <- mid %>% filter(StudyModel != BaseModel) %>% select(rowname,Rank,StudyModel,BaseModel) %>% filter(StudyModel == "PRJEB27928")
  mid2 <- mid %>% filter(StudyModel != BaseModel) %>% select(rowname,Rank,StudyModel,BaseModel) %>% filter(StudyModel == "PRJDB4176") %>%
    filter(BaseModel == "PRJEB27928") %>% rbind(mid2)
  mid3 <- mid2 %>% select(-StudyModel) %>% spread(BaseModel,Rank) %>% remove_rownames() %>% column_to_rownames("rowname")
  write.csv(mid3,paste("Study2Study.Top.",feature,".Integration.csv"))
}

# Every Study
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("Study2Study.Top",feature,".RF.model.Importance.csv",sep = ''))
  for (name in Study) {
    mid2 <- mid %>% filter(BaseModel == name) %>% select(-MeanDecreaseGini,-BaseModel) %>% spread(StudyModel,Rank) %>%
      remove_rownames() %>% column_to_rownames("rowname")
    
    pdf(paste("Study2Study.Top.",feature,".RF.model.Importance.BaseModel.",name,".pdf",sep = ''),height = 10,width = 8)
    pheatmap(mid2, 
             display_numbers = round(mid2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
             #annotation_row = annotation_row, 
             #annotation_colors = anno_color,
             color = colorRampPalette(c("red3","Yellow2"))(50),
             cluster_rows = F,cluster_cols = F,
             border_color = "black",
             cellwidth = 20,
             cellheight = 12,
             #fontsize_row = 5,
             fontsize_col = 10)
    dev.off()
  }
}



#### Sampling Intersect Meta ####
Meta <- read.csv("Diff-Pathway-MetaAnalysis.csv")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
DataList=list()
for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Humann2.Pathway-PairwilcoxonSign-res.csv",sep = ''))
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      Mid[i] <- length(intersect(mid2$Pathway,Meta$Pathway))
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  Data <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="IntersectPathway",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X"))
  
  p<-ggboxplot(Data,x="Sampling",y="IntersectPathway",color = "Sampling",add = "jitter")+
    labs(x="Sampling",y="Intersect Pathway Count")+
    theme_few() + theme(legend.position = "none")
  
  ggsave(p,filename = paste("Sampling.Top.",feature,".IntersectMeta.pdf",sep = ''),width = 3,height = 3)
  
}

#### Sampling Sample Pair ####
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
SamplePair <- data.frame()
for (SampleN in SampleNumber) {
  for (i in 1:10) {
    SampleData <- read.csv(paste("Sampling-",SampleN,".",i,"-Humann2.Pathway-Abundance.csv",sep = ''))
    PairData <- read.csv(paste("Sampling-",SampleN,".",i,".Pathway.pair.csv",sep = ''))
    
    SamplePair<-data.frame(SampleN = dim(SampleData)[1],PairN = dim(PairData)[1]) %>% mutate(RepeatTimes=i,Sampling=SampleN*2) %>% rbind(SamplePair)
  }
}

formula <- y ~ x 
p<-ggplot(SamplePair,aes(x=SampleN,y=PairN)) + 
  geom_point() + theme_bw() + #geom_smooth() + 
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #stat_fit_deviations(method = "lm", formula = formula, colour = "red")+ ##library(gginnards)
  geom_smooth(method = "lm", formula = formula) + 
  stat_poly_eq(aes(label =  paste(stat(eq.label), stat(adj.rr.label), sep = "*\", \"*")),
               formula = formula, parse = TRUE,size = 5)+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))

ggsave(p,filename = paste("Sampling.SamplePairNumber.pdf"),height = 4,width = 4)

#### Sampling Top 10-50 Feature Max Pvlue Distribution ####
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
AllData<-data.frame()
for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Humann2.Pathway-PairwilcoxonSign-res.csv",sep = ''))
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      Mid[i] <- max(mid2$adj.fdr.p)
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  AllData <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="MaxPvalue",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X")) %>% mutate(FeatureCount=feature) %>% rbind(AllData)
}
AllData$Sampling <- factor(AllData$Sampling,levels = SampleNumber*2)
AllData$FeatureCount <- factor(AllData$FeatureCount,levels = seq(10,50,10))
p<-ggplot(AllData,aes(x=Sampling,y=MaxPvalue,color=FeatureCount))+geom_boxplot()+
  geom_hline(yintercept = c(0.01))+
  facet_wrap(~FeatureCount)+theme_few()+theme(legend.position = "none")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme(axis.text=element_text(size=6),axis.title = element_text(size=12),axis.text.x = element_text(angle = 45))

ggsave(p,filename = paste("Sampling.Study.TopFeature.MaxPvalueDistribution.pdf"),height = 4,width = 4)

#### Sampling RF ####
# AUC
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
AUCdata <- read.csv("Sampling.RF.model.AUC.csv")

AUCdata$Sampling <- factor(AUCdata$Sampling,levels = SampleNumber*2)
AUCdata$FeatureCount <- factor(AUCdata$FeatureCount,levels = seq(10,50,10))
p<-ggplot(data=AUCdata,aes(Sampling,AUC,color=Sampling)) + geom_boxplot()+geom_jitter()+
  facet_wrap(~FeatureCount) + theme_few() + 
  theme(legend.position = "none")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15),axis.text.x = element_text(angle = 45))

ggsave(p,filename = "Sampling.RF.model.AUC.pdf",height = 6,width = 8)


######## @@ Figures @@ #############
###### @ Figure 1 ########
Fig1Data <- data.frame()
IntersectSpecies <- read_excel("MetaIntersectSpecies.xlsx")
IntersectSpecies$Species <- IntersectSpecies$Species %>% str_replace_all(" ","_")

Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

for (name in Study) {
  index <- which(Study == name)
  middata <- read.csv(paste(name,"-Metaphlan2-wilcoxonTest.csv",sep = ''))
  #Wilcon Test
  Fig1Data  <- na.omit(middata) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0)) %>% 
    filter(Species %in% IntersectSpecies$Species) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".1",sep = ''),Test="Wilcoxon") %>% rbind(Fig1Data )
  #Pair
  middata <- read.csv(paste(name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''))
  Fig1Data  <- na.omit(middata) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0)) %>% 
    filter(Species %in% IntersectSpecies$Species) %>% mutate(Study=name) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".2",sep = ''),Test="Pair") %>% rbind(Fig1Data )
}

#Meta Wilcoxon
MetaData <- read.csv("Species-8Study-20201010/All-8Study-Contained-Species-wilcoxon-test.csv")
Fig1Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0),Study="Meta") %>% 
  filter(Species %in% IntersectSpecies$Species) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.1",Test="Wilcoxon") %>% rbind(Fig1Data )

#Meta Pair
MetaData <- read.csv("Species-8Study-20201010/All-8Study-Contained-Species-pair-wilcoxonsign-res.csv")
Fig1Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0),Study="Meta") %>% 
  filter(Species %in% IntersectSpecies$Species) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.2",Test="Pair") %>% rbind(Fig1Data )
write.csv(Fig1Data,"Figure/Figure1.Data.csv")
#### Fig1a -> Pvalue Heatmap
Pdata <- Fig1Data %>% select(Species,Country,adj.fdr.p,Test) %>% arrange(Country) %>% arrange(desc(Test))
Pdata$Species <- factor(Pdata$Species,levels = IntersectSpecies$Species)
Pdata <- Pdata[order(Pdata$Species),] %>% data.frame() %>% select(-Test) %>% spread(Country,adj.fdr.p) %>% remove_rownames() %>%
  column_to_rownames("Species")
write.csv(Pdata,"Figure/Fig1a.csv")
#note: Manual Change Study Order
Pdata <- read.csv("Figure/Fig1a.csv",row.names = 1)
#Pdata <- t(Pdata)
rownames(Pdata) <- rownames(Pdata) %>% str_replace_all("_"," ")
pdf("Figure/Fig1a.pdf",width = 10,height = 5)
pheatmap(Pdata, 
         display_numbers = T, #round(Pdata,2), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
         number_format = "%.0e",
         number_color = "black",
         fontsize_number = 7.2,
         #annotation_row = annotation_row, 
         #annotation_colors = anno_color,
         color = colorRampPalette(c("white","red" ))(50),
         na_col = "yellow",
         cluster_rows = F,
         cluster_cols = F,
         border_color = "black",
         cellwidth = 25,
         cellheight = 15,
         fontsize_row = 12,
         fontsize_col = 12,
         legend = T,
         angle_col = 45,
         gaps_col = c(8,9,17))
dev.off()

#### Fig1b
Pdata2 <- Fig1Data %>% select(Species,Country,Rank,Test) %>% arrange(Country) %>% arrange(desc(Test))
Pdata2$Species <- factor(Pdata2$Species,levels = IntersectSpecies$Species)
Pdata2 <- Pdata2[order(Pdata2$Species),] %>% data.frame() %>% select(-Test) %>% spread(Country,Rank) %>% remove_rownames() %>%
  column_to_rownames("Species")
write.csv(Pdata2,"Figure/Fig1b.csv")
#note: Manual Change Study Order
Pdata2 <- read.csv("Figure/Fig1b.csv",row.names = 1)
#Pdata <- t(Pdata)
pdf("Figure/Fig1b.pdf",width = 10,height = 5)
pheatmap(Pdata2, 
         display_numbers = round(Pdata2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
         #number_format = "%.0e",
         number_color = "black",
         fontsize_number = 7.5,
         #annotation_row = annotation_row, 
         #annotation_colors = anno_color,
         color = c(colorRampPalette(c("#CD2626","#CD9B9B"))(10),
                   colorRampPalette(c("#CD9B9B","white"))(100)),
         na_col = "yellow",
         cluster_rows = F,
         cluster_cols = F,
         border_color = "black",
         cellwidth = 25,
         cellheight = 15,
         fontsize_row = 12,
         fontsize_col = 12,
         legend = T,
         angle_col = 45,
         gaps_col = c(8,9,17))
dev.off()


###### @ Figure 2 ######
#### Figure 2a => Intra Study CV 5 Repeat 10 Times ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")
NewStudy <- data.frame(Study,StudyPos)
## line plot
IntraData <- read.csv("Species-8Study-20201010/Res/IntraStudy.CV5Top10-50.RF.model.AUC.csv")
IntraData$FeatureCount <- factor(IntraData$FeatureCount,levels = seq(10,50,10))
IntraData2 <- IntraData %>% group_by(Study,FeatureCount) %>% summarise(MeanAUC = mean(AUC))
IntraData2 <- merge(IntraData2,NewStudy,by = "Study",all.x = T)

p1 <- ggplot(IntraData2,aes(FeatureCount,MeanAUC))+
  geom_point(aes(group=StudyPos,color=StudyPos),size=3)+
  geom_line(aes(group=StudyPos,color=StudyPos),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  labs(x="No. of features used",y="AUC")+ylim(0.6,1)+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
ggsave(p1,filename = "Figure/Fig2a.pdf",width = 5,height = 4)

## barplot
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")
NewStudy <- data.frame(Study,StudyPos)
IntraData <- read.csv("IntraStudy.CV5Top10-50.RF.model.AUC.csv")
IntraData$FeatureCount <- factor(IntraData$FeatureCount,levels = seq(10,50,10))
IntraData2 <- IntraData %>% group_by(Study,FeatureCount,RepeatTimes) %>% summarise(MeanAUC = mean(AUC))
IntraData2 <- merge(IntraData2,NewStudy,by = "Study",all.x = T)

p2<-ggbarplot(IntraData2, x = "StudyPos", y = "MeanAUC", 
          add = c("mean_se"),
          color = "FeatureCount",fill = "FeatureCount",
          palette = "npg",
          position = position_dodge(0.8))+theme_few()+
  labs(y="AUC",x="")+scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))+
  theme(legend.position = "top")
ggsave(p2,filename = "Figure/Fig2a-2.pdf",width = 6,height = 4)

#### Figure 2b => Study to Study Transfer ####
###
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # 计算长度
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # 重命名  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  # 计算标准偏差
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
###
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")
NewStudy <- data.frame(Study,StudyPos)
Fig2bData <- data.frame()
for (feature in seq(10,50,10)) {
  data <- read.csv(paste("Study2Study.Top",feature,".RF.model.AUC.csv",sep = ""))
  Fig2bData <- data %>% filter(Predict != BaseModel) %>% select(-Predict) %>% rename(AUC=1,Study=2) %>%
    mutate(Condition="Study2study",FeatureCount = feature) %>% rbind(Fig2bData)
  
  data <- read.csv(paste("LODO.Top",feature,".RF.model.AUC.csv",sep = ''))
  Fig2bData <- data %>% filter(Predict != "Self") %>% select(-Predict) %>% rename(AUC=1,Study=2) %>%
    mutate(Condition="LODOValidation", FeatureCount = feature) %>% rbind(Fig2bData)
}

Fig2bData <- merge(Fig2bData,NewStudy,by="Study",all.x = T)
write.csv(Fig2bData,"Fig2b.data.csv",row.names = F)
# setwd("D:/CRC-Pair")
Fig2bData <- read.csv("Figure/Fig2b.data.csv")
for (feature in seq(10,50,10)) {
  mid <- Fig2bData %>% filter(FeatureCount == feature)
  mid$Condition <- factor(mid$Condition,levels = c("Study2study","LODOValidation"))
  mid2 <- summarySE(mid, measurevar="AUC", groupvars=c("StudyPos","Condition"))
  mid2$Condition <- factor(mid2$Condition,levels = c("Study2study","LODOValidation"))
  p3<-ggplot(mid2, aes(x=StudyPos, y=AUC, fill=Condition)) + 
    geom_bar(position=position_dodge(),stat='identity',colour='black') +
    geom_jitter(data=mid,aes(StudyPos,AUC),color="grey60",position=position_dodge(width=0.85))+
    geom_errorbar(aes(x=StudyPos,ymin=AUC-se, ymax=AUC+se),
                  width=.3,size=1, # 设置误差线的宽度 
                  position=position_dodge(.8))+
    theme_few()+
    labs(y="AUC",x="")+scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))+
    theme(legend.position = "top",legend.title = element_blank())+
    scale_fill_manual(values=c( 'white','grey'))
  
  ggsave(p3,filename = paste("Figure/Fig2b.Features.",feature,".2.pdf",sep = ''),width = 6,height = 5)
}


#### Figure 2c => LODO Importance ####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

for (feature in c(10,20,40,50)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.Importance.csv",sep = ''))
  mid2 <- mid %>% filter(Predict != "Self") %>% select(-MeanDecreaseGini,-Predict) %>% spread(ModelExcludeStudy,Rank) %>%
    remove_rownames() %>% column_to_rownames("rowname")
  
  mid2 <- mid2[apply(mid2, 1, function(x){sum(is.na(x))}) <= 1,] %>% data.frame()
  
  DiffExplore <- data.frame()
  #confirm the species difference
  for (name in Study) {
    data<-read.csv(paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''))
    DiffExplore <- data[data$Species %in% rownames(mid2),] %>% data.frame() %>% select(Species,Enrichment) %>% 
      #mutate(Study=name) %>% 
      rbind(DiffExplore)
  }
  
  #explore enriched groups
  annodata <- DiffExplore %>% group_by(Species) %>% count()%>% merge(DiffExplore,by="Species") %>%top_n(1,n) %>% arrange(Enrichment) %>% unique()
  Number <- table(annodata$Enrichment) %>% data.frame()
  if (! "Ctl" %in% Number$Var1) {
    Number <- rbind(data.frame(Var1="Ctl",Freq=0),Number)
  }
  
  mid2 <- mid2 %>% rownames_to_column()
  mid2$rowname <- factor(mid2$rowname,levels = annodata$Species)
  mid2 <- mid2[order(mid2$rowname),] %>% data.frame() %>% remove_rownames() %>% column_to_rownames("rowname")
  
  #annotation_row = data.frame(EnrichGroup = factor(rep(c("Ctl","CRC"),c(Number$Freq[1],Number$Freq[2]))))
  #rownames(annotation_row) <- rownames(mid2)
  #anno_color=list(EnrichGroup = c(Ctl="#0066CC",CRC="#CC0000"))
  
  #colnames rename
  renameColname <- rep(NA,8)
  for (i in 1:8) {
    index <- which(Study == colnames(mid2)[i])
    renameColname[i] = StudyPos[index]
  }
  colnames(mid2) = renameColname
  
  ranksum <- apply(mid2,1,function(x){sum(na.omit(x))}) %>% sort()
  
  #mid2 <- mid2 %>% rownames_to_column() %>% mutate(rowname = factor(rowname,levels = c("Bifidobacterium_catenulatum",names(ranksum))))  %>%arrange(rowname) %>% remove_rownames() %>% column_to_rownames("rowname")
  
  ranksum2 <- data.frame(ranksum) %>% rownames_to_column() %>%
    merge(annodata,.,by.x = "Species",by.y = "rowname") %>% arrange(Enrichment,ranksum)
  
  mid2 <- mid2 %>% rownames_to_column() %>% mutate(rowname = factor(rowname,levels = ranksum2$Species)) %>%
    arrange(rowname) %>% remove_rownames() %>% column_to_rownames("rowname")
  
  annotation_row = data.frame(EnrichGroup = factor(rep(c("Ctl","CRC"),c(Number$Freq[1],Number$Freq[2]))))
  rownames(annotation_row) <- rownames(mid2)
  anno_color=list(EnrichGroup = c(Ctl="#0066CC",CRC="#CC0000"))
  
  pdf(paste("Figure/Fig2C.LODO-ExcludeStudy.Importance.",feature,".Re.pdf",sep = ''),height = 10,width = 8)
  pheatmap::pheatmap(mid2, 
           display_numbers = round(mid2,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
           annotation_row = annotation_row, 
           annotation_colors = anno_color,
           color = colorRampPalette(c("#8B1A1A","Red","RosyBrown", "white"))(100),cluster_rows = F,cluster_cols = F,
           border_color = "black",
           number_color = "#DCDCDC",
           cellwidth = 18,
           cellheight = 12,
           fontsize_row = 12,
           fontsize_col = 12,
           na_col = "yellow",
           legend = T,angle_col = 45,
           gaps_row = c(Number$Freq[1]))
  dev.off()
}

###### @ Figure 3  ######
#### Figure 3a => Sampling Intersect Meta Species ####
MetaSpecies <- read_excel("MetaIntersectSpecies.xlsx")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
DataList=list()
for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      Mid[i] <- length(intersect(mid2$Species,MetaSpecies$Species))
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  Data <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="IntersectSpecies",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X"))
  
  p<-ggboxplot(Data,x="Sampling",y="IntersectSpecies",color = "Sampling",add = "jitter")+
    labs(x="Sampling",y="Intersect Species Count")+
    theme_few() + theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45))
  
  ggsave(p,filename = paste("Figure/Figure3a.Sampling.Top.",feature,".Species.IntersectMeta.pdf",sep = ''),width = 3,height = 3)
  
}

#### Figure 3b => Sampling AUC ####
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
AUCdata <- read.csv("Species-8Study-20201010/Res/Sampling.RF.model.AUC.csv")

AUCdata$Sampling <- factor(AUCdata$Sampling,levels = SampleNumber*2)
AUCdata$FeatureCount <- factor(AUCdata$FeatureCount,levels = seq(10,50,10))

for (feature in seq(10,50,10)) {
  AUCdata2 <- AUCdata %>% filter(FeatureCount == feature)
  p<-ggplot(data=AUCdata2,aes(Sampling,AUC,color=Sampling)) + geom_boxplot()+geom_jitter(width = 0.1,height = 0.05,size=1)+ 
    theme_few() + 
    theme(legend.position = "none")+ylim(0,1)+
    theme(axis.text = element_text(size=13),axis.title = element_text(size=15),axis.text.x = element_text(angle = 45))
  
  ggsave(p,filename = paste("Figure/Figure3b.Sampling.RF.Top.",feature,".model0-1.AUC.pdf",sep = ""),height = 3,width = 3)
}

###### @ Figure 4  ######
#### Figure 4a ####
setwd("D:/CRC-Pair/Pathway-8Study-20201010/Res")
IntersectSpecies <- read.csv("Diff-Humann2Pathway-MetaAnalysis.csv")
Fig4Data <- data.frame()
IntersectSpecies$Pathway <- IntersectSpecies$Pathway %>% str_replace_all(" ","_")

Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

for (name in Study) {
  index <- which(Study == name)
  middata <- read.csv(paste(name,"-Humann2-wilcoxonTest.csv",sep = ''))
  #Wilcon Test
  Fig4Data  <- na.omit(middata) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0)) %>% 
    filter(Origin %in% IntersectSpecies$Pathway) %>% select(Origin,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".1",sep = ''),Test="Wilcoxon") %>% rbind(Fig4Data )
  #Pair
  middata <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''))
  Fig4Data  <- na.omit(middata) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0)) %>% 
    filter(Origin %in% IntersectSpecies$Pathway) %>% mutate(Study=name) %>% select(Origin,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".2",sep = ''),Test="Pair") %>% rbind(Fig4Data )
}


#Meta Wilcoxon
MetaData <- read.csv("../All-8Study-Contained-Pathway-wilcoxon-test.csv")
Fig4Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0),Study="Meta") %>% 
  filter(Origin %in% IntersectSpecies$Pathway) %>% select(Origin,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.1",Test="Wilcoxon") %>% rbind(Fig4Data )

#Meta Pair
MetaData <- read.csv("../All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")
Fig4Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p),0),Study="Meta") %>% 
  filter(Origin %in% IntersectSpecies$Pathway) %>% select(Origin,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.2",Test="Pair") %>% rbind(Fig4Data )
write.csv(Fig4Data,"../../Figure/Figure4.Data.csv",row.names = F)
## Pvalue Heatmap
Pdata <- Fig4Data %>% select(Origin,Country,adj.fdr.p,Test) %>% arrange(Country) %>% arrange(desc(Test))
Pdata$Origin <- factor(Pdata$Origin,levels = IntersectSpecies$Pathway)
Pdata <- Pdata[order(Pdata$Origin),] %>% data.frame() %>% select(-Test) %>% spread(Country,adj.fdr.p) %>% remove_rownames() %>%
  column_to_rownames("Origin")
write.csv(Pdata,"../../Figure/Fig4a.csv")
#note: Manual Change Study Order
Pdata <- read.csv("../../Figure/Fig4a.csv")
#Pdata <- t(Pdata)
Pdata <- merge(IntersectSpecies,Pdata,by.x = "Pathway",by.y = "X",all = T)
Pdata2 <- Pdata %>% data.frame() %>% arrange(Enrich) %>% select(-Pathway,-Enrich) %>% remove_rownames() %>% column_to_rownames("MetaPathway")

annotation_row = data.frame(EnrichGroup = factor(rep(c("CRC","Ctl"),c(8,12))))
rownames(annotation_row) <- rownames(Pdata2)
anno_color=list(EnrichGroup = c(Ctl="#0066CC",CRC="#CC0000"))

pdf("../../Figure/Fig4a.pdf",width = 12,height = 5)
pheatmap(Pdata2, 
         display_numbers = T, #round(Pdata,2), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
         number_format = "%.0e",
         number_color = "black",
         fontsize_number = 6,
         annotation_row = annotation_row, 
         annotation_colors = anno_color,
         color = colorRampPalette(c("white","red" ))(50),
         na_col = "yellow",
         cluster_rows = F,
         cluster_cols = F,
         border_color = "black",
         cellwidth = 20,
         cellheight = 15,
         fontsize_row = 12,
         fontsize_col = 12,
         legend = T,
         angle_col = 45,
         gaps_col = c(8,9,17),gaps_row = c(8))
dev.off()

#### Figure 4b ####
setwd("D:/CRC-Pair/Pathway-8Study-20201010/Res")
PairData <- read.csv("../All-8Study-Contained.Pathway.pair.csv")
AbunData <- read.table("../EightStudies-PathwayAbundance-Group.txt",row.names = 1,header = T,sep = '\t')

#PWY.5667..CDP.diacylglycerol.biosynthesis.I #Pathway223
#PWY0.1319..CDP.diacylglycerol.biosynthesis.II #Pathway467
#PWY.4984..urea.cycle #Pathway162
#PWY.6123..inosine.5..phosphate.biosynthesis.I #Pathway275
#PWY.6124..inosine.5..phosphate.biosynthesis.II #Pathway276
#PWY.7234..inosine.5..phosphate.biosynthesis.III #Pathway403
middata <- AbunData[rownames(AbunData) %in% c(unique(PairData$Ctl),unique(PairData$Disease)),]
middata2 <- middata[,c(1,2,225,164,277,278,405,469)]

p1<-ggplot(middata2,aes(x=study_condition,y=PWY.5667..CDP.diacylglycerol.biosynthesis.I,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "CDP-diacylglycerol biosynthesis I")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p1,filename = "../../Figure/Figure4b.CDP.pdf",width = 3.5,height = 3)

p6<-ggplot(middata2,aes(x=study_condition,y=PWY0.1319..CDP.diacylglycerol.biosynthesis.II,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "CDP-diacylglycerol biosynthesis II")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p6,filename = "../../Figure/Figure4b.CDPII.pdf",width = 3.5,height = 3)

p2<-ggplot(middata2,aes(x=study_condition,y=PWY.4984..urea.cycle,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "urea cycle")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p2,filename = "../../Figure/Figure4b.ureacycle.pdf",width = 3.5,height = 3)

p3<-ggplot(middata2,aes(x=study_condition,y=PWY.6123..inosine.5..phosphate.biosynthesis.I,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis I")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p3,filename = "../../Figure/Figure4b.IMP1.pdf",width = 3.5,height = 3)

p4<-ggplot(middata2,aes(x=study_condition,y=PWY.6124..inosine.5..phosphate.biosynthesis.II,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis II")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p4,filename = "../../Figure/Figure4b.IMP2.pdf",width = 3.5,height = 3)

p5<-ggplot(middata2,aes(x=study_condition,y=PWY.7234..inosine.5..phosphate.biosynthesis.III,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis III")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p5,filename = "../../Figure/Figure4b.IMP3.pdf",width = 3.5,height = 3)

## Sum the same pathway
middata2$"IMP biosynthesis I+II+III" <- middata2$PWY.6123..inosine.5..phosphate.biosynthesis.I+
  middata2$PWY.6124..inosine.5..phosphate.biosynthesis.II+
  middata2$PWY.7234..inosine.5..phosphate.biosynthesis.III
  
p8<-ggplot(middata2,aes(x=study_condition,y=`IMP biosynthesis I+II+III`,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis I+II+III")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))#+stat_compare_means(comparisons=c("control","CRC"),label = "p.value",hide.ns = F)

ggsave(p8,filename = "../../Figure/Figure4b.IMP1+2+3.pdf",width = 3.5,height = 3)

middata2$"CDP-diacylglycerol biosynthesis I+II" = middata2$PWY.5667..CDP.diacylglycerol.biosynthesis.I + 
  middata2$PWY0.1319..CDP.diacylglycerol.biosynthesis.II

p7<-ggplot(middata2,aes(x=study_condition,y=`CDP-diacylglycerol biosynthesis I+II`,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "CDP-diacylglycerol biosynthesis I+II")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p7,filename = "../../Figure/Figure4b.CDP1+2.pdf",width = 3.5,height = 3)


####  ######

######## @@ Figures - Revise @@ #############
#### Fig1a - Revise ####
Pdata.1 <- read.csv("Figure/Fig1a.csv",row.names = 1,stringsAsFactors = F)
Pdata.1 <- Pdata.1 %>% arrange(Meta.2,Meta.1)
#Pdata <- t(Pdata)
rownames(Pdata.1) <- rownames(Pdata.1) %>% str_replace_all("_"," ") 
Pdata <- Pdata %>% rownames_to_column() %>%
  gather(key = "Condition",value = "FDR",-rowname)

labels = c("<10^-5","10^-4","10^-3","10^-2","0.05","1")
breaks <- c(-1,10^-5,10^-4,10^-3,10^-2,0.05,1)
breaks2 <- c(10^-5,10^-4,10^-3,10^-2,0.05,1)
mid <- data.frame(labels=labels,breaks=breaks2)
Pdata$labels = cut(Pdata$FDR,breaks,labels,ordered_result = T)
Pdata<-merge(Pdata,mid,by= "labels")
Pdata$alpha <- Pdata$FDR/Pdata$breaks
Pdata$rowname <-factor(Pdata$rowname,levels = rownames(Pdata.1)%>% str_replace_all("_"," "))
interval.cols <- c("#8B1A1A","#FF3030","#FF6A6A","#FFC1C1","#CDC9C9","#EEE9E9")#brewer.pal(6,"Set2")
names(interval.cols) <- levels(Pdata$labels)
ggplot(Pdata,aes(x=Condition,y=rowname,fill=labels))+
  geom_tile(color=NA)+theme_few()+
  theme(strip.text.x=element_text(angle=90,vjust=1,hjust=0.5,size=6),
        panel.spacing=unit(0.00,"cm"),plot.margin=unit(c(1,1,1,1),"cm"),
        legend.key.size=unit(0.25,"cm"),
        panel.border=element_blank(),
        strip.background=element_blank(),
        axis.ticks.y=element_line(size=0.25))+
  labs(x="",y="")+
  scale_color_manual(drop=FALSE,values=interval.cols,labels=names(interval.cols),name="FDR")+
  scale_fill_manual(drop=FALSE,values=interval.cols,labels=names(interval.cols),name="FDR")


ggplot(Pdata,aes(x=rowname, y=1, fill=FDR)) +
  geom_tile() +
  facet_grid(Condition~.) +
  theme_few() +
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA, colour=NA),
        axis.ticks.x=element_blank(),
        strip.background = element_blank()) +
  scale_fill_gradientn(colours=c(paste0("red",1:3),paste0('grey', seq(20, 90, 5))),na.value = 'red') +
  xlab('') + ylab('') +
  geom_text(aes(label=FDR), colour='white')

##
Pdata.2 <- Pdata %>% select(labels,rowname,Condition) %>% spread(Condition,labels) %>% 
  mutate(rowname = factor(rowname,levels = rownames(Pdata.1)%>% str_replace_all("_"," "))) %>% arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>%
  select(colnames(Pdata.1))
## text note => False
TextFunc <- function(dat, col = "black", fontsize = 8, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}
###
#Pdata <- read.csv("Figure/Fig1a.Revise.csv",row.names = 1)
col_cat <- c("<10^-5"="#8B1A1A","10^-4"="#EE0000","10^-3"="#FF6A6A","10^-2"="#FFC1C1","0.05"="#CDC9C9","1"="#EEE9E9")
pdf("Figure/Fig1a.Re.pdf",width = 10,height = 4)
Heatmap(Pdata.2, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1)
        )
dev.off()

#### Fig1b - Revise ####
dPdata2 <- read.csv("Figure/Fig1b.csv",row.names = 1)

dPdata3 <- dPdata2 %>% rownames_to_column() %>%
  gather(key = "Condition",value = "FDR",-rowname)


labels = c("1-10","10-20","20-30","30-40","40-50",">50")
breaks <- c(0,10,20,30,40,50,500)
dPdata3$labels = cut(dPdata3$FDR,breaks,labels,ordered_result = T)

TextFunc <- function(dat, col = "black", fontsize = 8, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(round(dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(round(dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}

dPdata4 <- dPdata3 %>% select(labels,rowname,Condition) %>% spread(Condition,labels) %>% 
  mutate(rowname = factor(rowname,levels = rownames(Pdata.1))) %>% arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>%
  select(colnames(Pdata.1))

dPdata5 <- dPdata2 %>% rownames_to_column() %>% mutate(rowname = factor(rowname,levels = rownames(Pdata.1))) %>%
  arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>%
  select(colnames(Pdata.1))

#Pdata <- t(Pdata)
#pdf("Figure/Fig1b.pdf",width = 10,height = 5)
col_cat <- c("1-10"="#8B1A1A","10-20"="#EE0000","20-30"="#FF6A6A","30-40"="#FFC1C1","40-50"="#CDC9C9",">50"="#EEE9E9")
pdf("Figure/Fig1b.Re.pdf",width = 8,height = 3)
Heatmap(dPdata4, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(dPdata5)
)
dev.off()





#### Sampling Figure Intersect with 7 Species ####
Species7 <- c("Parvimonas_unclassified","Gemella_morbillorum","Peptostreptococcus_stomatis",
              "Fusobacterium_nucleatum","Parvimonas_micra","Porphyromonas_asaccharolytica",
              "Clostridium_symbiosum")

SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
DataList=list()
for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      mid <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      mid2 <- mid %>% top_n(-feature,adj.fdr.p) %>% arrange(adj.fdr.p) %>% arrange(Enrichment)
      Mid[i] <- length(intersect(mid2$Species,Species7))
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  Data <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="IntersectSpecies",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X"))
  
  p<-ggboxplot(Data,x="Sampling",y="IntersectSpecies",color = "Sampling",add = "jitter")+
    labs(x="Sampling",y="Intersect Species Count")+
    theme_few() + theme(legend.position = "none")
  
  ggsave(p,filename = paste("Sampling.Top.",feature,".IntersectMeta7.pdf",sep = ''),width = 3,height = 3)
  
}

#############

#### Fig4a ValconoPlot ####
data <- read.csv("Pathway-8Study-20201010/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv",check.names = F,
                 stringsAsFactors = F)
data$Label <- str_replace_all(data$Label,"\\?","\\-")
data$log2FC <- log2(data$Dismean/data$Ctlmean)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
color <- sample(color, 13)

p<-ggplot(data,aes(x=log2FC,y=-log10(adj.fdr.p),shape=Color,color=Color))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(1e-22)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "top")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'FDR'))

ggsave(p,filename = "Fig4a-Left2.pdf",height = 5.6,width = 5)  

# revise
data$Size <- factor(data$Size,levels = c(1,3))
p<-ggplot(data,aes(x=log2FC,y=-log10(adj.fdr.p),color=Color,shape=Size))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  scale_x_continuous(limits = c(-2,4))+
  scale_shape_manual(values=c(17,19))+
  geom_hline(aes(yintercept=-log10(1e-22)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "top")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'FDR'))

ggsave(p,filename = "Fig4a-Left2.pdf",height = 5.6,width = 5)  

###
ggplot(data,mapping=aes(x=log2FC,y=-log10(adj.fdr.p),shape=Color,color=Label,size=Size))+geom_point(size=2,alpha=0.8)+theme_few()+
  scale_color_manual(values = c("grey",color))+
  scale_size_continuous(range = c(1,10))+
  #geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),alpha=1)+
  scale_x_continuous(limits = c(-2,2))+
  geom_hline(aes(yintercept=-log10(1e-22)),linetype=2,color="black")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'FDR'))

###$
data <- read.csv("Pathway-8Study-20201010/All-8Study-Contained-Pathway-wilcoxon-test.csv",check.names = F,
                 stringsAsFactors = F)

data$Label <- str_replace_all(data$Label,"\\?","\\-")
data$log2FC <- log2(data$CRCmean/data$Ctlmean)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
color <- sample(color, 13)

p<-ggplot(data,aes(x=log2FC,y=-log10(adj.fdr.p),shape=Color,color=Color))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.001)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "top")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'FDR'))

ggsave(p,filename = "Fig4a-Right.pdf",height = 5.6,width = 5) 

# revise
p<-ggplot(data,aes(x=log2FC,y=-log10(adj.fdr.p),color=Color))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#BC3C28","#0072B5","grey60"))+
  scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.001)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "top")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'FDR'))

ggsave(p,filename = "Fig4a-Right2.pdf",height = 5.6,width = 5)

### Fig4a Rev ####
setwd("D:/CRC-Pair/Pathway-8Study-20201010/Res")
Pdata.1 <- read.csv("../../Figure/Fig4a.csv",row.names = 1,stringsAsFactors = F)
Pdata.1<-Pdata.1[,-1]
#Pdata.1 <- Pdata.1 %>% arrange(Meta.2,Meta.1)
#Pdata <- t(Pdata)
#rownames(Pdata.1) <- rownames(Pdata.1) %>% str_replace_all("_"," ")
Pdata <- read.csv("../../Figure/Figure4.Data.csv",check.names = F,stringsAsFactors = F)
Pdata <- Pdata%>%select(Origin,adj.fdr.p,Country) #%>% gather(key = "Condition",value = "FDR",-Origin)

labels = c("<10^-5","10^-4","10^-3","10^-2","0.05","1")
breaks <- c(-1,10^-5,10^-4,10^-3,10^-2,0.05,1)
breaks2 <- c(10^-5,10^-4,10^-3,10^-2,0.05,1)
mid <- data.frame(labels=labels,breaks=breaks2)
Pdata$adj.fdr.p <- as.numeric(as.character(Pdata$adj.fdr.p))
Pdata$labels = cut(Pdata$adj.fdr.p,breaks,labels,ordered_result = T)
Pdata<-merge(Pdata,mid,by= "labels")
Pdata$Origin <-factor(Pdata$Origin,levels = rownames(Pdata.1)%>% str_replace_all("_"," "))
interval.cols <- c("#8B1A1A","#FF3030","#FF6A6A","#FFC1C1","#CDC9C9","#EEE9E9")#brewer.pal(6,"Set2")
names(interval.cols) <- levels(Pdata$labels)

Pdata.2 <- Pdata %>% select(labels,Origin,Country) %>% spread(Country,labels) %>% 
  mutate(Origin = factor(Origin,levels = rownames(Pdata.1)%>% str_replace_all("_"," "))) %>% arrange(Origin) %>%
  remove_rownames() %>% column_to_rownames("Origin") %>%
  select(colnames(Pdata.1))
## text note => False
TextFunc <- function(dat, col = "black", fontsize = 8, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}
###
#Pdata <- read.csv("Figure/Fig1a.Revise.csv",row.names = 1)
col_cat <- c("<10^-5"="#8B1A1A","10^-4"="#EE0000","10^-3"="#FF6A6A","10^-2"="#FFC1C1","0.05"="#CDC9C9","1"="#EEE9E9")
setwd("D:/CRC-Pair")
pdf("Fig4a.Rev.pdf",width = 10,height = 4)
Heatmap(Pdata.2, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1)
)
dev.off()


####### Pathway









































##### ValconoPlot Meta #####
setwd("D:/CRC-Pair/Species-8Study-20201010")
data <- read.csv("All-8Study-Contained-Species-wilcoxon-test.csv")
data$log2FC <- log2(data$CRCmean/data$Ctlmean)
p<-ggplot(data,aes(x=log2FC,y=-log10(adj.fdr.p),color=Color))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#BC3C28","grey60"))+
  scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.00001)),linetype=2,color="black")+
  geom_vline(aes(xintercept=0),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),
        legend.position = "top")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'FDR'))
ggsave(p,filename = "FigS3-Left.pdf",height = 5.6,width = 5)

####
data <- read.csv("Species-8Study-20201010/All-8Study-Contained-Species-pair-wilcoxonsign-res.csv")
data$log2FC <- log2(data$Dismean/data$Ctlmean)
data$Size <- factor(data$Size,levels = c(1,3))
p<-ggplot(data,aes(x=log2FC,y=-log10(adj.fdr.p),color=Color,shape=Size))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#BC3C28","#0072B5","grey60"))+
  scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  scale_shape_manual(values=c(17,19))+
  geom_hline(aes(yintercept=-log10(1e-75)),linetype=2,color="black")+
  geom_vline(aes(xintercept=0),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),
        legend.position = "top")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'FDR'))

ggsave(p,filename = "FigS3-Final.pdf",height = 5.6,width = 5)

######### Species Pair Transfer Pathways ############
PairData <- read.csv("Species-8Study-20201010/All-8Study-Contained.Speices.pair.csv")

PathwayData2 <- read.table("Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
#write.csv(NameData,file = "Pathway-8Study-20201010/Pathway.Rename.csv",row.names = F)

colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

PathwayData2$study_condition <- factor(PathwayData2$study_condition,levels = c("control","CRC"))
PathwayData2 <- PathwayData2[order(PathwayData2$study_condition),] %>% data.frame()

name=c("PathwayName","Origin","Rank","Pvalue","CtrlMean","CtrlMedian","CRCmean","CRCmedian")
res <- matrix(ncol = 8,nrow = 0)
for (j in 3:dim(PathwayData2)[2]) {
  Ctrl <- rep(NA,dim(PairData)[1])
  CRC <- rep(NA,dim(PairData)[1])
  for (i in 1:dim(PairData)[1]) {
    ctrindex <- which(rownames(PathwayData2) == PairData$Ctl[i])
    Ctrl[i] = as.numeric(as.character(PathwayData2[ctrindex,j]))
    
    crcindex <- which(rownames(PathwayData2) == PairData$Disease[i])
    CRC[i] = as.numeric(as.character(PathwayData2[crcindex,j]))
  }

  ctrdata <- PathwayData2 %>% filter(study_condition == "control")
  crcdata <- PathwayData2 %>% filter(study_condition == "CRC")
  
  test <- wilcox.test(Ctrl,CRC,paired = T,exact = F)
  if (sum(Ctrl > CRC) > sum(Ctrl < CRC)) {
    Rank = "Ctrl"
  }else if(sum(Ctrl > CRC) < sum(Ctrl < CRC)){
    Rank = "CRC"
  }else {
    Rank = "No"
  }
  res <- rbind(res,c(colnames(PathwayData2)[j],NameData$Origin[j-2],Rank,test$p.value,mean(ctrdata[,j]),median(ctrdata[,j]),mean(crcdata[,j]),median(crcdata[,j]),median(crcdata[,j])))
}
colnames(res) = c("PathwayName","Origin","Rank","Pvalue","CtrlMean","CtrlMedian","CRCmean","CRCmedian")
res <- as.data.frame(res)
res$Pvalue <- as.numeric(as.character(res$Pvalue))
res <- res %>% arrange(Pvalue) %>% mutate(FDR = p.adjust(Pvalue)) %>% mutate(Enrich = ifelse(FDR <= 0.001,ifelse(Rank == "Ctrl","Ctrl","CRC"),"N.S."))
write.csv(res,"SpeciesPair-TransferPathway.csv",row.names = F)




















