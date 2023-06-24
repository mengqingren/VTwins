## VTwins
The goal of `VTwins` is to perform phenotype-enriched features including species or pathways for metagenomic sequencing or 16S sequencing.

if you want to use standalone version in Linux, you can check VTwins.Linux from [mengqingren/VTwins.Linux](https://github.com/mengqingren/VTwins.Linux).

You can also find the R script of this paper [mengqingren/VTwins.Linux](https://github.com/mengqingren/VTwins.Linux).

if you want to use R package, you can keep reading. 

### Installation
Now `VTwins` is not on cran, You can install the development version of
VTwins from [GitHub](https://github.com/) with:
  
  ``` r
# install.packages("remotes")
remotes::install_github("mengqingren/VTwins")
```
### Dependence
- R version 4.0.5
- tidyverse_1.3.1

### Basic Usage
``` r
library(VTwins)
pair_find(data=YourRelativeAbundanceDataframe, # must be a data frame with columns representing features and rows representing samples.
          phenodata=YourPhenotypeDataframe, # must be a data frame with two columns, and colnames are `id` (column 1) and `grp` (column 2). column id represent the sample id, column grp consist of `grp1` and `grp2`, representing the ctrl and disease, repsectively.
          k="euclidean", # distance calculating method. it must bu consist with the method in `dist` function.
          SavePath = NULL, # output directory. Default: ./
          ShuffleWstat = NULL,  # output W stat
          BoundarySample = NULL, # output boundary samples with distance
          BoundaryPair=NULL, # output final pairs 
          ShuffleTime=10000, # shuffle time
          DownPercent = 0.2, # lower percentage of shuffle pairs
          Uppercent=0.8,
          PvalueCutoff = 0.05) # p value cutoff of `Incre.aveRank.P` and `Decre.aveRank.P`
```
### Parameter Description
- data: must be a data frame with columns representing features and rows representing samples.
- phenodata: must be a data frame with two columns, and colnames are `id` (column 1) and `grp` (column 2). column id represent the sample id, column grp consist of `grp1` and `grp2`, representing the ctrl and disease, repsectively.
- k: distance calculating method. it must bu consist with the method in `dist` function.
- SavePath: filename of output directory. Default: ./
- Cut_pair: the cutoff of redundant pairs to perform permutation or wilcox rank paired test. Defult 25
- method_choose: necessary parameter of greated than 25 pair, choose from :"Wilcox","Permutation",
- ShuffleWstat: filename of shuffle W stats
- BoundarySample: filename of output boundary samples with distance
- BoundaryPair: filename of output final pairs 
- ShuffleTime: shuffle time
- DownPercent: lower percentage of shuffle pairs
- Uppercent: higher percentage of shuffle pairs
- PvalueCutoff: p value cutoff of `Incre.aveRank.P` and `Decre.aveRank.P`

**Note:** if the sample pairs are less than 10, it will return nothing.

### Output
- **Decre.aveRank.P** : P value based on averange rank. Generally, to confirm the significant features, we usually use variable `Decre.aveRank.P` to evaluate the disease-enriched features and `Incre.aveRank.P` to evaluate the control enriched features. 
- **Incre.aveRank.P** : P value based on averange rank. Generally, to confirm the significant features, we usually use variable `Decre.aveRank.P` to evaluate the disease-enriched features and `Incre.aveRank.P` to evaluate the control enriched features. 
- Decre.minRank.* : P value based on min rank.
- Incre.minRank.P : P value based on min rank.
- Decre.maxRank.P : P value based on max rank.
- Incre.maxRank.P : P value based on max rank.
- **Species** : Feature
- **Enriched** : phenotype-enriched groups with provided p value cutoff
- Decre.maxRank.P.FDR : P.adjust of Decre.maxRank.P
- Decre.minRank.P.FDR : P.adjust of Decre.minRank.P
- **Decre.aveRank.P.FDR** : P.adjust of Decre.aveRank.P
- **Incre.aveRank.P.FDR** : P.adjust of Incre.aveRank.P
- Decreasing.Rank.Max : shuffle W stats rank based on max method by decreasing
- Increasing.Rank.Max : shuffle W stats rank based on max method by increasing
- Decreasing.Rank.Min : shuffle W stats rank based on min method by decreasing
- Increasing.Rank.Min : shuffle W stats rank based on min method by increasing
- **Decreasing.Rank.Average** : shuffle W stats rank based on average method by decreasing
- **Increasing.Rank.Average** : shuffle W stats rank based on average method by increasing
- **Ctlmean**: mean value of relativa abundance for paired control samples
- **Dismean**: mean value of relativa abundance for paired disease samples

**Important** : Must keep the `Same Order` of samples in  relative abundance dataframe and phenotype dataframe. And we also keep the samples `clustered and placed` according the `grp1 and grp2` of variable `grp` in phenotype data. 

**You can refer to the following example data format.**
  
### Usage -> Mock Data 
  
**if you want to skip the download step, you can find the mock and real datasets in [mengqingren/VTwins.Linux](https://github.com/mengqingren/VTwins.Linux)**
  
``` r
library(tidyverse)
library(VTwins)
library(vegan)
# relative abundance dataframe
set.seed(12345)
dataset <- data.frame(matrix(runif(1200, min = 1e-5, max = 1),nrow = 120,ncol = 10))
colnames(dataset) <- paste("Feature",1:10,sep = '')
rownames(dataset) <- paste("Sample",1:120,sep = '')
dataset.normalized <- decostand(dataset,method = "total",1) #normalization for fetures like species's relative abundance 
#write.table(dataset.normalized,file = "test.data.txt",sep = '\t',quote = F)

# phenotype dataframe
phe_data <- data.frame(id = paste("Sample",1:120,sep = ''),grp=rep(c("grp1","grp2"),c(60,60)))

### dataset.normalized <- read.table("test.data.txt",header = T,row.names = 1,sep = '\t')
### phe_data <- read.table("test.phenodata.txt",header = T,sep = "\t")

# Run
res <- pair_find(data=dataset.normalized,
                 phenodata=phe_data,
                 k="euclidean",
                 Cut_pair=25, 
                 method_choose="Permutation",
                 SavePath = "./",
                 ShuffleWstat = "ShuffleWstat", 
                 BoundarySample = "BoundarySample",
                 BoundaryPair="BoundaryPair",
                 ShuffleTime=10000,
                 DownPercent = 0.2,
                 Uppercent=0.8,PvalueCutoff=0.01)
res$log2FC = log2(as.numeric(res$Dismean)/as.numeric(res$Ctlmean))
write.csv(res,"Results.csv",row.names=F)
```

### Usage -> Real Data

``` r
library(curatedMetagenomicData) #bioconductor
library(phyloseq) #bioconductor
library(tidyverse)

# download dataset -> object
ZellerGData <- curatedMetagenomicData("ZellerG_2014.metaphlan_bugs_list*",
                                      dryrun = FALSE,
                                      counts = TRUE,
                                      bugs.as.phyloseq = TRUE)
psAB <- subset_samples(ZellerGData$ZellerG_2014.metaphlan_bugs_list.stool, study_condition != "adenoma")
psAB <- prune_samples(sample_sums(psAB) >= 10^3, psAB)
CRC_WMS <- filter_taxa(psAB,function(x) sum(x>0)>0,1)

# phenodata frame
grp1_name = "control"
grp2_name = "CRC"
variable_name = "study_condition"
# transitional phenotype dataframe
phe_data <- CRC_WMS@sam_data %>% data.frame() %>% dplyr::select("study_condition") %>% rownames_to_column() %>% dplyr::rename(id=1,grp=2) %>%
  mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2")) %>% arrange(grp)

# transitional relative abundance dataframe
dataset <- t(CRC_WMS@otu_table@.Data)%>% data.frame() %>% rownames_to_column() %>% merge(phe_data,.,by.x="id",by.y="rowname") %>%
  arrange(grp)

# Final phenotype dataframe
phe_data2 <- dataset %>% dplyr::select(id,grp)

# Final relative abundance dataframe
dataset2 <- dataset %>% dplyr::select(-grp) %>% remove_rownames() %>% column_to_rownames("id")

### dataset2 <- readRDS("RealData.dataset2.RDS")
### phe_data2 <- readRDS("RealData.phe_data2.RDS")

# Run
res <- pair_find(data=dataset2,
                 phenodata=phe_data2,
                 k="euclidean",
                 Cut_pair=25, 
                 method_choose="Permutation",
                 SavePath = "./",
                 ShuffleWstat = "RealData.ShuffleWstat", 
                 BoundarySample = "RealData.BoundarySample",
                 BoundaryPair="RealData.BoundaryPair",
                 ShuffleTime=10000,
                 DownPercent = 0.2,
                 Uppercent=0.8,PvalueCutoff=0.01)
res$log2FC = log2(as.numeric(res$Dismean)/as.numeric(res$Ctlmean))
write.csv(res,"RealData.Results.csv",row.names=F)
```
