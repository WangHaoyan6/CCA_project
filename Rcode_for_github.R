setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28")
##library load
{library(singscore);library("scatterplot3d")
  library("CMScaller")
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(data.table)
  library(tidyverse)
  library(reshape2)
  library(finalfit)
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(psych)
  library(pheatmap)
  library(org)
  library(DoubletFinder)
  library(phangorn)
  library(ape)
  #library(SingleR)
  #library(scRNAseq);
  #library(scater);
  #library(xlsx);
  library(ggpubr);
  library(reshape2)
  library(cowplot)
  library(stringi)
  library(stringr)
  `%notin%` <- Negate(`%in%`)
  library(reflect)
  library(ComplexHeatmap)
  library(circlize)
  library("GSVA")
  library("dplyr")
  library("survival")
  library("ggplot2")
  library("survminer")
  library("ggpubr")
  library(limma)
  library("scales")
  library(RColorBrewer)
  library(pheatmap)
  library(msigdbr)
  library(readr)
  library(readxl)
  library("DESeq2")
  library(vsn)
  library(readr)
  library(NMF)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library("GSVA")
  library("dplyr")
  library("survival")
  library("ggplot2")
  library("survminer")
  library("ggpubr")
  library(limma)
  library("scales")
  library(RColorBrewer)
  library(pheatmap)
  library(msigdbr)
  library(readr)
  library(readxl)
  library("DESeq2")
  library(vsn);library(estimate);library(readxl)
  library(GSVA)
  library("apeglm")
  library(readr)
  library(NMF)
  library(vcd)
  library(dplyr)
  library(pheatmap)
  library("GSVA")
  library("dplyr")
  library("survival")
  library("ggplot2")
  library("survminer")
  library("ggpubr")
  library(limma)
  library("scales")
  library(RColorBrewer)
  library(pheatmap)
  library(msigdbr)
  library(readr)
  library(readxl)
  library("DESeq2")
  library(vsn)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library("clusterProfiler")
  library(MASS)
  library(nnet)
  library("apeglm")
  library("dplyr")
  library(ggrepel)
  library(ggplot2)
  library(grid)
  library(ggpubr)
  library(ggthemes)
  library(gmodels)
  library(export)
  library(ggstatsplot)
}

# Mutation
gene_fusion <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/gene_fusion.xlsx",sheet = "positive")

# RNA_sample_filter_ QC合格 + removing_samples
# input
#QC合格 444例
QC <- read_excel("QC.xlsx", sheet = "RNA (2)")
QC<-QC[which(QC$QC=="合格"),]
#474 生存合格
Survival_outcome <- read_excel("Survival_outcome.xlsx")
Survival_outcome<-Survival_outcome[which(Survival_outcome$Sample_excluding!="排除样本" | is.na(Survival_outcome$Sample_excluding)),]
unique(QC$Identity)
Merged_data<-merge(QC,Survival_outcome,all.x = T,all.y = T,by.x="Identity",by.y ="Identity")
Merged_data<-Merged_data[which(!is.na(Merged_data$SampleID)),]
Merged_data<-Merged_data[which(!is.na(Merged_data$Name.y)),]
samples <- read.delim("~/Desktop/GenePlus/project/胆管癌/2022_12_01/1201_WGCNA_group.txt")
samples<-samples$ID
Merged_data<-Merged_data[which(Merged_data$SampleID%in%samples),]

#共440样本
genecount <- read_csv("genecount.csv")

## 439 sample RNA count data
RNA_Sample<-colnames(genecount)
RNA_Sample<-RNA_Sample[-1]
RNA_Sample<-RNA_Sample[RNA_Sample%in%Merged_data$SampleID]
genecount<-genecount[,c("...1",RNA_Sample)]
#write.table(Merged_data,file = "~/Desktop/GenePlus/project/胆管癌/2022_12_01/Merged_data.txt",sep="\t",quote=F)

### Mutation merge
gene_list<-merge(gene_fusion,Merged_data,all.x=T,all.y=F,by.x="sampleID",by.y="TestID")

# 439 sample survival analysis
Survival_analysis<-Merged_data[which(Merged_data$SampleID%in%RNA_Sample),]
Survival_analysis[is.na(Survival_analysis$loss_information),]$"loss_information"<-"No"
Survival_analysis<-Survival_analysis[which(Survival_analysis$loss_information=="No"),]
mydata = Survival_analysis
mydata<-as.data.frame(mydata)
groups=c("histological_type")
for(vals in groups){
  dfDFS = mydata[which(!is.na(mydata[,vals]) & 
                         !is.na(mydata["DFS"]) &
                         mydata["DFS"] !="Unknown" &
                         !is.na(mydata["Relapse"])),]
  dfDFS$Relapse = as.numeric(as.character(dfDFS$Relapse))
  dfDFS$DFS = as.numeric(as.character(dfDFS$DFS))
  dfDFS["group"] = dfDFS[,vals]
  fitDFS <- survfit(Surv(DFS,Relapse) ~ group,data =dfDFS)
  # coxDFS=coxph(Surv(DFS_MONTHS, DFS_Event) ~group+Stage,data = dfDFS)
  #DFS
  ggsurvplot(fitDFS,pval = TRUE, conf.int = FALSE,linetype = 1,
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "jco",
             legend.title="",
             ylab="DFS",
             xlab="Time(DAYS)",
             pval.size= 8,pval.coord=c(6,0.3),
             size=1.2
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1),
                    legend.text = c(14,"bold", "darkblue")
    )
  name1 = paste(vals,"DFS.png",sep="_")
  name2 = paste(vals,"DFS.pdf",sep="_")
  ggsave(name1,width=6.5,height = 4)
  ggsave(name2,width=6.5,height = 4)
  #OS
  dfOS = mydata[which(!is.na(mydata[,vals]) & 
                        !is.na(mydata["Survival_status"]) &
                        !is.na(mydata["OS"])&
                        (mydata["OS"] !="Unknown")),]
  dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
  dfOS$OS = as.numeric(as.character(dfOS$OS))
  dfOS["group"] = dfOS[,vals]
  fitOS <- survfit(Surv(OS,Survival_status) ~ group,data =dfOS) 
  ggsurvplot(fitOS,pval = TRUE, conf.int = FALSE,linetype = 1,
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "jco",
             legend.title="",
             ylab="OS",
             xlab="Time(DAYS)",
             pval.size= 8,pval.coord=c(6,0.3)
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1)
    )
  name1 = paste(vals,"OS.png",sep="_")
  name2 = paste(vals,"OS.pdf",sep="_")
  ggsave(name1,width=6.5,height = 4)
  ggsave(name2,width=6.5,height = 4)
}

############### pre Chapter: sample selection ######################################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/0_Sample selection")
Merged_data$Sum_normal<-NA
Merged_data$Sum_normal<-(Merged_data$Normal_duodenum_percentage+Merged_data$Normal_liver_percentage+Merged_data$Normal_Lymphoid_percentage+Merged_data$Normal_Neuron_percentage+Merged_data$Normal_pancreas_percentage)
Purified_merged_data<-Merged_data[which(Merged_data$Sum_normal==0),]
#共158 纯肿瘤样本
#extrahepatic        hilar   intrahepatic 
#56                   84           18 

#随机抽样,取一半做探索集
#extra_sa<-sample(x=Purified_merged_data[which(Purified_merged_data$histological_type=="extrahepatic"),]$SampleID, size=28, replace = FALSE, prob = NULL)
#hilar_sa<-sample(x=Purified_merged_data[which(Purified_merged_data$histological_type=="hilar"),]$SampleID, size=42, replace = FALSE, prob = NULL) 
#intra_sa<-sample(x=Purified_merged_data[which(Purified_merged_data$histological_type=="intrahepatic"),]$SampleID, size=9, replace = FALSE, prob = NULL) 
#explore_set<-Purified_merged_data[which(Purified_merged_data$SampleID%in%c(extra_sa,hilar_sa,intra_sa)),]

randomly_selected_samples<-sample(x=Merged_data[which(Merged_data$Normal_liver_percentage>50 | Merged_data$Normal_pancreas_percentage>50),]$SampleID, size=40, replace = FALSE, prob = NULL)
#randomly_selected_samples<-sample(x=Merged_data[which(Merged_data$Sum_normal>0),]$SampleID, size=50, replace = FALSE, prob = NULL)

explore_set<-Merged_data[which(Merged_data$SampleID%in%randomly_selected_samples),]
#write.table(explore_set,file="explore_set.txt",sep = "\t",quote=F,row.names = F)
###################### 3 classifier by DGE with limma ############################
#######################  CMScaller for selecting samples
raw_TPM <- read_csv("raw_TPM.csv")
X<-as.data.frame(table(raw_TPM$...1))
raw_TPM_1<-raw_TPM[raw_TPM$...1%in%as.character(X[which(X$Freq==1),1]),]
raw_TPM_2<-raw_TPM[raw_TPM$...1%in%as.character(X[which(X$Freq>1),1]),]
temperate<-matrix(NA,nrow=length(unique(raw_TPM_2$...1)),ncol = ncol(raw_TPM)-1)
for (i in 1:length(unique(raw_TPM_2$...1))) {
  temperate[i,]<-colMeans(raw_TPM_2[which(raw_TPM_2$...1==unique(raw_TPM_2$...1)[i]),-1])
}
temperate<-as.data.frame(temperate)
colnames(temperate)<-colnames(raw_TPM)[-1]
temperate$...1<-unique(raw_TPM_2$...1)
raw_TPM<-rbind(raw_TPM_1,temperate)
raw_TPM_gene<-raw_TPM$...1
raw_TPM<-raw_TPM[,-1]
raw_TPM<-as.data.frame(raw_TPM)
rownames(raw_TPM)<-raw_TPM_gene

#Nearest_template <- read.delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/0_Sample selection/Nearest_template.tsv")
Nearest_template <- read.delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/0_Sample selection/Nearest_template_2.tsv")
Nearest_template<-Nearest_template[which(Nearest_template$Gene%in%rownames(raw_TPM)),]
#Nearest_template<-Nearest_template[which(Nearest_template$sig=="up"),]
#Nearest_template<-Nearest_template[which(abs(Nearest_template$logFC)>3),]
#Nearest_template<-Nearest_template[which(abs(Nearest_template$logFC>2)&Nearest_template$P.Val<0.05),]
table(Nearest_template$Classification)
templates<-data.frame(class=Nearest_template$Classification,probe=Nearest_template$Gene)

##改样本！
#emat<-tpm_exp[,explore_set$SampleID]
#emat<-as.data.frame(emat)
#rownames(emat)<-rownames(tpm_exp)
emat<-raw_TPM[,explore_set$SampleID]
emat<-as.data.frame(emat)
rownames(emat)<-rownames(raw_TPM)
##########################################
emat<-emat[unique(templates$probe),]
rownames(emat)<-unique(templates$probe)
ematAdjust<-ematAdjust(emat = emat,center = T,scale = T)
test<-ntp(
  emat=ematAdjust,
  templates=templates,
  nPerm = 500,
  distance = "cosine",
  nCores = 1,
  seed = NULL,
  verbose = getOption("verbose"),
  doPlot = T)

table(test$FDR<0.05)
test<-test[explore_set$SampleID,]
explore_set$predicted_result<-test$prediction
explore_set$predict_FDR<-test$FDR
explore_set$d.liver<-test$d.Liver
explore_set$d.Pancreas<-test$d.Pancreas
table(explore_set$Normal_liver_percentage,explore_set$predicted_result)
table(explore_set$Normal_pancreas_percentage,explore_set$predicted_result)

ggplot(data=explore_set, aes(y=d.Pancreas, x=Normal_pancreas_percentage)) +
  geom_point(size=1) +
  geom_smooth(method="lm",linetype=1)+
  labs(title="Dist:pancreas vs Pancreas infiltration", 
       x="Normal_pancreas_percentage%", 
       y="Distance of samples regarding to pan_tempelate")+
  stat_regline_equation(label.y = 1,label.x = 5, aes(label = after_stat(..eq.label..))) +stat_regline_equation(label.y = 0.95,label.x = 5, aes(label = ..rr.label..))

ggplot(data=explore_set, aes(y=d.liver, x=Normal_liver_percentage)) +
  geom_point(size=1) +
  geom_smooth(method="lm",linetype=1)+
  labs(title="Dist:liver vs  liver infiltration", 
       x="Normal_liver_percentage%", 
       y="Distance of samples regarding to liver tempelate")+
  stat_regline_equation(label.y = 1,label.x = 5, aes(label = after_stat(..eq.label..))) +stat_regline_equation(label.y = 0.95,label.x = 5, aes(label = ..rr.label..))
########################## 
#总筛:所有样本
emat<-raw_TPM[,Merged_data$SampleID]
emat<-as.data.frame(emat)
rownames(emat)<-rownames(raw_TPM)

##########################################
emat<-emat[unique(templates$probe),]
rownames(emat)<-unique(templates$probe)
ematAdjust<-ematAdjust(emat = emat,center = T,scale = T)
test<-ntp(
  emat=ematAdjust,
  templates=templates,
  nPerm = 500,
  distance = "cosine",
  nCores = 1,
  seed = NULL,
  verbose = getOption("verbose"),
  doPlot = T)
table(test$FDR<0.05)
test$SampleID<-rownames(test)
test<-as.data.frame(test)
Merged_data<-merge(Merged_data,test,by.x="SampleID",by.y="SampleID")
table(Merged_data$Normal_liver_percentage,Merged_data$predicted_result)
table(Merged_data$Normal_pancreas_percentage,Merged_data$predicted_result)

#筛选标准：Sum_normal<20 & FDR>0.1
Merged_data$Sig_Normal_infiltration<-"Significant"
Merged_data[which(Merged_data$Sum_normal<30 & Merged_data$FDR>0.1),]$Sig_Normal_infiltration<-"Not Sig"

ggplot(data=Merged_data, aes(x=FDR, y=Sum_normal,colour=Sig_Normal_infiltration)) +
  geom_point(alpha=0.4,size=3) + 
  geom_hline(yintercept =25,lty=1,col="grey",lwd=1)+
  geom_vline(xintercept=c(0.1),lty=1,col="grey",lwd=1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme_classic(base_size = 16)+
  labs(title = "Normal tissue infiltration vs NTP-FDR for selecting samples unaffected",x="FDR",y="Overall normal tissue infiltration (%)")+scale_colour_manual(values = c("red", "grey"))

table(Merged_data$Sig_Normal_infiltration)
Sample_analysising<-Merged_data[which(Merged_data$Sig_Normal_infiltration=="Not Sig"),]$"SampleID"
table(Merged_data[which(Merged_data$Sig_Normal_infiltration=="Not Sig"),]$"histological")
#extrahepatic        hilar intrahepatic 
#53           84           27 

##############Chapter 1. Concensus clustering ###################################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/")
library(ConsensusClusterPlus)
NMF_mat<-data_for_nmf
NMF_mat<-NMF_mat[,c("...1",Sample_analysising)]

#NMF_mat<-NMF_mat[which(rowSums(NMF_mat[,-1])>0),]
Gene_tmp<-NMF_mat$...1
NMF_mat<-NMF_mat[,-1]
rownames(NMF_mat)<-Gene_tmp
#res.multirun_extraH <- nmf(NMF_mat, rank=2:6, nrun=5)

mads<-apply(NMF_mat, 1, mad)
NMF_mat<-NMF_mat[rev(order(mads)),]
NMF_mat= sweep(NMF_mat,1, apply(NMF_mat,1,median,na.rm=T))
library(ConsensusClusterPlus)

############################### To linux ###############################
#title="~/Desktop/GenePlus/project/胆管癌/2022_12_28/1_NMF/"
#results = ConsensusClusterPlus(as.matrix(NMF_mat),maxK=6,reps=800,pItem=0.8,pFeature=1, title=title,clusterAlg="km",distance="euclidean",plot="png")
#paper:Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data
#In this case, the CDF will be horizontal on the interval [0, 1) as in panel B. In a less-than-ideal case, the CDF will rise more or less gradually to the right, as in panel A. We also want to identify as many clusters as possible without sacrificing the goal that most consensus indices should be near 0 or 1.

load("/Users/zouxinchen/Desktop/GenePlus/project/胆管癌/2022_12_28/1_NMF/Consensus_clustering.RData")
icl = calcICL(results,title=title,plot="png")
icl[["clusterConsensus"]]
k=2
consensus_clustering<-icl[["itemConsensus"]]
consensus_clustering<-consensus_clustering[which(consensus_clustering$k==2),]
consensus_clustering<-consensus_clustering[which(consensus_clustering$itemConsensus>0.5),]

cluster_information<-data.frame(sample_name=consensus_clustering$item,clustering=consensus_clustering$cluster)
cluster_information<-unique(cluster_information)
Merged_data<-merge(Merged_data,cluster_information,by.x = "SampleID",by.y="sample_name",all.x=T,all.y=T)
#write.table(cluster_information,file="~/Desktop/GenePlus/project/胆管癌/2022_12_28/1_NMF/cluster_information.txt",sep="\t",quote=F)
Merged_data<-Merged_data[which(!is.na(Merged_data$clustering)),]

########## Chapter 2. 胆管癌亚型的生存分析 & 不同location 亚型生存分析#########
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/2_survival_outcome/")
#intrahepatic,perihilar,Distal
mydata = Merged_data[which(Merged_data$histological_type=="intrahepatic"),]
mydata<-as.data.frame(mydata)
groups=c("clustering")
for(vals in groups){
  dfDFS = mydata[which(!is.na(mydata[,vals]) & 
                         mydata$SampleID%notin% mydata[which(mydata$Death_reason==" 0"),]$SampleID&
                         !is.na(mydata["DFS"]) &
                         mydata["DFS"] !="Unknown" &
                         !is.na(mydata["Relapse"])),]
  dfDFS$Relapse = as.numeric(as.character(dfDFS$Relapse))
  dfDFS$DFS = as.numeric(as.character(dfDFS$DFS))
  dfDFS["group"] = dfDFS[,vals]
  fitDFS <- survfit(Surv(DFS,Relapse) ~ group,data =dfDFS)
  # coxDFS=coxph(Surv(DFS_MONTHS, DFS_Event) ~group+Stage,data = dfDFS)
  #DFS
  ggsurvplot(fitDFS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="DFS",
             xlab="Time(DAYS)",
             pval.size= 8,pval.coord=c(6,0.3),
             size=1.2
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1),
                    legend.text = c(14,"bold", "darkblue")
    )
  name2 = paste(vals,"DFS.pdf",sep="_")
  ggsave(name2,width=12,height = 8)
  #OS
  dfOS = mydata[which(!is.na(mydata[,vals]) &  
                        mydata$SampleID%notin% mydata[which(mydata$Death_reason==" 0"),]$SampleID &
                        !is.na(mydata["Survival_status"]) &
                        !is.na(mydata["OS"])&
                        (mydata["OS"] !="Unknown")),]
  dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
  dfOS$OS = as.numeric(as.character(dfOS$OS))
  dfOS["group"] = dfOS[,vals]
  fitOS <- survfit(Surv(OS,Survival_status) ~ group,data =dfOS) 
  ggsurvplot(fitOS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = F,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="OS",
             xlab="Time(DAYS)",
             pval.size= 8,pval.coord=c(6,0.3)
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1)
    )
  name2 = paste(vals,"OS.pdf",sep="_")
  ggsave(name2,width=12,height = 8)
}


################################# Chapter 3. PCA plot  #######
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/3_PCA_plot/")
tpm_exp<-data_for_nmf
tpm_exp<-tpm_exp[,c("...1",cluster_information[order(cluster_information$clustering),]$sample_name)]
gene_list<-tpm_exp$...1
tpm_exp<-tpm_exp[,-1]
rownames(tpm_exp)<-gene_list
pca.info <- prcomp(tpm_exp, scale = T,
                   center = TRUE, retx = T)
screeplot(pca.info, npcs = 10, type = "lines")
head(pca.info)
summary(pca.info) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x) #样本得分score
rownames(Merged_data)<-Merged_data$SampleID
pca.data<-data.frame(sample=rownames(pca.info$rotation),Type=as.character(c(Merged_data[rownames(pca.info$rotation),]$"clustering")), pca.info$rotation)

ggscatter(data = pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, 
          size=2,repel=T, main="PCA plot"+theme_base(),palette="ucscgb")
ggsave(filename = "./molecular_subtype_PCA.png")

pca.data<-data.frame(sample=rownames(pca.info$rotation),Type=as.character(c(Merged_data[rownames(pca.info$rotation),]$"histological_type")), pca.info$rotation)

ggscatter(data = pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, 
          size=2,repel=T, main="PCA plot"+theme_base(),palette="ucscgb")
ggsave(filename = "./histological_PCA.png")


#####################  Chapter 4. Hallmarker score for class #################

setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/4_Hallmarker_GSEA/")
all_gene_sets<-msigdbr(species = "Homo sapiens")
all_gene_sets<-as.data.frame(all_gene_sets)
Hallmarker_gene_set<-all_gene_sets[which(all_gene_sets$gs_cat=="H"),]
length(table(Hallmarker_gene_set$gs_description))
Hallmarker_gene_set
geneset_description<-unique(Hallmarker_gene_set$gs_description)
geneset_description
gene_list<-rownames(tpm_exp)
gene_list<-gene_list[gene_list%in%Hallmarker_gene_set$gene_symbol]
sample_hallmarker_expre<-tpm_exp[gene_list,]
rownames(sample_hallmarker_expre)<-gene_list
Hallmarker_gene_set<-Hallmarker_gene_set[which(Hallmarker_gene_set$gene_symbol %in% row.names(sample_hallmarker_expre)),]
Hallmarker_gene_set<-Hallmarker_gene_set[,c(4,15)]
Hallmarker_gene_set<-as.data.frame(Hallmarker_gene_set)
mat<-sample_hallmarker_expre
gen<-Hallmarker_gene_set
head(gen)
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="ssgsea", verbose=FALSE, parallel.sz=2)

######## a, Overall distribution #####################
#standarization by z-score
GSEA_sd<-rowSds(es)
es<-as.data.frame(es)
GSEA_mean<-rowMeans(es)
es$GSEA_mean<-GSEA_mean
es$GSEA_sd<-GSEA_sd
es[,-c(165,166)]<-(es[,-c(165,166)]-es$GSEA_mean)/es$GSEA_sd
geneset_description <- read.csv("~/Desktop/GenePlus/project/胆管癌/2022_12_28/4_Hallmarker_GSEA//geneset_description.csv", row.names=1)
geneset_description<-geneset_description[rownames(es),]
rownames(es)<-geneset_description$function.
final_enriched_scores<-as.data.frame(matrix(NA,nrow=50))

for (i in 1:2) {
  tmp_sample<-Merged_data[which(Merged_data$clustering==i),]$SampleID
  Mean_cluster_score<-rowMeans(as.matrix((es[,tmp_sample])))
  final_enriched_scores[,i]<-Mean_cluster_score
}
rownames(final_enriched_scores)<-rownames(es)
colnames(final_enriched_scores)<-c("cluster1","cluster2")
colnames(final_enriched_scores)<-c("mesenchymal and immunosupressive","Metabolic and proliferative")
#final_enriched_scores<-final_enriched_scores[,c(1:2)]
bk<-c(seq(min(final_enriched_scores),0,length.out = 50),seq(0.00001,max(final_enriched_scores),length.out = 50))
###先z-score 后 scale
pdf(file="./Hallmark_GSEA.pdf",width = 8,height = 10)
pheatmap::pheatmap(as.matrix(final_enriched_scores),cluster_cols = F,cluster_rows = T,color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk,main="pheatmap of hallmarker enrichment",scale ="row")
dev.off()

###comment: 
Merged_data$Anno_cluster<-NA
Merged_data[which(Merged_data$clustering==1),]$Anno_cluster<-"mesenchymal and immunosupressive"
Merged_data[which(Merged_data$clustering==2),]$Anno_cluster<-"Metabolic and proliferative"

#additional exploration:
#C1 免疫抑制+间质类
#C2 代谢+分泌

################### Chapter 5、 免疫浸润，CibersortX ################################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/5_immune_score/")
CIBERSORT_score <- read_delim("CIBERSORT_score.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
CIBERSORT_score<-as.data.frame(CIBERSORT_score)
rownames(CIBERSORT_score)<-CIBERSORT_score$Mixture
CIBERSORT_score<-CIBERSORT_score[Merged_data$SampleID,c(1:10,12:23)]
CIBERSORT_score$Anno_cluster<-Merged_data$Anno_cluster
#2:22
for (i in 2:22) {
  A<-ggboxplot(CIBERSORT_score,x = "Anno_cluster", y = colnames(CIBERSORT_score)[i], color = "Anno_cluster",ylab = paste0(colnames(CIBERSORT_score)[i]," proportion", xlab = "Molecular subtype")+stat_summary(fun="mean",color="black"))+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
    stat_compare_means(label.y = 0.5)+
    ggtitle(label=paste0("CIBERSORT result: ",colnames(CIBERSORT_score)[i]))
  
  pdf(file = paste0("./immune cell/",colnames(CIBERSORT_score)[i],".pdf"),width = 10,height=8)
  print(A)
  dev.off()
}


##ESTIMATE: stromal and immune score
library(tidyestimate)
data_for_nmf <- read_csv("~/Desktop/GenePlus/project/胆管癌/2022_12_01/data_for_nmf.csv")
tpm_exp<-data_for_nmf[,-1]
rownames(tpm_exp)<-data_for_nmf$...1
scores <-filter_common_genes(df =tpm_exp, id = "hgnc_symbol", tidy = FALSE, tell_missing = T, find_alias = F)
scores<-estimate_score(df = scores,is_affymetrix = F)
Merged_data<-merge(Merged_data,scores,all.x = T,all.y = T,by.x = "SampleID",by.y="sample")
Merged_data<-Merged_data[which(!is.na(Merged_data$clustering)),]


pdf(file = paste0("./immune cell/stromal score",".pdf"),width = 10,height=8)
A=ggboxplot(Merged_data, x = "Anno_cluster",y = "stromal",color = "Anno_cluster",
            ylab = "Stromal score", xlab = "Molecular subtype")+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means(label.y = 2000)+
  ggtitle(label="ESTIMATE:stromal score vs molecular subtype")
print(A)
dev.off()


pdf(file = paste0("./immune cell/immune score",".pdf"),width = 10,height=8)
A<-ggboxplot(Merged_data, 
             x = "Anno_cluster", y = "immune", 
             color = "Anno_cluster", 
             ylab = "Immune score", xlab = "Molecular subtype")+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means(label.y = 3000)+
  ggtitle(label="ESTIMATE:immune score vs molecular subtype")
print(A)
dev.off()

ggboxplot(Merged_data, 
          x = "Anno_cluster", y = "estimate.y", 
          color = "Anno_cluster", 
          ylab = "estimate score", xlab = "Anno_cluster")+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means(label.y = 3000)+ggtitle(label="ESTIMATE:estimate score vs molecular subtype")
ggsave("./ESTIMATE_SCORE.png",width=8,height=8)

######## TME_exploration ###########################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/5_immune_score/TME_exploration/")
TME_subtype_Exploration <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/TME_subtype_Exploration.xlsx",sheet = "Gene Description")
gen<-data.frame(gene_symbol=TME_subtype_Exploration$Gene,gs_description=TME_subtype_Exploration$`Gene signature`,gs_subtype=TME_subtype_Exploration$TME_subtype)
gene_list<-rownames(raw_TPM)
gene_list<-gene_list[gene_list%in%gen$gene_symbol]
sample_hallmarker_expre<-raw_TPM[gene_list,Merged_data$SampleID]
rownames(sample_hallmarker_expre)<-gene_list
gen<-gen[which(gen$gene_symbol %in% row.names(sample_hallmarker_expre)),]
gen<-as.data.frame(gen)
mat<-sample_hallmarker_expre
head(gen)
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="ssgsea",verbose=FALSE, parallel.sz=2)

results_pancancer <- read_delim("MFP-master/results_pancancer.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
cluster_information<-data.frame(sample_name=results_pancancer$...1,TME_clustering=results_pancancer$MFP)
cluster_information<-unique(cluster_information)
Merged_data<-merge(Merged_data,cluster_information,by.x = "SampleID",by.y="sample_name",all.x=T,all.y=F)

es<-t(scale(t(es),scale=TRUE, center=TRUE))
es<-as.data.frame(es)
Merged_data<-Merged_data[order(Merged_data$Anno_cluster,Merged_data$TME_clustering.y,decreasing=T),]
annotation_col<-data.frame(molecular_cluster=Merged_data$Anno_cluster,TME_clustering=as.character(Merged_data$TME_clustering.y))
rownames(annotation_col)<-Merged_data$SampleID
annotation_row<-data.frame(Category=gen$gs_subtype,Name=gen$gs_description)
annotation_row<-unique(annotation_row)
rownames_annotation_row<-annotation_row$Name
annotation_row<-data.frame(Category=annotation_row$Category)
rownames(annotation_row)<-rownames_annotation_row
es<-es[rownames(annotation_row),Merged_data$SampleID]
bk<-c(seq(min(es),0,length.out = 100),seq(0.00001,max(es),length.out = 100))
###先z-score 后 scale
library(ComplexHeatmap)

columnAnnotation<-columnAnnotation(
  molecular_cluster=annotation_col$molecular_cluster,
  TME_subtype=annotation_col$TME_clustering,
  col =list(molecular_cluster = c("mesenchymal and immunosupressive" = "red", 
                                  "Metabolic and priliferative" ="blue"),
            TME_subtype= c("IE"="royalblue","IE/F"="blue","F"="purple","D"="gold")))

rowAnnotation<-rowAnnotation(pathways=factor(annotation_row$Category,levels = c("Angiogenesis_Fibroblasts","Pro-tumor_Immune infiltration","Anti-tumor_Immune infilatrate","EMT Signature_proliferation rate")))
pdf(file="./MFP-master/TME_subtype_164.pdf",width = 20,height = 12)

Heatmap(as.matrix(es),
        name="TME_subtype vs signatures (z-score of ssgsea)",
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        column_gap=unit(2, "mm"),
        row_gap=unit(2, "mm"),
        right_annotation = rowAnnotation,
        row_title_gp = gpar(fontsize = 10),
        row_title_side = "right",
        row_title_rot = 0,
        row_names_side = "left",
        top_annotation = columnAnnotation,
        show_row_dend = F,
        row_split = annotation_row$Category,
        column_split = annotation_col$TME_clustering)
dev.off()


########## sankey plot
sankey_input<-Merged_data[,c('TME_clustering.y','histological_type','Anno_cluster')]
mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)
#格式转换
UCB_lodes <- to_lodes_form(sankey_input,
                           axes = 1:3,
                           id = "Cohort")
ggplot(UCB_lodes,
       aes(x = x, 
           stratum = stratum, 
           alluvium = Cohort,
           fill = stratum, 
           label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/8,aes.flow = "backward") +
  geom_stratum(alpha = .9,width = 1/8) +
  geom_text(stat = "stratum", size = 3.3,color="black",angle="0") +
  scale_fill_discrete()+
  xlab("") + ylab("") +
  theme_bw() + coord_flip()+
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  theme(axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(face = "bold",color = "black"),
        axis.text.y = element_text(face = "bold",color = "black"))+ #显示分组名字
  ggtitle("") +
  guides(fill = FALSE)
ggsave("./MFP-master/sankey3.pdf",width = 10,height = 10)


A<-as.data.frame(table(Merged_data$Anno_cluster,Merged_data$TME_clustering.y))
A<-A[order(A$Var1),]
A$proportion<-NA
A[c(1:4),]$proportion<-A[c(1:4),]$Freq/sum(A[c(1:4),]$Freq)
A[c(5:8),]$proportion<-A[c(5:8),]$Freq/sum(A[c(5:8),]$Freq)
colnames(A)<-c("Molecular subtype","TME subtype","Freq","Proportion")
# Create a new data frame with the sum of the proportions for each TME subtype

#A<-as.data.frame(table(Merged_data$Anno_cluster,Merged_data$TME_clustering.y))
#A_chisq<-chisq.test(A)
#A_chisq$p.value

ggplot(data = A,aes(x=`Molecular subtype`,y=Proportion,fill =`TME subtype`))+
  geom_col(position = "stack")+
  labs(x = NULL,y = "Fraction (%)")+  theme_classic()+
  scale_fill_manual(values=c("gold",
                             "#9933FF",
                             "royalblue",
                             "darkblue"))+
  theme(axis.text.x = element_text(colour = "black",size = rel(1.7)))+
  theme(axis.text.y = element_text(colour = "black",size = rel(1.7)))+
  theme(axis.title.y = element_text(size=17))+
  ggtitle(subtitle ="Chi-square test:p : 0.03122",label = "Distribution of TME subtypes in each molecular subtype")
ggsave("./MFP-master/TME_subtype_Distribution.pdf",width = 10,height = 8)



############# Chapter 6、 每组免疫治疗相关RNA表达  #############################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/6_Actionable_RNA_target/")
RNA_actionable <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/RNA_changed_actionable_drug.xlsx")
RNA_Act<-RNA_actionable$GeneID
RNA_matrix_actionable<-tpm_exp[RNA_Act,]
rownames(RNA_matrix_actionable)<-RNA_Act
RNA_matrix_actionable<-t(RNA_matrix_actionable)
RNA_matrix_actionable<-RNA_matrix_actionable[Merged_data[which(!is.na(Merged_data$Anno_cluster)),]$"SampleID",]
#RNA_matrix_actionable<-log2(RNA_matrix_actionable+1)
Merged_data_tmp<-cbind(Merged_data[which(!is.na(Merged_data$Anno_cluster)),],RNA_matrix_actionable)

Groups<-unique(Merged_data_tmp$Anno_cluster)
Genes_actionable<-colnames(Merged_data_tmp)[-c(1:51)]
Data_frame_tmp<-as.data.frame(matrix(data = NA,nrow = 89,ncol = 2))
colnames(Data_frame_tmp)<-Groups
rownames(Data_frame_tmp)<-Genes_actionable
Merged_data_tmp<-Merged_data_tmp[,-c(1:29,31:51)]

for (i in unique(Merged_data_tmp$Anno_cluster)) {
  X<-colMeans(Merged_data_tmp[which(Merged_data_tmp$Anno_cluster==i),-1])
  Data_frame_tmp[,i]<-X
}
pdf(file = "actionable_RNA.pdf",width = 10,height=15)
pheatmap(Data_frame_tmp,cluster_cols = T,cluster_rows = T,main="Heatmap of 89 actionable target RNA expression",scale = "row")
dev.off()
##comment: 预后倾向更差的分子分型 actionable target 表达越高！！！
################# Oncoprint for drug prediction #############################
###linux 跑 Rscript 
#### predicts a phenotype (drug sensitivity score) when provided with microarray or bulk RNAseq gene expression data of different platforms
Drug_Prediction <- read_csv("calcPhenotype_Output/DrugPredictions.csv")
Drug_Prediction<-as.data.frame(Drug_Prediction)
rownames(Drug_Prediction)<-Drug_Prediction$...1
Drug_Prediction<-Drug_Prediction[,-1]
Drug_Prediction<-Drug_Prediction[Merged_data$SampleID,]

statistical_result<-data.frame(drug=colnames(Drug_Prediction),t_test_p=NA)
Drug_Prediction$molecular_cluster<-Merged_data$Anno_cluster
Drug_Prediction<-Drug_Prediction[,c(statistical_result$drug,"molecular_cluster")]

for (i in 1:nrow(statistical_result)) {
  X<-t.test(Drug_Prediction[,i]~molecular_cluster, data = Drug_Prediction)
  statistical_result[i,]$t_test_p<-X$p.value
}
table(statistical_result$t_test_p<0.05)
##筛出89个2组显著差异sensitivity的drug
Drug_Prediction<-Drug_Prediction[,statistical_result[which(statistical_result$t_test_p<0.05),]$drug]
Drug_Prediction<-t(Drug_Prediction)

Drug_Prediction<-log2(Drug_Prediction)

Drug_Prediction<-t(scale(t(Drug_Prediction),scale=T,center = T))
anno
anno_col<-data.frame(Molecular_subtype=Merged_data$Anno_cluster)
rownames(anno_col)<-Merged_data$SampleID
pheatmap(Drug_Prediction,show_colnames =F,cluster_cols = F,annotation_col = anno_col)



screened_compunds <- read_csv("DataFiles/DataFiles/GLDS/GDSCv2/screened_compunds_rel_8.2.csv")
annotation_row<-data.frame(drug_name=screened_compunds$DRUG_NAME,target=screened_compunds$TARGET_PATHWAY,targeted_DNA=screened_compunds$TARGET)
annotation_row<-unique(annotation_row)
match_drug2pathway <- read_delim("calcPhenotype_Output/match_drug2pathway.txt", 
                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)
match_drug2pathway$target_name<-NA
annotation_row<-annotation_row[which(annotation_row$drug_name%in%match_drug2pathway$Mathced),]
annotation_row<-annotation_row[which(annotation_row$drug_name!="Dasatinib" | annotation_row$target!="Other, kinases"),]
for (i in 1:nrow(match_drug2pathway)) {
  for (l in 1:nrow(annotation_row)) {
    if(match_drug2pathway[i,]$Mathced==annotation_row[l,]$drug_name){
      match_drug2pathway[i,]$target<-annotation_row[l,]$target
      match_drug2pathway[i,]$target_name<-annotation_row[l,]$targeted_DNA
    }
  }
}

match_drug2pathway<-as.data.frame(match_drug2pathway)
rownames(match_drug2pathway)<-match_drug2pathway$drug_name
library(ggalluvial)
df<-match_drug2pathway[,-c(1,2)]
mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)
UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                           axes = 1:ncol(df),
                           id = "Cohort")
{ggplot(UCB_lodes,
        aes(x = x, stratum = stratum, alluvium = Cohort,
            fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0, 0.2)) + 
    geom_flow(width = 1/10) + #线跟方块间空隙的宽窄
    geom_stratum(alpha = .9,width = 1/8) + #方块的透明度、宽度
    geom_text(stat = "stratum", size = 2,color="black") + 
    xlab("") + ylab("") +
    theme_bw() + #去除背景色
    theme(panel.grid =element_blank()) + #去除网格线
    theme(panel.border = element_blank()) + #去除外层边框
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
    ggtitle("")+
    guides(fill = FALSE) 
}
annotation_row<-data.frame(Target=match_drug2pathway$target)
rownames(annotation_row)<-match_drug2pathway$drug_name
Drug_Prediction <- read_csv("calcPhenotype_Output/DrugPredictions.csv")
Drug_Prediction<-as.data.frame(Drug_Prediction)
rownames(Drug_Prediction)<-Drug_Prediction$...1
Drug_Prediction<-Drug_Prediction[,-1]
Drug_Prediction<-Drug_Prediction[Merged_data$SampleID,]
Drug_Prediction<-t(Drug_Prediction)
Drug_Prediction<-as.data.frame(Drug_Prediction)
Drug_Prediction<-Drug_Prediction[statistical_result[which(statistical_result$t_test_p<0.05),]$drug,]
GSEA_sd<-rowSds(as.matrix(Drug_Prediction))
GSEA_mean<-rowMeans(Drug_Prediction)
Drug_Prediction$GSEA_mean<-GSEA_mean
Drug_Prediction$GSEA_sd<-GSEA_sd
Drug_Prediction[,-c(165,166)]<-(Drug_Prediction[,-c(165,166)]-Drug_Prediction$GSEA_mean)/Drug_Prediction$GSEA_sd
final_enriched_scores<-as.data.frame(matrix(NA,nrow=length(statistical_result[which(statistical_result$t_test_p<0.05),]$drug)))
for (i in 1:2) {
  tmp_sample<-Merged_data[which(Merged_data$clustering==i),]$SampleID
  Mean_cluster_score<-rowMeans(as.matrix((Drug_Prediction[,tmp_sample])))
  final_enriched_scores[,i]<-Mean_cluster_score
}
rownames(final_enriched_scores)<-rownames(Drug_Prediction)
colnames(final_enriched_scores)<-c("mesenchymal and immunosupressive","Metabolic and proliferative")
bk<-c(seq(min(final_enriched_scores),0,length.out = 50),seq(0.00001,max(final_enriched_scores),length.out = 50))
X<-pheatmap::pheatmap(as.matrix(final_enriched_scores),cluster_cols = T,cluster_rows = T,main="Heatmap of Mean(z-scored drug sensitivity) for molecular clusters",color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk,scale ="row",annotation_row = annotation_row)

summary_target<-final_enriched_scores[X$tree_row$order,]
summary_target$enriched_group<-"mesenchymal and immunosupressive"
summary_target[which(summary_target$`mesenchymal and immunosupressive`<summary_target$`Metabolic and secretory`),]$enriched_group<-"metabolic and proliferative"
annotation_row<-annotation_row[rownames(summary_target),]
summary_target$target<-annotation_row
X<-as.data.frame(table(summary_target$enriched_group,summary_target$target))
#write.table(summary_target,file = "./target_supplementary_table.txt",sep="\t",quote=F,row.names = T)
#write.table(X,file = "./target_summary.txt",sep="\t",quote=F,row.names = F)
#X<-X[which(X$Var1=="mesenchymal and immunosupressive"),]
X<-X[which(X$Var1=="metabolic and proliferative"),]
X<-X[which(X$Freq!=0),]
X<-X[order(X$Freq,decreasing=F),]
X$Var2<-factor(X$Var2,levels = X$Var2)

ggplot(X,aes(x=Var2,y=Freq,fill=Freq))+geom_col() +
  # geom_bar(stat = "identity",position=position_dodge())+
  xlab("Target of drug")+ ylab("drug count")+
  theme(axis.text.x =element_text(face="bold",color="Black", size=10, angle=25)
        ,axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=8))+
  ggtitle("Distribution of drug target:metabolic and proliferative")+ coord_flip()+viridis::scale_fill_viridis(direction = -1)
ggsave(filename = "metabolic and proliferative.png",width= 8,height = 6)




######################## Chapter 7、每组的病理学特征 #############################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/7_Pathological_features/")
Merged_data <- Merged_data[which(Merged_data$SampleID%in%samples),]
Merged_data<-Merged_data[order(Merged_data$Anno_cluster),]
Merged_data<-as.matrix(Merged_data)
Merged_data[Merged_data=="NA"]<-NA
Merged_data<-as.data.frame(Merged_data)
unique(Merged_data$histological_type)
Merged_data[which(Merged_data$histological_type=="extrahepatic"),]$histological_type<-"Distal"
Merged_data[which(Merged_data$histological_type=="hilar"),]$histological_type<-"perihilar"

ann<-data.frame(
  Anno_cluster=Merged_data$Anno_cluster,
  histological_type2=Merged_data$histological_type,
  Stage=Merged_data$Stage,
  Differentiation=Merged_data$Differentiation,
  Pathological=Merged_data$Pathological,
  Sex=Merged_data$Sex,
  Resection=Merged_data$Resection,
  Normal_liver_percentage=as.numeric(Merged_data$Normal_liver_percentage),
  Normal_duodenum_percentage=as.numeric(Merged_data$Normal_duodenum_percentage),
  Normal_Lymphoid_percentage=as.numeric(Merged_data$Normal_Lymphoid_percentage),
  Normal_Neuron_percentage=as.numeric(Merged_data$Normal_Neuron_percentage),
  Normal_pancreas_percentage=as.numeric(Merged_data$Normal_pancreas_percentage),
  stringsAsFactors=F)

#type: Allowed values are one of: “div”, “qual”, “seq”, or “all”
#Seq: Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu YlOrBr, YlOrRd.
#Qual: Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3.
#Diverging: BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
#display.brewer.pal(n = 8, name = "Pastel2")
#brewer.pal(n = 8, name = "Pastel2")

#colors
{
  colours <- list(
    Stage=c('I'= "#FFF7FB", 'IA'="#ECE2F0", 'IB'="#D0D1E6",'II'="#A6BDDB",'IIA'="#67A9CF",'IIB'="#3690C0",'IIIA'="#2171B5",'IIIB'="#084594",'IIIC'= '#02818A', 'IV'= "#016450",'IVA'="#016450" ,'IVB'="#016450"),
    Sex=c('Male'='#386CB0','female'='#F0027F'),
    Resection = c('0' = '#F7FCF5', '1' = '#C7E9C0', '2' = '#74C476', '3' = '#238B45', '4' = '#00441B'),
    histological_type2=c('Distal' = '#7FC97F', 'perihilar' = '#BEAED4', 'intrahepatic'='#FDC086'),
    Differentiation=c('Low'= '#00441B','Medium_low'='#238B45','Medium'='#FB8072','High_medium'='#D94801','High'='#7F2704'),
    Pathological=c('adenocarcinoma'= '#8DD3C7', 'hepatocellular carcinoma'='#FFFFB3' ,'Lymphoepithelioma'= '#BEBADA', 'Mixed_adenocarcinoma_and_squamousCarcinoma'='#FB8072' ,'mucinous adenocarcinoma'='#80B1D3' , 'Other'= '#FDB462' ,'sarcomatoid carcinoma'='#D9D9D9'),Anno_cluster=c('mesenchymal and immunosupressive'='#7FC97F',
                                                                                                                                                                                                                                                                                          'Metabolic and secretory'='#386CB0'))
}
#plot
{
  colAnn <- HeatmapAnnotation(
    df=ann,
    simple_anno_size = unit(8, "mm"),
    na_col = 'white',
    col = colours,
    annotation_height = 2,
    annotation_width = unit(0.5, 'mm'),
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      Anno_cluster=list(nrow = 4, # number of rows across which the legend will be arranged
                        title = 'Anno_cluster',
                        title_position = 'topcenter',
                        legend_direction = 'vertical',
                        title_gp = gpar(fontsize = 6, fontface = 'bold'),
                        labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Sex=list(nrow = 2, # number of rows across which the legend will be arranged
               title = 'Sex',
               title_position = 'topcenter',
               legend_direction = 'vertical',
               title_gp = gpar(fontsize = 6, fontface = 'bold'),
               labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Resection=list(nrow = 5, # number of rows across which the legend will be arranged
                     title = 'Resection',
                     title_position = 'topcenter',
                     legend_direction = 'vertical',
                     title_gp = gpar(fontsize = 6, fontface = 'bold'),
                     labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      histological_type2=list(nrow = 3,# number of rows across which the legend will be arranged
                              title = 'histological_type2',
                              title_position = 'topcenter',
                              legend_direction = 'vertical',
                              title_gp = gpar(fontsize = 6, fontface = 'bold'),
                              labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Stage=list(nrow = 12, # number of rows across which the legend will be arranged
                 title = 'Stage',
                 title_position = 'topcenter',
                 legend_direction = 'vertical',
                 title_gp = gpar(fontsize = 6, fontface = 'bold'),
                 labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Differentiation=list(nrow = 5, 
                           title = 'Differentiation',
                           title_position = 'topcenter',
                           legend_direction = 'vertical',
                           title_gp = gpar(fontsize = 6, fontface = 'bold'),
                           labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Pathological=list(nrow = 7, # number of rows across which the legend will be arranged
                        title = 'Cancer type',
                        title_position = 'topcenter',
                        legend_direction = 'vertical',
                        title_gp = gpar(fontsize = 6, fontface = 'bold'),
                        labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Normal_liver_percentage=list(nrow = 20, # number of rows across which the legend will be arranged
                                   title = 'liver %',
                                   title_position = 'topcenter',
                                   legend_direction = 'vertical',
                                   title_gp = gpar(fontsize = 6, fontface = 'bold'),
                                   labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Normal_duodenum_percentage=list(nrow = 10, # number of rows across which the legend will be arranged
                                      title = 'duodenum %',
                                      title_position = 'topcenter',
                                      legend_direction = 'vertical',
                                      title_gp = gpar(fontsize = 6, fontface = 'bold'),
                                      labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Normal_Lymphoid_percentage=list(nrow = 20, # number of rows across which the legend will be arranged
                                      title = 'lymphoid %',
                                      title_position = 'topcenter',
                                      legend_direction = 'vertical',
                                      title_gp = gpar(fontsize = 6, fontface = 'bold'),
                                      labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Normal_Neuron_percentage=list(nrow = 20, # number of rows across which the legend will be arranged
                                    title = 'neuron %',
                                    title_position = 'topcenter',
                                    legend_direction = 'vertical',
                                    title_gp = gpar(fontsize = 6, fontface = 'bold'),
                                    labels_gp = gpar(fontsize = 6, fontface = 'bold')),
      Normal_pancreas_percentage=list(nrow = 20, # number of rows across which the legend will be arranged
                                      title = 'pancreas %',
                                      title_position = 'topcenter',
                                      legend_direction = 'vertical',
                                      title_gp = gpar(fontsize = 6, fontface = 'bold'),
                                      labels_gp = gpar(fontsize = 6, fontface = 'bold')) 
    ))
}
#
pdf(file = "./pathological_features.pdf",width = 20,height=15)
plot(colAnn)
dev.off()


table(Merged_data$Sex)
summary(Merged_data$Age)
table(Merged_data$Resection)
table(Merged_data$histological_type)

table(Merged_data$Stage)
table(Merged_data$Stage)/438*100

table(Merged_data$Differentiation)
table(Merged_data$Differentiation)/438*100

table(Merged_data$Pathological)
table(Merged_data$Pathological)/438*100
table(Merged_data$Normal_liver_percentage)
table(Merged_data$Normal_pancreas_percentage)
table(Merged_data$Normal_duodenum_percentage)
table(Merged_data$Normal_Lymphoid_percentage)
table(Merged_data$Normal_Neuron_percentage)
selectedsampleID<-(Merged_data$SampleID)




###p value
######### categorical vs continuous #########################################
table(Merged_data$Anno_cluster,Merged_data$Age)
table(Merged_data$Anno_cluster,Merged_data$Normal_liver_percentage)
table(Merged_data$Anno_cluster,Merged_data$Normal_pancreas_percentage)
table(Merged_data$Anno_cluster,Merged_data$Normal_duodenum_percentage)
table(Merged_data$Anno_cluster,Merged_data$Normal_Lymphoid_percentage)
table(Merged_data$Anno_cluster,Merged_data$Normal_Neuron_percentage)
Merged_data$Age<-as.numeric(Merged_data$Age)
Merged_data$Normal_liver_percentage<-as.numeric(Merged_data$Normal_liver_percentage)
Merged_data$Normal_pancreas_percentage<-as.numeric(Merged_data$Normal_pancreas_percentage)
Merged_data$Normal_duodenum_percentage<-as.numeric(Merged_data$Normal_duodenum_percentage)
Merged_data$Normal_Lymphoid_percentage<-as.numeric(Merged_data$Normal_Lymphoid_percentage)
Merged_data$Normal_Neuron_percentage<-as.numeric(Merged_data$Normal_Neuron_percentage)

catVScon<-c("Age","Normal_liver_percentage","Normal_pancreas_percentage","Normal_duodenum_percentage","Normal_Lymphoid_percentage","Normal_Neuron_percentage")
Merged_data[which(is.na(Merged_data$Normal_liver_percentage)),]$Normal_liver_percentage<-as.numeric(0)
Merged_data[which(is.na(Merged_data$Normal_pancreas_percentage)),]$Normal_pancreas_percentage<-as.numeric(0)
Merged_data[which(is.na(Merged_data$Normal_duodenum_percentage)),]$Normal_duodenum_percentage<-as.numeric(0)
Merged_data[which(is.na(Merged_data$Normal_Lymphoid_percentage)),]$Normal_Lymphoid_percentage<-as.numeric(0)
Merged_data[which(is.na(Merged_data$Normal_Neuron_percentage)),]$Normal_Neuron_percentage<-as.numeric(0)


for (i in catVScon) {
  ggboxplot(Merged_data, x = "Anno_cluster", y = i, 
            color = "Anno_cluster", 
            ylab = i, xlab = "Anno_cluster")+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
    stat_compare_means(label.y = 100)+ggtitle(label =paste0(i," vs molecular subtype"))
  ggsave(paste0("./7_Pathological_features/",i,".png"),width=8,height=8)
}


######### categorical vs categorical #################################
X<-table(Merged_data$Anno_cluster,Merged_data$Sex)
X<-table(Merged_data$Anno_cluster,Merged_data$Resection)
X<-table(Merged_data$Anno_cluster,Merged_data$Stage)
X<-table(Merged_data$Anno_cluster,Merged_data$Differentiation)
X<-table(Merged_data$Anno_cluster,Merged_data$Pathological)

test <- fisher.test(X,alternative = "two.sided",simulate.p.value = T,B=1e6)
mosaicplot(X,
           main = "Mosaic plot",
           color = T,
           xlab = "Molecular subtype",ylab = "Pathological",
           sub = paste0("two-sided Fisher's Exact Test with simulated p-value:",test$p.value,sep=" "))

####################### Driver Mutation distribution ##########################
gene_fusion <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/7_Pathological_features/gene_fusion.xlsx")
gene_fusion<-gene_fusion[,c(1,2,6)]
gene_list<-merge(gene_fusion,Merged_data,all.x=F,all.y=F,by.x="sampleID",by.y="TestID")
gene_list<-gene_list[,c("SampleID","driver gene")]
unique(gene_list$`driver gene`)
FGFR2_mut_sample<-gene_list[which(gene_list$`driver gene`=="FGFR2"),]$SampleID
KRAS_mut_sample<-gene_list[which(gene_list$`driver gene`=="KRAS"),]$SampleID
PIK3CA_mut_sample<-gene_list[which(gene_list$`driver gene`=="PIK3CA"),]$SampleID
NTRK1_mut_sample<-gene_list[which(gene_list$`driver gene`=="NTRK1"),]$SampleID
BRAF_mut_sample<-gene_list[which(gene_list$`driver gene`=="BRAF"),]$SampleID

Merged_data$FGFR_mut<-""
Merged_data[which(Merged_data$SampleID%in%FGFR2_mut_sample),]$FGFR_mut<-"FGFR+_"
Merged_data$KRAS_mut<-""
Merged_data[which(Merged_data$SampleID%in%KRAS_mut_sample),]$KRAS_mut<-"KRAS+_"
Merged_data$PIK3CA_mut<-""
Merged_data[which(Merged_data$SampleID%in%PIK3CA_mut_sample),]$PIK3CA_mut<-"PIK3CA+_"
Merged_data$NTRK1_mut<-""
Merged_data[which(Merged_data$SampleID%in%NTRK1_mut_sample),]$NTRK1_mut<-"NTRK1+_"
Merged_data$BRAF_mut<-""
Merged_data[which(Merged_data$SampleID%in%BRAF_mut_sample),]$BRAF_mut<-"BRAF+_"
Merged_data$mut_5in1<-""
Merged_data$mut_5in1<-paste0(Merged_data$FGFR_mut,Merged_data$KRAS_mut,Merged_data$PIK3CA_mut,Merged_data$NTRK1_mut,Merged_data$BRAF_mut)

Merged_data[which(Merged_data$mut_5in1==Merged_data[1,39]),]$mut_5in1<-"wild_type"
X<-Merged_data[,c("Anno_cluster","mut_5in1")]
colnames(X)<-c("molecular_type","Mut_type")
X<-as.data.frame(table(X))
count<-as.data.frame(table(Merged_data$Anno_cluster))
for (i in unique(X$molecular_type)) {
  for (l in unique(X$Mut_type)) {
    X[which(X$molecular_type==i & X$Mut_type==l),]$Freq<-X[which(X$molecular_type==i & X$Mut_type==l),]$Freq/count[which(count$Var1==i),]$Freq
  }
}
X$Mut_type<-factor(X$Mut_type,levels = c("wild_type","BRAF+_", "FGFR+_", "KRAS+_","PIK3CA+_"))

significance_p_value<-fisher.test(table(Merged_data$Anno_cluster, Merged_data$mut_5in1),simulate.p.value = T,B=2e6)
ggplot(X,aes(x=molecular_type,y=Freq,fill=Mut_type))+
  geom_bar(stat = "identity",position="stack")+
  xlab("Molecular subtype")+ ylab("Proportion of mut type")+
  theme(axis.text.x =element_text(face="bold",color="Black", size=10, angle=25)
        ,axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))+
  ggtitle("Mut_type sample distribution in each molecular type",subtitle = "Fisher's Exact Test with simulated p-value: 0.1179 (based on 2e+06 replicates)")+scale_fill_manual("legend", values = tol21rainbow[c(1,3,5,7,9,11,13,15,17,19)])
ggsave("./7_Pathological_features//Mutation/mutation_distribution.png",width = 12,height=8)

##每亚型组的样本proportion vs anatomical
table(Merged_data$Anno_cluster,Merged_data$histological_type)

#                                  Distal intrahepatic perihilar
#mesenchymal and immunosupressive     13           10        35
#Metabolic and proliferative              40           17        49



########## Chapter 8、DGE & classifier #############################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/")
tmp<-raw_TPM[Protein_coding_genes,Merged_data$SampleID]
#write.table(tmp,file = "./raw_tpm.txt",sep="\t",quote=F,row.names = T)
#write.table(Merged_data[,c("SampleID","Anno_cluster")],file = "./class.txt",sep="\t",quote=F)
#write.table(Merged_data,file ="cluster_information.txt",sep="\t",quote=F )
data.clin <- read.table("样本分组.txt",header = T,sep = "\t")
data.count <- read.csv("~/Desktop/GenePlus/project/胆管癌/2022_12_28/genecount.csv")
colnames(data.count)[2:ncol(data.count)] <- gsub("X","",colnames(data.count)[2:ncol(data.count)])
data.count <- aggregate(data.count[,2:ncol(data.count)],by = list(data.count$X),FUN = median)
rownames(data.count) <- data.count[,1]
data.count <- data.count[,-1]
data.clin <- dplyr::arrange(data.clin,Anno_cluster)
# 批量差异分析
data.count.type <- data.count[,data.clin$SampleID] # count值
list.anno.sub <- unique(data.clin$Anno_cluster)
# 自定义
padjcutoff = 0.01 
Log2FCcutoff= 1.5
for (anno.sub in list.anno.sub) {
  anno.case <- data.clin[data.clin$Anno_cluster == anno.sub,]
  anno.contr <- data.clin[data.clin$Anno_cluster != anno.sub,]
  
  anno.case$group <- "case"
  anno.contr$group <- "control"
  anno.cc <- rbind(anno.case,anno.contr)
  
  ## DESeq2
  condition = factor(anno.cc$group)
  coldata <- data.frame(row.names = factor(anno.cc$SampleID), condition)
  data.count.type[is.na(data.count.type)] <- 0
  data.order <- anno.cc$SampleID
  data.count.type <- data.count.type[,data.order]
  
  dds <- DESeqDataSetFromMatrix(countData = round(data.count.type),
                                colData = coldata,
                                design = ~condition)
  
  dds$condition<- relevel(dds$condition, ref = "control") # 指定哪一组作为对照组
  dds <- DESeq(dds)  
  res <- results(dds, contrast = c('condition', 'case', 'control')) #提取分析结果res。注意，需将 tumor 在前，normal 在后，意为 tumor 相较于normal中哪些基因上调/下调
  res1 <- as.data.frame(results(dds))
  
  res1=subset(res1,padj!="NA") 
  res1=subset(res1,log2FoldChange!="NA") 
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  
  res1[which(res1$log2FoldChange > Log2FCcutoff & res1$padj <padjcutoff ),'sig'] <- 'up'
  res1[which(res1$log2FoldChange < (-Log2FCcutoff) & res1$padj <padjcutoff),'sig'] <- 'down'
  res1[which(abs(res1$log2FoldChange) <= Log2FCcutoff | res1$padj >= padjcutoff),'sig'] <- 'none'
  res1_updown <- subset(res1, sig %in% c('up', 'down'))
  
  path.all <- paste("all_DESeq2.",anno.sub,".txt")
  path.updown <- paste("updown_005_DESeq2.",anno.sub,".txt")
  write.table(res1, file=path.all, sep = '\t', col.names = NA, quote = FALSE)
  write.table(res1_updown, path.updown, sep = '\t', col.names = NA, quote = FALSE)
}

#write.table(res1, file="./all.deseq.txt", sep = '\t', col.names = NA, quote = FALSE)
#write.table(res1_updown, file="./deseq_up_down.txt", sep = '\t', col.names = NA, quote = FALSE)

################  volcano plot###############
Summary <- read.delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/DGE/all.deseq.txt")
Summary<-Summary[order(Summary$log2FoldChange,decreasing = T),]
Gene2Biotype <- read.delim("~/Desktop/GenePlus/project/胆管癌/genecount/Gene2Biotype.txt", header=FALSE)
Protein_coding_genes<-unique(Gene2Biotype[which(Gene2Biotype$V5=="protein_coding"),]$"V3")
Summary<-Summary[which(Summary$X%in%Protein_coding_genes),]
Summary$sig<-NA
Summary[which(abs(Summary$log2FoldChange)<1 | Summary$padj>0.01),]$sig<-"none"
Summary[which(Summary$log2FoldChange>1 & Summary$padj<0.01),]$sig<-"up"
Summary[which(Summary$log2FoldChange<(-1) & Summary$padj<0.01),]$sig<-"down"
DGE_data_tmp<-Summary

ggplot(data=DGE_data_tmp, aes(x=log2FoldChange, y=-log10(padj),colour=sig, fill=sig)) +geom_point(alpha=0.4,size=2) + 
  geom_hline(yintercept = -log10(0.01),lty=1,col="grey",lwd=1)+
  geom_vline(xintercept=c(-1,1),lty=1,col="grey",lwd=1)+
  theme_bw() +xlim(c(-5,5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme_classic(base_size = 16)+
  labs(x="Log2(Fold change)",y="-log10 (p-value)",title ="mesenchymal vs proliferative: DGE" )+ 
  scale_color_manual(values=c("down"="royalblue", "up"="red","none"="grey"))
ggsave(filename = "./DGE_volcano_plot.png")

Summary <- read.delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/all.deseq.txt")
Summary<-Summary[order(Summary$log2FoldChange,decreasing = T),]
Gene2Biotype <- read.delim("~/Desktop/GenePlus/project/胆管癌/genecount/Gene2Biotype.txt", header=FALSE)
Protein_coding_genes<-unique(Gene2Biotype[which(Gene2Biotype$V5=="protein_coding"),]$"V3")
Summary$sig<-NA
Summary[which(abs(Summary$log2FoldChange)<1 | Summary$padj>0.01),]$sig<-"none"
Summary[which(Summary$log2FoldChange>1 & Summary$padj<0.01),]$sig<-"up"
Summary[which(Summary$log2FoldChange<(-1) & Summary$padj<0.01),]$sig<-"down"

genes.df<-bitr(Summary[which(Summary$sig=="up"),]$X, fromType = as.character("SYMBOL"), toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego_up<-enrichGO(gene = genes.df$ENTREZID,OrgDb = org.Hs.eg.db, 
                 ont="BP", pvalueCutoff=0.05, 
                 pAdjustMethod = "BH",
                 qvalueCutoff =0.1,readable = TRUE)

genes.df<-bitr(Summary[which(Summary$sig=="down"),]$X, fromType = as.character("SYMBOL"), toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego_down<-enrichGO(gene = genes.df$ENTREZID,OrgDb = org.Hs.eg.db, 
                   ont="BP", pvalueCutoff=0.05, 
                   pAdjustMethod = "BH", 
                   qvalueCutoff =0.1)

pdf(file ="~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/c1vsc2_up.pdf",width = 10,height = 10)
dotplot(ego_up,showCategory=20, title ="upreg in mesenchymal and immunosupressive")
dev.off()
pdf(file ="~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/c1vsc2_down.pdf",width = 10,height = 10)
dotplot(ego_down,showCategory=20, title ="upreg in metabolic and proliferative")
dev.off()
#write.table(Summary,file = "Summary.txt",sep="\t",quote = F,row.names = F)
#################

##Classifier:用 class neighbor tools from GenePattern based on signal-to-noise distance##
SNR_gene_seletion <- read_excel("class_neighbors/486581/Protein_coding.genelist.xlsx")
SNR_gene_seletion<-as.data.frame(SNR_gene_seletion)
#SNR_gene_seletion<-SNR_gene_seletion[which(SNR_gene_seletion$Score>0.1),]
Gene_selected<-SNR_gene_seletion[c(1:20,121:130),]
SNR_gene_seletion$Classifier<-NA
SNR_gene_seletion[c(1:20),]$Classifier<-"Mesenchymal_immunosupressive"
SNR_gene_seletion[c(21:120),]$Classifier<-"Not_selected"
SNR_gene_seletion[c(121:130),]$Classifier<-"Metabolic_proliferative"
SNR_gene_seletion[c(131:240),]$Classifier<-"Not_selected"
pdf(file = "./class_neighbors/486581/SNR_score vs gene rank.pdf",width=10,height=8)
ggscatter(SNR_gene_seletion, x = "#", y = "Score",
          color ="Classifier",palette =c("red","blue","grey") ,
          alpha = 0.6, pch = 21,size = 4,
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
          ,title = "SNR score vs gene rank",xlab = "gene int",ylab = "SNR score")+
  geom_vline(xintercept=c(a<-c(1,20,121,130)),linetype = "dashed", colour = "red")
dev.off()
ggsave(filename = "SNR_score vs gene rank.png")


### complexheatmap
gen<-data.frame(gene_symbol=Gene_selected$Feature,description=c(rep("Mesenchymal_classifier",times=20),rep("Metabolic_classifier",times=10)))
#write.table(gen,file = "./NTP_genelist_genepattern.txt",sep="\t",quote=F,row.names = T)
mat<-tmp
mat<-as.data.frame(mat)
rownames(mat)<-rownames(tmp)
mat<-mat[gen$gene_symbol,]
mat<-mat[which(!is.na(rowSums(mat))),]
Merged_data<-Merged_data[order(Merged_data$Anno_cluster,decreasing=T),]
mat<-mat[gen$gene_symbol,Merged_data$SampleID]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-data.frame(molecular_cluster=Merged_data$Anno_cluster)
rownames(annotation_col)<-Merged_data$SampleID
annotation_row<-data.frame(classifier=gen$description)
rownames(annotation_row)<-gen$gene_symbol
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))

###先z-score 后 scale
library(ComplexHeatmap)
pdf(file="./class_neighbors/486581/heatmap_zscore_classifier.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = T,cluster_rows = T,
         main="heatmap_zscore_classifier (z-score)",
         annotation_col = annotation_col,
         annotation_row=annotation_row,
         row_split=annotation_row$classifier,
         show_colnames = F,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()

#### NTP_Genepattern_prediction
NTP_prediction_result <- read_delim("class_neighbors/486581/NTP/486601/NTP_prediction_result.xls",  delim = "\t", escape_double = FALSE,trim_ws = TRUE)
NTP_prediction_result<-as.data.frame(NTP_prediction_result)
NTP_prediction_result$predict_molecular<-"Metabolic and priliferative"
NTP_prediction_result[which(NTP_prediction_result$predict.label==1),]$predict_molecular<-"mesenchymal and immunosupressive"
#Merged_data[which(Merged_data$Anno_cluster=="Metabolic and proliferative"),]$Anno_cluster<-"Metabolic and priliferative"
rownames(NTP_prediction_result)<-NTP_prediction_result$sample.names
NTP_prediction_result<-NTP_prediction_result[Merged_data$SampleID,]
Merged_data<-cbind(Merged_data,NTP_prediction_result)

table(Merged_data$BH.FDR<0.05)
#FALSE  TRUE 
#38   126 
##high confidence (false discovery rate < 0.05)in 76.8% of our cohort (126/164)

table(Merged_data[which(Merged_data$BH.FDR<0.05),]$Anno_cluster,Merged_data[which(Merged_data$BH.FDR<0.05),]$predict_molecular)
#precision: Precision refers to positive predictive value 
#(number of true positives/number of positive calls)
## precision: 121/126 (96.03%)

#total
#                   mesenchymal and immunosupressive   Metabolic and priliferative
#mesenchymal and immunosupressive     54                         4
#Metabolic and priliferative          19                        87
#Total                                73                        91

#positive calls
#                   mesenchymal and immunosupressive   Metabolic and priliferative
#mesenchymal and immunosupressive     50                         2
#Metabolic and priliferative          3                         71
#Total                                53                        73

## accuracy per class (true positives + true negatives / samples) 
## mesenchymal and immunosupressive= (50+4-2)/(54+4)=0.89655
## Metabolic and priliferative = (71+19-3)/(87+19)=0.82075

Merged_data$Corrected_molecular_prection<-NA
Merged_data$Corrected_molecular_prection<-"unclassified (FDR>0.05)"
Merged_data[which(Merged_data$BH.FDR<0.05 & Merged_data$predict_molecular=="Metabolic and priliferative"),]$Corrected_molecular_prection<-"Metabolic and priliferative"
Merged_data[which(Merged_data$BH.FDR<0.05 & Merged_data$predict_molecular=="mesenchymal and immunosupressive"),]$Corrected_molecular_prection<-"Mesenchymal and immunosupressive"
table(Merged_data$Corrected_molecular_prection)

annotation_col<-HeatmapAnnotation(
  FDR=anno_empty(border = T,height = unit(3, "cm")),
  molecular_cluster=Merged_data$Anno_cluster,
  predicted=Merged_data$Corrected_molecular_prection,
  #FDR=anno_points(Merged_data$BH.FDR,border = T,height = unit(3, "cm")),
  col =list(molecular_cluster = c("mesenchymal and immunosupressive" = "red", "Metabolic and priliferative" ="blue"),
            predicted= c("unclassified (FDR>0.05)" = "grey", "Mesenchymal and immunosupressive" ="red","Metabolic and priliferative" ="blue")))

annotation_row<-rowAnnotation(classifier=gen$description,col =list(
  classifier = c("Mesenchymal_classifier" = "red", "Metabolic_classifier" ="blue")))

ht<-Heatmap(as.matrix(mat),row_km = 2,
            name="heatmap_zscore_classifier (z-score)",
            cluster_rows = T,
            cluster_columns = F,
            show_column_names = F,
            right_annotation = annotation_row,
            top_annotation = annotation_col,
            column_order =Merged_data[order(Merged_data$Anno_cluster,Merged_data$Corrected_molecular_prection),]$SampleID)


pdf(file="./class_neighbors/486581/heatmap_final.pdf",width = 12,height = 10)
ht
co = column_order(ht)
value=Merged_data$BH.FDR
decorate_annotation("FDR", {
  # value on x-axis is always 1:ncol(mat)
  x = 1:ncol(mat)
  # while values on y-axis is the value after column reordering
  value = value[co]
  pushViewport(viewport(xscale = c(0.5, 164.5),yscale = c(0, 1)))
  grid.lines(c(0.5, 163.5),c(0.05, 0.05), gp = gpar(lty = 10),
             default.units = "native")
  grid.points(x, value, pch = 3, size = unit(1.8, "mm"),
              gp = gpar(col = ifelse(value > 0.05, "grey", "orange")), default.units = "native")
  grid.yaxis(at = c(0, 0.05, 1))
  popViewport()
})
dev.off()

#######外部数据验证 classifier 
######################################### TCGA-CHOL #################################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/TCGA_CHOL/")
TCGA_CHOL_tmp <- read_csv("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/TCGA_CHOL/TCGA.CHOL_tmp_36.csv")

library(ConsensusClusterPlus)
NMF_mat<-TCGA_CHOL_tmp
NMF_mat<-NMF_mat[which(rowSums(NMF_mat[,-1])>0),]
Gene_tmp<-NMF_mat$...1
NMF_mat<-NMF_mat[,-1]
rownames(NMF_mat)<-Gene_tmp
mads<-apply(NMF_mat, 1, mad)
NMF_mat<-NMF_mat[rev(order(mads)),]
NMF_mat= sweep(NMF_mat,1, apply(NMF_mat,1,median,na.rm=T))
library(ConsensusClusterPlus)

############ To linux 
#title="~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/TCGA_CHOL/"
#results = ConsensusClusterPlus(as.matrix(NMF_mat),maxK=2,reps=800,pItem=0.8,pFeature=1, title=title,clusterAlg="km",distance="euclidean",plot="png")
#paper:Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data
#In this case, the CDF will be horizontal on the interval [0, 1) as in panel B. In a less-than-ideal case, the CDF will rise more or less gradually to the right, as in panel A. We also want to identify as many clusters as possible without sacrificing the goal that most consensus indices should be near 0 or 1.

load("/Users/zouxinchen/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/TCGA_CHOL/Consensus_clustering.RData")
icl = calcICL(results,title=title,plot="png")
icl[["clusterConsensus"]]
k=2
consensus_clustering<-icl[["itemConsensus"]]
consensus_clustering<-consensus_clustering[which(consensus_clustering$k==2),]
consensus_clustering<-consensus_clustering[which(consensus_clustering$itemConsensus>0.5),]
cluster_information<-data.frame(sample_name=consensus_clustering$item,clustering=consensus_clustering$cluster)
cluster_information<-unique(cluster_information)
TCGA_CHOL_NTP_prediction_result <- read_delim("TCGA_CHOL_NTP_prediction_result.xls", 
                                              delim = "\t", escape_double = FALSE,trim_ws = TRUE)

Merged_TCGA<-merge(TCGA_CHOL_NTP_prediction_result,cluster_information,by.x = "sample_names",by.y="sample_name",all.x=T,all.y=T)
TCGA_survival <- read_delim("data_clinical_patient.txt",  delim = "\t", escape_double = FALSE,trim_ws = TRUE)
Merged_TCGA<-cbind(Merged_TCGA,TCGA_survival)
mat<-TCGA_CHOL_tmp
mat<-mat[which(rowSums(mat[,-1])>0),]
Gene_tmp<-mat$...1
mat<-mat[,-1]
rownames(mat)<-Gene_tmp
mat<-as.data.frame(mat)
rownames(mat)<-Gene_tmp
mat<-mat[gen$gene_symbol,]
mat<-mat[which(!is.na(rowSums(mat))),]
Merged_TCGA<-Merged_TCGA[order(Merged_TCGA$clustering,Merged_TCGA$predict.label,decreasing=T),]
mat<-mat[gen$gene_symbol,Merged_TCGA$sample_names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-HeatmapAnnotation(
  FDR=anno_empty(border = T,height = unit(3, "cm")),
  molecular_cluster=Merged_TCGA$clustering,
  predicted=Merged_TCGA$predict.label,
  col =list(molecular_cluster = c("1" = "red", "2" ="blue"),
            predicted= c("1"="red","2"="blue")))
annotation_row<-rowAnnotation(classifier=gen$description,col =list(
  classifier = c("Mesenchymal_classifier" = "red", "Metabolic_classifier" ="blue")))
ht<-Heatmap(as.matrix(mat),
            name="heatmap_zscore_classifier (z-score)",
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = F,
            right_annotation = annotation_row,
            top_annotation = annotation_col)

pdf(file="./TCGA_NTP_heatmap.pdf",width = 12,height = 10)
ht
co = column_order(ht)
value=Merged_TCGA$BH.FDR
decorate_annotation("FDR", {
  # value on x-axis is always 1:ncol(mat)
  x = 1:ncol(mat)
  # while values on y-axis is the value after column reordering
  value = value[co]
  pushViewport(viewport(xscale = c(0.5, 36.5),yscale = c(0, 0.1)))
  grid.lines(c(0.5, 36.5),c(0.05, 0.05), gp = gpar(lty = 10),
             default.units = "native")
  grid.points(x, value, pch = 3, size = unit(1.8, "mm"),
              gp = gpar(col = ifelse(value > 0.05, "grey", "orange")), default.units = "native")
  grid.yaxis(at = c(0, 0.05, 1))
  popViewport()
})
dev.off()
#write.table(Merged_TCGA,file = "Merged_TCGA.txt",sep="\t",quote=F,row.names = T)
#TCGA DGE分析

#################### 不放结果 survival outcome #########################
mydata = Merged_TCGA
mydata[which(mydata$DFS_STATUS=="1:Recurred/Progressed"),]$DFS_STATUS<-1
mydata[which(mydata$DFS_STATUS=="0:DiseaseFree"),]$DFS_STATUS<-0
mydata[which(mydata$OS_STATUS=="1:DECEASED"),]$OS_STATUS<-1
mydata[which(mydata$OS_STATUS=="0:LIVING"),]$OS_STATUS<-0
mydata<-as.data.frame(mydata)
groups=c("predict.label","clustering")
for(vals in groups){
  dfDFS = mydata[which(!is.na(mydata[,vals]) & 
                         !is.na(mydata["DFS_MONTHS"]) &
                         mydata["DFS_STATUS"] !="Unknown" &
                         !is.na(mydata["DFS_STATUS"])),]
  dfDFS$DFS_STATUS = as.numeric(as.character(dfDFS$DFS_STATUS))
  dfDFS$DFS_MONTHS = as.numeric(as.character(dfDFS$DFS_MONTHS))
  dfDFS["group"] = dfDFS[,vals]
  fitDFS <- survfit(Surv(DFS_MONTHS,DFS_STATUS) ~ group,data =dfDFS)
  # coxDFS=coxph(Surv(DFS_MONTHS, DFS_Event) ~group+Stage,data = dfDFS)
  #DFS
  ggsurvplot(fitDFS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="DFS",
             xlab="Time(Month)",
             pval.size= 8,pval.coord=c(6,0.3),
             size=1.2
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1),
                    legend.text = c(14,"bold", "darkblue")
    )
  name1 = paste(vals,"DFS.png",sep="_")
  name2 = paste(vals,"DFS.pdf",sep="_")
  ggsave(name1,width=12,height = 8)
  ggsave(name2,width=12,height = 8)
  #OS
  dfOS = mydata[which(!is.na(mydata[,vals]) & 
                        !is.na(mydata["OS_MONTHS"]) &
                        !is.na(mydata["OS_STATUS"])&
                        (mydata["OS_STATUS"] !="Unknown")),]
  dfOS$OS_STATUS = as.numeric(as.character(dfOS$OS_STATUS))
  dfOS$OS_MONTHS = as.numeric(as.character(dfOS$OS_MONTHS))
  dfOS["group"] = dfOS[,vals]
  fitOS <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ group,data =dfOS) 
  ggsurvplot(fitOS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="OS",
             xlab="Time(Month)",
             pval.size= 8,pval.coord=c(6,0.3)
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1)
    )
  name1 = paste(vals,"OS.png",sep="_")
  name2 = paste(vals,"OS.pdf",sep="_")
  ggsave(name1,width=12,height = 8)
  ggsave(name2,width=12,height = 8)
}


############################# molecular traits ##################################
all_gene_sets<-msigdbr(species = "Homo sapiens")
all_gene_sets<-as.data.frame(all_gene_sets)
Hallmarker_gene_set<-all_gene_sets[which(all_gene_sets$gs_cat=="H"),]
geneset_description<-unique(Hallmarker_gene_set$gs_description)
mat<-TCGA_CHOL_tmp
mat<-mat[which(rowSums(mat[,-1])>0),]
Gene_tmp<-mat$...1
mat<-mat[,-1]
rownames(mat)<-Gene_tmp
mat<-as.data.frame(mat)
rownames(mat)<-Gene_tmp
gene_list<-Gene_tmp
gene_list<-gene_list[gene_list%in%Hallmarker_gene_set$gene_symbol]
sample_hallmarker_expre<-mat[gene_list,]
Hallmarker_gene_set<-Hallmarker_gene_set[,c(4,15)]
Hallmarker_gene_set<-as.data.frame(Hallmarker_gene_set)
mat<-sample_hallmarker_expre
gen<-Hallmarker_gene_set
head(gen)
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="ssgsea", verbose=FALSE, parallel.sz=2)

#standarization by z-score
es<-t(scale(t(es), scale=TRUE, center=TRUE))
geneset_description <- read.csv("~/Desktop/GenePlus/project/胆管癌/2022_12_28/4_Hallmarker_GSEA//geneset_description.csv", row.names=1)
geneset_description<-geneset_description[rownames(es),]
rownames(es)<-geneset_description$function.
final_enriched_scores<-as.data.frame(matrix(NA,nrow=50))

for (i in 1:2) {
  tmp_sample<-Merged_TCGA[which(Merged_TCGA$predict.label==i),]$sample_names
  Mean_cluster_score<-rowMeans(as.matrix((es[,tmp_sample])))
  final_enriched_scores[,i]<-Mean_cluster_score
}
rownames(final_enriched_scores)<-rownames(es)
colnames(final_enriched_scores)<-c("cluster1","cluster2")
colnames(final_enriched_scores)<-c("mesenchymal and immunosupressive","Metabolic and proliferative")
bk<-c(seq(min(final_enriched_scores),0,length.out = 50),seq(0.00001,max(final_enriched_scores),length.out = 50))
###先z-score 后 scale
pdf(file="./Hallmark_GSEA_TCGA_validation.pdf",width = 8,height = 10)
pheatmap::pheatmap(as.matrix(final_enriched_scores),cluster_cols = F,cluster_rows = T,color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk,main="pheatmap of hallmarker enrichment_TCGA_enrichment",scale ="row")
dev.off()

#################### Ferroptosis resistance ###############
Ferroptosis_related <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/11_additional_exploration/ferroptosis resistance/Ferroptosis_related.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description="Ferroptosis-related")
{
  mat<-TCGA_CHOL_tmp
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$...1
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(rownames(mat)%in%gen$gene_symbol),]
  mat<-mat[which(!is.na(rowSums(mat))),]}
mat
gen
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="zscore",verbose=FALSE, parallel.sz=2)
es<-t(es)
es<-as.data.frame(es)
es<-es[Merged_TCGA$sample_names,]
explor_info<-cbind(Merged_TCGA,es)
my_comparisons = list(c("1","2"))
label.y = c(1)

ggplot(explor_info,aes(x=factor(explor_info$predict.label,levels = c("1","2"))
                       ,y=explor_info$es,
                       color=predict.label))+
  geom_boxplot(notch = F)+
  labs(x = "molecular subtype",
       y = paste0("TCGA_validation:log2 enrichment score: Ferroptosis"),
       title =paste0("enrichment of Ferroptosis in each molecular subtype"))+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
  stat_compare_means(label.y = 4)+
  stat_compare_means(comparisons = my_comparisons, label.y = 5, method ="t.test")+
  theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
ggsave(filename ="Ferroptosis_z_score_gsva.png",width = 10,height=8)

Merged_TCGA<-Merged_TCGA[order(Merged_TCGA$predict.label,decreasing=F),]
{mat<-TCGA_CHOL_tmp
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$...1
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(rownames(mat)%in%gen$gene_symbol),]
  mat<-mat[which(!is.na(rowSums(mat))),]}
mat<-mat[,Merged_TCGA$sample_names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
Merged_TCGA$predict.label<-factor(Merged_TCGA$predict.label,levels = c(1,2))
annotation_col<-data.frame(molecular_cluster=as.character(Merged_TCGA$predict.label))
rownames(annotation_col)<-Merged_TCGA$sample_names
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))
###先z-score 后 scale
pdf(file="./zscore_expression_Ferroptosis.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = T,cluster_rows = T,
         main="Ferroptosis-related gene expression (z-score)",
         annotation_col = annotation_col,
         show_colnames = F,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()

################### MHC molecules ###################
Ferroptosis_related <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/11_additional_exploration/MHC_molecules/MHC molecules.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description=Ferroptosis_related$Category)

{
  mat<-TCGA_CHOL_tmp
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$...1
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(rownames(mat)%in%gen$gene_symbol),]
  mat<-mat[which(!is.na(rowSums(mat))),]
}
mat<-mat[gen$gene_symbol,Merged_TCGA$sample_names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-data.frame(molecular_cluster=Merged_TCGA$predict.label)
rownames(annotation_col)<-Merged_TCGA$sample_names
annotation_row<-data.frame(Category=gen$gs_description)
rownames(annotation_row)<-gen$gene_symbol
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))
###先z-score 后 scale
library(ComplexHeatmap)
pdf(file="./zscore_expression_molecular_MHCs_antigens.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = T,cluster_rows = T,
         main="MHC_moleculars/co-inhibitor/Co-simulator expression (z-score)",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_colnames = F,
         row_split=annotation_row$Category,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()






######################################################################### cohort 3 ############################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/Cohort3/")
cohort_3_tpm <- read_xlsx("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/Cohort3/Cohort3_TPM.xlsx")
cohort3_protein <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/Cohort3/protein_expression.xlsx")

{
  mat<-cohort_3_tpm
  mat[,-1]<- (2^cohort_3_tpm[,-1])-1
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(!is.na(rowSums(mat))),]
  colnames(mat)<-paste0("s",colnames(mat))
}
#write.table(mat,file="cohort3_TPM.txt",sep="\t",quote=F,row.names = T)

library(ConsensusClusterPlus)
Cohort3_NTP_prediction_result <- read_delim("489341/Cohort3_prediction_result.xls", 
                                            delim = "\t", escape_double = FALSE,trim_ws = TRUE)
Merged_Cohort3<-as.data.frame(Cohort3_NTP_prediction_result)
Merged_Cohort3[which(Merged_Cohort3$BH.FDR>0.05),]$predict.label<-"unclassified"
table(Merged_Cohort3$predict.label)
clinical_information_cohort3 <- read_excel("clinical_information.xlsx")
Merged_Cohort3<-merge(Merged_Cohort3,clinical_information_cohort3,by.x="sample.names",by.y="Patient_ID",all.x=T,all.y=F)


{
  mat<-cohort_3_tpm
  mat[,-1]<- (2^cohort_3_tpm[,-1])-1
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[gen$gene_symbol,]
  mat<-mat[which(!is.na(rowSums(mat))),]
  colnames(mat)<-paste0("s",colnames(mat))
}
gen<-data.frame(gene_symbol=Gene_selected$Feature,description=c(rep("Mesenchymal_classifier",times=20),rep("Metabolic_classifier",times=10)))
Merged_Cohort3<-Merged_Cohort3[order(Merged_Cohort3$predict.label),]
Merged_Cohort3<-Merged_Cohort3[which(Merged_Cohort3$predict.label!="unclassified"),]
mat<-mat[gen$gene_symbol,Merged_Cohort3$sample.names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-HeatmapAnnotation(
  predicted=Merged_Cohort3$predict.label,
  col =list(predicted= c("1"="red","2"="blue","unclassified"="grey")))
annotation_row<-rowAnnotation(classifier=gen$description,col =list(
  classifier = c("Mesenchymal_classifier" = "red", "Metabolic_classifier" ="blue")))

ht1<-Heatmap(as.matrix(mat),
             name="heatmap_zscore_classifier (z-score)",
             cluster_rows = T,
             cluster_columns = F,
             show_column_names = F,
             right_annotation = annotation_row,
             top_annotation = annotation_col,
             row_split = gen$description,
             column_split = Merged_Cohort3$predict.label)
pdf(file="./cohort3_NTP_heatmap.pdf",width = 12,height = 10)
ht1
dev.off()

####protein expression
{
  mat<-cohort3_protein
  mat<-cohort3_protein[which(cohort3_protein$Protein_NAME%in%gen$gene_symbol),]
  Gene_tmp<-mat$Protein_NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat <- as.data.frame(sapply(mat, as.numeric))
  rownames(mat)<-Gene_tmp
  #mat<-mat[which(!is.na(rowSums(mat))),]
}
gen<-data.frame(gene_symbol=Gene_selected$Feature,description=c(rep("Mesenchymal_classifier",times=20),rep("Metabolic_classifier",times=10)))
gen<-gen[which(gen$gene_symbol%in%Gene_tmp),]
#Merged_Cohort3<-Merged_Cohort3[which(Merged_Cohort3$sample.names%in%colnames(mat)&Merged_Cohort3$predict.label!="unclassified"),]
Merged_Cohort3<-Merged_Cohort3[order(Merged_Cohort3$predict.label),]
additional_sample<-Merged_Cohort3$sample.names[Merged_Cohort3$sample.names%notin%colnames(mat)]
additional_sample_mat<-matrix(NA,nrow = 20,ncol =length(additional_sample))
colnames(additional_sample_mat)<-additional_sample
mat<-cbind(mat,additional_sample_mat)
mat<-mat[gen$gene_symbol,Merged_Cohort3$sample.names]

#mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-HeatmapAnnotation(
  predicted=Merged_Cohort3$predict.label,
  col =list(predicted= c("1"="red","2"="blue","unclassified"="grey")))
annotation_row<-rowAnnotation(classifier=gen$description,col =list(
  classifier = c("Mesenchymal_classifier" = "red", "Metabolic_classifier" ="blue")))
ht2<-Heatmap(as.matrix(mat),
             name="protein expression of classifier genes (log2 trans median normalized)",
             cluster_rows = T,
             cluster_columns = F,
             show_column_names = F,
             right_annotation = annotation_row,
             row_split = gen$description,
             column_split = Merged_Cohort3$predict.label)

pdf(file="./cohort3_protein_exp_classifier_gene_heatmap.pdf",width = 12,height = 10)
ht2
dev.off()

pdf(file="./complexheatmap.pdf",width = 14,height = 10)
ht1 %v% ht2
dev.off()


################################## 不放结果 survival outcome ####
mydata = Merged_Cohort3[which(Merged_Cohort3$predict.label!="unclassified"),]
mydata<-as.data.frame(mydata)
groups=c("predict.label")
for(vals in groups){
  dfOS = mydata[which(!is.na(mydata[,vals]) & 
                        !is.na(mydata["OS_DAYS"]) &
                        !is.na(mydata["OS_STATUS"])&
                        (mydata["OS_STATUS"] !="Unknown")),]
  dfOS$OS_STATUS = as.numeric(as.character(dfOS$OS_STATUS))
  dfOS$OS_DAYS = as.numeric(as.character(dfOS$OS_DAYS))
  dfOS["group"] = dfOS[,vals]
  fitOS <- survfit(Surv(OS_DAYS,OS_STATUS) ~ group,data =dfOS) 
  ggsurvplot(fitOS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="OS",
             xlab="Time(Month)",
             pval.size= 8,pval.coord=c(6,0.3)
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1)
    )
  name1 = paste(vals,"OS.png",sep="_")
  name2 = paste(vals,"OS.pdf",sep="_")
  ggsave(name1,width=12,height = 8)
  ggsave(name2,width=12,height = 8)
}


############################# molecular traits ##################################
all_gene_sets<-msigdbr(species = "Homo sapiens")
all_gene_sets<-as.data.frame(all_gene_sets)
Hallmarker_gene_set<-all_gene_sets[which(all_gene_sets$gs_cat=="H"),]
geneset_description<-unique(Hallmarker_gene_set$gs_description)
{
  mat<-cohort_3_tpm
  mat[,-1]<- (2^cohort_3_tpm[,-1])-1
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(!is.na(rowSums(mat))),]
  colnames(mat)<-paste0("s",colnames(mat))
}
gene_list<-Gene_tmp
gene_list<-gene_list[gene_list%in%Hallmarker_gene_set$gene_symbol]
sample_hallmarker_expre<-mat[gene_list,]
Hallmarker_gene_set<-Hallmarker_gene_set[,c(4,15)]
Hallmarker_gene_set<-as.data.frame(Hallmarker_gene_set)
mat<-sample_hallmarker_expre
gen<-Hallmarker_gene_set
head(gen)
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="ssgsea", verbose=FALSE, parallel.sz=2)
#standarization by z-score
es<-t(scale(t(es), scale=TRUE, center=TRUE))
geneset_description <- read.csv("~/Desktop/GenePlus/project/胆管癌/2022_12_28/4_Hallmarker_GSEA//geneset_description.csv", row.names=1)
geneset_description<-geneset_description[rownames(es),]
rownames(es)<-geneset_description$function.
final_enriched_scores<-as.data.frame(matrix(NA,nrow=50))

for (i in 1:2) {
  tmp_sample<-Merged_Cohort3[which(Merged_Cohort3$predict.label==i),]$sample.names
  Mean_cluster_score<-rowMeans(as.matrix((es[,tmp_sample])))
  final_enriched_scores[,i]<-Mean_cluster_score
}
rownames(final_enriched_scores)<-rownames(es)
colnames(final_enriched_scores)<-c("cluster1","cluster2")
colnames(final_enriched_scores)<-c("mesenchymal and immunosupressive","Metabolic and proliferative")
bk<-c(seq(min(final_enriched_scores),0,length.out = 50),seq(0.00001,max(final_enriched_scores),length.out = 50))
###先z-score 后 scale
pdf(file="./Hallmark_GSEA_cohort3_validation.pdf",width = 8,height = 10)
pheatmap::pheatmap(as.matrix(final_enriched_scores),cluster_cols = F,cluster_rows = T,color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk,main="pheatmap of hallmarker enrichment_cohort3_enrichment",scale ="row")
dev.off()

#################### Ferroptosis resistance ###############
Ferroptosis_related <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/11_additional_exploration/ferroptosis resistance/Ferroptosis_related.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description="Ferroptosis-related")
{
  mat<-cohort_3_tpm
  mat[,-1]<- (2^cohort_3_tpm[,-1])-1
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(!is.na(rowSums(mat))),]
  colnames(mat)<-paste0("s",colnames(mat))
}
mat<-mat[gen$gene_symbol,Merged_Cohort3$sample.names]
mat
gen
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="zscore",verbose=FALSE, parallel.sz=2)
es<-t(es)
es<-as.data.frame(es)
es<-es[Merged_Cohort3$sample.names,]
explor_info<-cbind(Merged_Cohort3,es)
my_comparisons = list(c("1","2"))
label.y = c(1)
explor_info<-explor_info[which(explor_info$predict.label!="unclassified"),]
ggplot(explor_info,aes(x=factor(explor_info$predict.label,levels = c("1","2"))
                       ,y=explor_info$es,
                       color=predict.label))+
  geom_boxplot(notch = F)+
  labs(x = "molecular subtype",
       y = paste0("Cohort3_validation:log2 enrichment score: Ferroptosis"),
       title =paste0("enrichment of Ferroptosis in each molecular subtype"))+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
  stat_compare_means(label.y = 4)+
  stat_compare_means(comparisons = my_comparisons, label.y = 5, method ="t.test")+
  theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
ggsave(filename ="Ferroptosis_z_score_gsva.png",width = 10,height=8)

Merged_Cohort3<-Merged_Cohort3[order(Merged_Cohort3$predict.label,decreasing=F),]
{
  mat<-cohort_3_tpm
  mat[,-1]<- (2^cohort_3_tpm[,-1])-1
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(!is.na(rowSums(mat))),]
  colnames(mat)<-paste0("s",colnames(mat))
}
Merged_Cohort3<-Merged_Cohort3[which(Merged_Cohort3$predict.label!="unclassified"),]
Merged_Cohort3$predict.label<-factor(Merged_Cohort3$predict.label,levels = c("1","2"))
mat<-mat[gen$gene_symbol,Merged_Cohort3$sample.names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-data.frame(molecular_cluster=as.character(Merged_Cohort3$predict.label))
rownames(annotation_col)<-Merged_Cohort3$sample.names
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))
###先z-score 后 scale
pdf(file="./Cohort3_zscore_expression_Ferroptosis.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = T,cluster_rows = T,
         main="Ferroptosis-related gene expression (z-score)",
         annotation_col = annotation_col,
         show_colnames = F,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()

################### MHC molecules ###################
Ferroptosis_related <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/11_additional_exploration/MHC_molecules/MHC molecules.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description=Ferroptosis_related$Category)
{
  mat<-cohort_3_tpm
  mat[,-1]<- (2^cohort_3_tpm[,-1])-1
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(!is.na(rowSums(mat))),]
  colnames(mat)<-paste0("s",colnames(mat))
}

mat<-mat[gen$gene_symbol,Merged_Cohort3$sample.names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-HeatmapAnnotation(
  predicted=Merged_Cohort3$predict.label,
  col =list(predicted= c("1"="red","2"="blue","unclassified"="grey")))

annotation_row<-rowAnnotation(classifier=gen$gs_description,col =list(
  classifier = c("MHC-I" = "red", "MHC-II" ="blue","Co-simulator"="black","Co-inhibitor"="green")))

bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))
###先z-score 后 scale
library(ComplexHeatmap)
ht3<-Heatmap(as.matrix(mat),
             name="MHC_moleculars/co-inhibitor/Co-simulator expression (z-score)",
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             right_annotation = annotation_row,top_annotation = annotation_col,
             row_split = gen$gs_description,
             column_split =Merged_Cohort3$predict.label)
pdf(file="./zscore_expression_molecular_MHCs_antigens.pdf",width = 12,height = 10)
ht3
dev.off()

####protein expression
{
  mat<-cohort3_protein
  mat<-cohort3_protein[which(cohort3_protein$Protein_NAME%in%gen$gene_symbol),]
  Gene_tmp<-mat$Protein_NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat <- as.data.frame(sapply(mat, as.numeric))
  rownames(mat)<-Gene_tmp
  #mat<-mat[which(!is.na(rowSums(mat))),]
}
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description=Ferroptosis_related$Category)
gen<-gen[which(gen$gene_symbol%in%Gene_tmp),]

Merged_Cohort3<-Merged_Cohort3[order(Merged_Cohort3$predict.label),]
additional_sample<-Merged_Cohort3$sample.names[Merged_Cohort3$sample.names%notin%colnames(mat)]
additional_sample_mat<-matrix(NA,nrow = 26,ncol =length(additional_sample))
colnames(additional_sample_mat)<-additional_sample
mat<-cbind(mat,additional_sample_mat)
mat<-mat[gen$gene_symbol,Merged_Cohort3$sample.names]

#mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-HeatmapAnnotation(
  predicted=Merged_Cohort3$predict.label,
  col =list(predicted= c("1"="red","2"="blue","unclassified"="grey")))
annotation_row<-rowAnnotation(classifier=gen$gs_description,col =list(
  classifier = c("MHC-I" = "red", "MHC-II" ="blue","Co-simulator"="black","Co-inhibitor"="green")))

ht4<-Heatmap(as.matrix(mat),
             name="protein expression of MHC_moleculars/co-inhibitor/Co-simulator",
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             right_annotation = annotation_row,
             row_split = gen$gs_description,
             column_split = Merged_Cohort3$predict.label)

pdf(file="./cohort3_protein_molecular_MHCs_antigens.pdf",width = 12,height = 10)
ht4
dev.off()

pdf(file="./complexheatmap_MHCs_antigens.pdf",width = 14,height = 14)
ht3 %v% ht4
dev.off()


############################################### our_cohort_sample_left ###########################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/class_neighbors/486581/NTP/our_cohort_validation/")
our_cohort_sample_left<-RNA_Sample[RNA_Sample%notin%Sample_analysising]
validation_cohort_prediction_result <- 
  read_delim("489550/validation_cohort_prediction_result.xls",delim = "\t", escape_double = FALSE,trim_ws = TRUE)
validation_cohort_prediction_result<-validation_cohort_prediction_result[which(validation_cohort_prediction_result$sample.names%in%our_cohort_sample_left),]
mat<-raw_TPM
validation_cohort_prediction_result
Merged_validation_cohort<-as.data.frame(validation_cohort_prediction_result)
Merged_validation_cohort[which(Merged_validation_cohort$BH.FDR>0.05),]$predict.label<-"unclassified"
Merged_validation_cohort<-merge(Merged_validation_cohort,Merged_data[which(Merged_data$SampleID%in%our_cohort_sample_left),],by.x="sample.names",by.y="SampleID",all.x=T,all.y=F)
table(Merged_validation_cohort$predict.label)

gen<-data.frame(gene_symbol=Gene_selected$Feature,description=c(rep("Mesenchymal_classifier",times=20),rep("Metabolic_classifier",times=10)))
Merged_validation_cohort<-Merged_validation_cohort[order(Merged_validation_cohort$predict.label),]
mat<-mat[gen$gene_symbol,Merged_validation_cohort$sample.names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-HeatmapAnnotation(
  predicted=Merged_validation_cohort$predict.label,
  Normal_liver_percentage=Merged_validation_cohort$Normal_liver_percentage,
  Normal_Neuron_percentage=Merged_validation_cohort$Normal_Neuron_percentage,
  Normal_Lymphoid_percentage=Merged_validation_cohort$Normal_Lymphoid_percentage,
  Normal_pancreas_percentage=Merged_validation_cohort$Normal_pancreas_percentage,
  Normal_duodenum_percentage=Merged_validation_cohort$Normal_duodenum_percentage,
  col =list(predicted= c("1"="red","2"="blue","unclassified"="grey")))

gen<-data.frame(gene_symbol=Gene_selected$Feature,description=c(rep("Mesenchymal_classifier",times=20),rep("Metabolic_classifier",times=10)))
annotation_row<-rowAnnotation(classifier=gen$description,col =list(
  classifier = c("Mesenchymal_classifier" = "red", "Metabolic_classifier" ="blue")))

ht<-Heatmap(as.matrix(mat),
            name="heatmap_zscore_classifier (z-score)",
            cluster_rows = T,
            cluster_columns = F,
            show_column_names = F,right_annotation = annotation_row,
            top_annotation = annotation_col)
pdf(file="./validation_cohort_NTP_heatmap.pdf",width = 12,height = 10)
ht
dev.off()

############################# molecular traits ##################################
all_gene_sets<-msigdbr(species = "Homo sapiens")
all_gene_sets<-as.data.frame(all_gene_sets)
Hallmarker_gene_set<-all_gene_sets[which(all_gene_sets$gs_cat=="H"),]
geneset_description<-unique(Hallmarker_gene_set$gs_description)
mat<-raw_TPM
gene_list<-rownames(raw_TPM)
gene_list<-gene_list[gene_list%in%Hallmarker_gene_set$gene_symbol]
sample_hallmarker_expre<-mat[gene_list,]
Hallmarker_gene_set<-Hallmarker_gene_set[,c(4,15)]
Hallmarker_gene_set<-as.data.frame(Hallmarker_gene_set)
mat<-sample_hallmarker_expre
gen<-Hallmarker_gene_set
head(gen)
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="ssgsea", verbose=FALSE, parallel.sz=2)
#standarization by z-score
es<-t(scale(t(es), scale=TRUE, center=TRUE))
geneset_description <- read.csv("~/Desktop/GenePlus/project/胆管癌/2022_12_28/4_Hallmarker_GSEA//geneset_description.csv", row.names=1)
geneset_description<-geneset_description[rownames(es),]
rownames(es)<-geneset_description$function.
final_enriched_scores<-as.data.frame(matrix(NA,nrow=50))

for (i in 1:2) {
  tmp_sample<-Merged_validation_cohort[which(Merged_validation_cohort$predict.label==i),]$sample.names
  Mean_cluster_score<-rowMeans(as.matrix((es[,tmp_sample])))
  final_enriched_scores[,i]<-Mean_cluster_score
}
rownames(final_enriched_scores)<-rownames(es)
colnames(final_enriched_scores)<-c("cluster1","cluster2")
colnames(final_enriched_scores)<-c("mesenchymal and immunosupressive","Metabolic and proliferative")
bk<-c(seq(min(final_enriched_scores),0,length.out = 50),seq(0.00001,max(final_enriched_scores),length.out = 50))
###先z-score 后 scale
pdf(file="./Hallmark_GSEA_Merged_validation_cohort.pdf",width = 8,height = 10)
pheatmap::pheatmap(as.matrix(final_enriched_scores),cluster_cols = F,cluster_rows = T,color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk,main="pheatmap of hallmarker enrichment_Merged_validation_cohort",scale ="row")
dev.off()

#################### Ferroptosis resistance ###############
Ferroptosis_related <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/11_additional_exploration/ferroptosis resistance/Ferroptosis_related.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description="Ferroptosis-related")
mat<-raw_TPM
mat<-mat[gen$gene_symbol,Merged_validation_cohort$sample.names]
mat
gen
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="zscore",verbose=FALSE, parallel.sz=2)
es<-t(es)
es<-as.data.frame(es)
es<-es[Merged_validation_cohort$sample.names,]
explor_info<-cbind(Merged_validation_cohort,es)
my_comparisons = list(c("1","2"))
label.y = c(1)
explor_info<-explor_info[which(explor_info$predict.label!="unclassified"),]
ggplot(explor_info,aes(x=factor(explor_info$predict.label,levels = c("1","2"))
                       ,y=log2(explor_info$es),
                       color=predict.label))+
  geom_boxplot(notch = F)+
  labs(x = "molecular subtype",
       y = paste0("Merged_validation_cohort:log2 enrichment score: Ferroptosis"),
       title =paste0("enrichment of Ferroptosis in each molecular subtype"))+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
  stat_compare_means(label.y = 4)+
  stat_compare_means(comparisons = my_comparisons, label.y = 5, method ="t.test")+
  theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
ggsave(filename ="Ferroptosis_z_score_gsva.png",width = 10,height=8)

mat<-raw_TPM
Merged_validation_cohort<-Merged_validation_cohort[order(Merged_validation_cohort$predict.label,decreasing=F),]
Merged_validation_cohort<-Merged_validation_cohort[which(Merged_validation_cohort$predict.label!="unclassified"),]
Merged_validation_cohort$predict.label<-factor(Merged_validation_cohort$predict.label,levels = c("1","2"))
mat<-mat[gen$gene_symbol,Merged_validation_cohort$sample.names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-data.frame(
  molecular_cluster=as.character(Merged_validation_cohort$predict.label),Normal_liver_percentage=Merged_validation_cohort$Normal_liver_percentage,
  Normal_Neuron_percentage=Merged_validation_cohort$Normal_Neuron_percentage,
  Normal_Lymphoid_percentage=Merged_validation_cohort$Normal_Lymphoid_percentage,
  Normal_pancreas_percentage=Merged_validation_cohort$Normal_pancreas_percentage,
  Normal_duodenum_percentage=Merged_validation_cohort$Normal_duodenum_percentage)
rownames(annotation_col)<-Merged_validation_cohort$sample.names
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))
###先z-score 后 scale
pdf(file="./Cohort3_zscore_expression_Ferroptosis.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = T,cluster_rows = T,
         main="Ferroptosis-related gene expression (z-score)",
         annotation_col = annotation_col,
         show_colnames = F,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()

################### MHC molecules ###################
Ferroptosis_related <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/11_additional_exploration/MHC_molecules/MHC molecules.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description=Ferroptosis_related$Category)
mat<-raw_TPM
mat<-mat[gen$gene_symbol,Merged_validation_cohort$sample.names]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-data.frame(
  molecular_cluster=as.character(Merged_validation_cohort$predict.label),Normal_liver_percentage=Merged_validation_cohort$Normal_liver_percentage,
  Normal_Neuron_percentage=Merged_validation_cohort$Normal_Neuron_percentage,
  Normal_Lymphoid_percentage=Merged_validation_cohort$Normal_Lymphoid_percentage,
  Normal_pancreas_percentage=Merged_validation_cohort$Normal_pancreas_percentage,
  Normal_duodenum_percentage=Merged_validation_cohort$Normal_duodenum_percentage)
rownames(annotation_col)<-Merged_validation_cohort$sample.names
annotation_row<-data.frame(Category=gen$gs_description)
rownames(annotation_row)<-gen$gene_symbol
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))
###先z-score 后 scale
library(ComplexHeatmap)
pdf(file="./zscore_expression_molecular_MHCs_antigens.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = T,cluster_rows = T,
         main="MHC_moleculars/co-inhibitor/Co-simulator expression (z-score)",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_colnames = F,
         row_split=annotation_row$Category,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()


#





########################### Chapter 9、prognostics ######
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/9_prognostic/")
Merged_data_tmp_2<-Merged_data
Merged_data_tmp_2<-as.data.frame(Merged_data_tmp_2)
Merged_data_tmp_2$Stage_new<-"IIIA + IIIB + IIIC + IV + IVB"
Merged_data_tmp_2[which(Merged_data_tmp_2$Stage%in%c("I","IA","IB","II","IIA","IIB")),]$Stage_new<-"I + IA/B + II + IIA/B"
Merged_data_tmp_2$Stage_new<-factor(Merged_data_tmp_2$Stage_new,levels = c("I + IA/B + II + IIA/B","IIIA + IIIB + IIIC + IV + IVB"))
Merged_data_tmp_2$Differentiation_new<-"Medium + High"
Merged_data_tmp_2[which(Merged_data_tmp_2$Differentiation%in%c("Medium_low","Low")),]$Differentiation_new<-"Low + Medium_low"
Merged_data_tmp_2$Differentiation_new<-factor(Merged_data_tmp_2$Differentiation_new,levels = c("Low + Medium_low","Medium + High"))
Merged_data_tmp_2$Anno_cluster<-factor(Merged_data_tmp_2$Anno_cluster,levels=c("Metabolic and secretory","mesenchymal and immunosupressive"))


valid_column_names <- make.names(names=names(Merged_data_tmp_2), unique=TRUE, allow_ = TRUE)
names(Merged_data_tmp_2) <- valid_column_names
dfOS = Merged_data_tmp_2[which( !is.na(Merged_data_tmp_2["Survival_status"]) &
                                  !is.na(Merged_data_tmp_2["OS"])&
                                  (Merged_data_tmp_2["OS"] !="Unknown")),]
dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
dfOS$OS = as.numeric(as.character(dfOS$OS))
dependent_os<-"Surv(OS,Survival_status)"
explanatory <- c("Age","Sex","histological_type","Stage_new","Differentiation_new","Anno_cluster")
M<-finalfit(.data =dfOS,dependent =dependent_os,explanatory = explanatory)
write.table(M,file="./COX_model_all.xls",sep = "\t",quote=F)

#univariate 得到 histological type， stage，anno_cluster 显著
#多因素
explanatory <- c("histological_type","Stage_new","Anno_cluster")
M<-finalfit(.data =dfOS,dependent =dependent_os,explanatory = explanatory)
colnames(M)<-c("Overall survival","", "","HR (univariable)","HR (multivariable)")
write.table(M,file="./COX_model_selected.xls",sep = "\t",quote=F)
pdf("./Hazard ratio.pdf",width = 12,height=8)
hr_plot(.data =dfOS,dependent =dependent_os,explanatory = explanatory)
dev.off()


####################### 10. DGE for prognistic signature#######################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/")

####applying singscores
library(readr)
classifier_genes <- read_delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/8_DGE/NTP/filtered_classifier_genes.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
classifier_genes<-classifier_genes[which(classifier_genes$Gene%in%rownames(TCGA_CHOL_RNA_array)),]
classifier_genes<-classifier_genes[c(1:25,150:161),]
table(classifier_genes$Classification)

rankData<-rankGenes(raw_TPM)
GSEA_classifier_genes<-data.frame(gene_symbol=classifier_genes$Gene,gs_description=classifier_genes$Classification)
tgfb_gs_down<-as.character(GSEA_classifier_genes[which(GSEA_classifier_genes$gs_description=="Metabolic and secretory"),]$gene_symbol)
tgfb_gs_up<-as.character(GSEA_classifier_genes[which(GSEA_classifier_genes$gs_description=="mesenchymal and immunosupressive"),]$gene_symbol)
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up,downSet = tgfb_gs_down,knownDirection = T,centerScore = F,dispersionFun = mad)

scoredf<-scoredf[Merged_data$SampleID,]
#scoredf<-scoredf[RNA_Sample[RNA_Sample%notin%Merged_data_tmp_2$SampleID],]
scoredf$DownScore<- -(scoredf$DownScore)
explor_info<-cbind(Merged_data,scoredf)

#explor_info<-cbind(Merged_data[which(Merged_data$SampleID%in%RNA_Sample[RNA_Sample%notin%Merged_data_tmp_2$SampleID]),],scoredf[Merged_data[which(Merged_data$SampleID%in%RNA_Sample[RNA_Sample%notin%Merged_data_tmp_2$SampleID]),]$SampleID,])

ggscatter(explor_info, x = "DownScore", y = "UpScore",
          shape="Anno_cluster",
          color = "TotalScore"  # Add regressin line
) + ggtitle(label="C2_1_score")+ 
  gradient_color(c("blue", "lightgrey", "red"))+
  labs(x = "metabolic_proliferative_score",
       y ="mesenchymal_immunosupressive_socre",
       title ="C2_1_score")
ggsave(filename = "./signature_vs_signature_molecular_cluster.pdf",width = 12,height = 8)

{ggscatter(explor_info, x = "DownScore", y = "UpScore",
           shape="histological_type",
           color = "TotalScore") + ggtitle(label="C2_1_score_for the rest samples")+ 
    gradient_color(c("blue", "lightgrey", "red"))+
    labs(x = "metabolic_proliferative_score",
         y ="mesenchymal_immunosupressive_socre",
         title ="C2_1_score")
  ggsave(filename = "./signature_vs_signature_the_rest_samples.pdf",width = 12,height = 8)
}

ggscatter(explor_info, x = "DownScore", y = "UpScore",
          color = "TotalScore",
          shape="histological_type") +
  labs(x = "metabolic_proliferative_score",
       y ="mesenchymal_immunosupressive_socre",
       title ="C2_1_score")+ 
  gradient_color(c("blue", "lightgrey", "red"))
ggsave(filename = "./signature_vs_signature_histological_cluster.pdf",width = 12,height = 8)

#前25% 和后 25% 做预后
Quantile_score<-quantile(scoredf$TotalScore)
Quantile1_sample<-rownames(scoredf[which(scoredf$TotalScore<=Quantile_score[2]),])
Quantile4_sample<-rownames(scoredf[which(scoredf$TotalScore>=Quantile_score[4]),])
mydata = Merged_data[which(Merged_data$SampleID%in%c(Quantile1_sample,Quantile4_sample)),]
#mydata = explor_info[which(explor_info$SampleID%in%c(Quantile1_sample,Quantile4_sample)),]

mydata$C2_1_score_group<-"Metabolic_proliferative-like (bottom 25%)"
mydata[which(mydata$SampleID%in%Quantile4_sample),]$C2_1_score_group<-"Mesenchymal_immunosupressive-like (top 25%)"
mydata<-as.data.frame(mydata)
groups=c("C2_1_score_group")
for(vals in groups){
  dfDFS = mydata[
    which(
      mydata$SampleID%notin% mydata[which(mydata$Death_reason==" 0"),]$SampleID &
        !is.na(mydata[,vals]) & 
        !is.na(mydata["DFS"]) &
        mydata["DFS"] !="Unknown" &
        !is.na(mydata["Relapse"])),]
  dfDFS$Relapse = as.numeric(as.character(dfDFS$Relapse))
  dfDFS$DFS = as.numeric(as.character(dfDFS$DFS))
  dfDFS["group"] = dfDFS[,vals]
  fitDFS <- survfit(Surv(DFS,Relapse) ~ group,data =dfDFS)
  # coxDFS=coxph(Surv(DFS_MONTHS, DFS_Event) ~group+Stage,data = dfDFS)
  #DFS
  ggsurvplot(fitDFS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="DFS",
             xlab="Time(DAYS)",
             pval.size= 8,pval.coord=c(6,0.3),
             size=1.2
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1),
                    legend.text = c(14,"bold", "darkblue")
    )
  name1 = paste(vals,"DFS.png",sep="_")
  name2 = paste(vals,"DFS.pdf",sep="_")
  ggsave(name1,width=12,height = 8)
  ggsave(name2,width=12,height = 8)
  #OS
  dfOS = mydata[
    which(mydata$SampleID%notin% mydata[which(mydata$Death_reason==" 0"),]$SampleID &
            !is.na(mydata[,vals]) & 
            !is.na(mydata["Survival_status"]) &
            !is.na(mydata["OS"])&
            (mydata["OS"] !="Unknown")),]
  dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
  dfOS$OS = as.numeric(as.character(dfOS$OS))
  dfOS["group"] = dfOS[,vals]
  fitOS <- survfit(Surv(OS,Survival_status) ~ group,data =dfOS) 
  ggsurvplot(fitOS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="OS",
             xlab="Time(DAYS)",
             pval.size= 8,pval.coord=c(6,0.3)
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1)
    )
  name1 = paste(vals,"OS.png",sep="_")
  name2 = paste(vals,"OS.pdf",sep="_")
  ggsave(name1,width=12,height = 8)
  ggsave(name2,width=12,height = 8)
}

#Quantile_score
#    0%       25%       50%       75%      100% 
#0.6519478 1.0579628 1.1708596 1.2830581 1.7222166 

scoredf <- simpleScore(rankData, upSet = tgfb_gs_up,downSet = tgfb_gs_down,knownDirection = T,centerScore = F,dispersionFun = mad)
scoredf$Score_scaled_level<-"low_quantile"
scoredf[which(scoredf$TotalScore< 1.2830581 & scoredf$TotalScore>1.0579628),]$Score_scaled_level<-"median_quantile"
scoredf[which(scoredf$TotalScore>1.2830581),]$Score_scaled_level<-"high_quantile"
scoredf$Score_scaled_level<-factor(scoredf$Score_scaled_level,levels = c("low_quantile","median_quantile","high_quantile"))
scoredf<-scoredf[Merged_data$SampleID,]

explor_info<-cbind(Merged_data,scoredf)
mydata = explor_info
mydata<-as.data.frame(mydata)
dfOS = mydata[which(mydata$SampleID%notin% mydata[which(mydata$Death_reason==" 0"),]$SampleID &!is.na(mydata["Survival_status"]) &
                      !is.na(mydata["OS"])&
                      (mydata["OS"] !="Unknown")),]
dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
dfOS$OS = as.numeric(as.character(dfOS$OS))
fit.coxph <- coxph(Surv(OS,Survival_status)~Sex+Age+histological_type+Score_scaled_level,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest.pdf")


fit.coxph <- coxph(Surv(OS,Survival_status)~Sex+Age+histological_type+TotalScore,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest2.pdf")

############################### 外部数据验证预后 ########################################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/")
tgfb_gs_down<-as.character(GSEA_classifier_genes[which(GSEA_classifier_genes$gs_description=="Metabolic and secretory"),]$gene_symbol)
tgfb_gs_up<-as.character(GSEA_classifier_genes[which(GSEA_classifier_genes$gs_description=="mesenchymal and immunosupressive"),]$gene_symbol)
###################### a, TCGA ###
TCGA_CHOL_clinical <- read.delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/TCGA_CHOL/data_clinical_sample.txt", comment.char="#")
TCGA_CHOL_survival_outcome <- read_delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/TCGA_CHOL/data_clinical_patient.txt", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)
TCGA_CHOL_survival_outcome<-merge(TCGA_CHOL_survival_outcome,TCGA_CHOL_clinical,all.x=T,all.y=T,by.x="PATIENT_ID",by.y="PATIENT_ID")

TCGA_CHOL_RNA <- read_delim("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/TCGA_CHOL/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
TCGA_CHOL_RNA<-TCGA_CHOL_RNA[which(!is.na(TCGA_CHOL_RNA$Hugo_Symbol!="")),-2]
TCGA_CHOL_RNA_array<-TCGA_CHOL_RNA[which(TCGA_CHOL_RNA$Hugo_Symbol%notin%c("ELMOD1", "FGF13", "NKAIN3", "PALM2AKAP2", "QSOX1", "SNAP47", "TMEM8B")),-1]
TCGA_CHOL_RNA_array<-as.data.frame(TCGA_CHOL_RNA_array)
rownames(TCGA_CHOL_RNA_array)<-TCGA_CHOL_RNA[which(TCGA_CHOL_RNA$Hugo_Symbol%notin%c("ELMOD1", "FGF13", "NKAIN3", "PALM2AKAP2", "QSOX1", "SNAP47", "TMEM8B")),]$Hugo_Symbol

rankData<-rankGenes(TCGA_CHOL_RNA_array)
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up[tgfb_gs_up%in%rownames(TCGA_CHOL_RNA_array)],downSet = tgfb_gs_down[tgfb_gs_down%in%rownames(TCGA_CHOL_RNA_array)],knownDirection = T,centerScore = F,dispersionFun = mad)
scoredf$TotalScore
Quantile_score<-quantile(scoredf$TotalScore)
scoredf$Score_scaled_level<-"low_quantile"
scoredf[which(scoredf$TotalScore< 1.1273469 & scoredf$TotalScore>0.9020276),]$Score_scaled_level<-"median_quantile"
scoredf[which(scoredf$TotalScore>1.1273469),]$Score_scaled_level<-"high_quantile"
scoredf$Score_scaled_level<-factor(scoredf$Score_scaled_level,levels = c("low_quantile","median_quantile","high_quantile"))

TCGA_CHOL_survival_outcome$Relapse<-NA
TCGA_CHOL_survival_outcome[which(TCGA_CHOL_survival_outcome$DFS_STATUS=="1:Recurred/Progressed"),]$Relapse<-1
TCGA_CHOL_survival_outcome[which(TCGA_CHOL_survival_outcome$DFS_STATUS=="0:DiseaseFree"),]$Relapse<-0
TCGA_CHOL_survival_outcome$Survival_status<-NA
TCGA_CHOL_survival_outcome[which(TCGA_CHOL_survival_outcome$OS_STATUS=="1:DECEASED"),]$Survival_status<-1
TCGA_CHOL_survival_outcome[which(TCGA_CHOL_survival_outcome$OS_STATUS=="0:LIVING"),]$Survival_status<-0
TCGA_CHOL_clinical$SAMPLE_ID
scoredf<-scoredf[TCGA_CHOL_clinical$SAMPLE_ID,]
TCGA_CHOL_clinical<-cbind(TCGA_CHOL_clinical,scoredf)
rownames(TCGA_CHOL_survival_outcome)<-TCGA_CHOL_survival_outcome$PATIENT_ID
TCGA_CHOL_survival_outcome<-TCGA_CHOL_survival_outcome[TCGA_CHOL_clinical$PATIENT_ID,]
TCGA_CHOL_survival_outcome$Score_scaled_level<-TCGA_CHOL_clinical$Score_scaled_level

mydata = TCGA_CHOL_survival_outcome
mydata<-as.data.frame(mydata)
dfOS = mydata[which(!is.na(mydata["Survival_status"]) &
                      !is.na(mydata["OS_MONTHS"])&
                      (mydata["OS_MONTHS"] !="Unknown")),]
dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
dfOS$OS_MONTHS = as.numeric(as.character(dfOS$OS_MONTHS))
dfOS$TUMOR_TYPE<-factor(dfOS$TUMOR_TYPE,levels=c("Cholangiocarcinoma, Distal","Cholangiocarcinoma, Hilar/Perihilar","Cholangiocarcinoma, Intrahepatic"))
dfOS$AJCC_PATHOLOGIC_TUMOR_STAGE<-factor(dfOS$AJCC_PATHOLOGIC_TUMOR_STAGE,levels=c("STAGE I","STAGE II","STAGE III", "STAGE IVA", "STAGE IVB", "STAGE IV"))
dfOS$SEX<-factor(dfOS$SEX,levels=c("Male","Female"))
dfOS$GRADE<-factor(dfOS$GRADE,levels=c("G1","G2","G3"))
dependent_os<-"Surv(OS_MONTHS,Survival_status)"
explanatory <- c("Score_scaled_level")
M<-finalfit(.data =dfOS,dependent =dependent_os,explanatory = explanatory)
colnames(M)<-c("Overall survival","", "","HR (univariable)","HR (multivariable)")
#write.table(M,file="./COX_model_selected.xls",sep = "\t",quote=F)
hr_plot(.data =dfOS,dependent =dependent_os,explanatory = explanatory)

library(survival)
X<-coxph(Surv(OS_MONTHS,Survival_status) ~  AGE+SEX+Score_scaled_level, data = dfOS)
#%>% summary()
ggforest(X,data=dfOS)


################################ cancer discovery #############################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/Cancer Discovery/")
CD_Clinical_information <- read_excel("Clinical_information.xlsx",sheet = "A. Clinical information")
CD_expression <- read_excel("Clinical_information.xlsx",sheet = "B. RNA-seq")
CD_expression<-CD_expression[,c("Gene",CD_Clinical_information$Sample_ID)]
mat<-CD_expression
mat<-as.data.frame(mat)
mat<-mat[which(mat$Gene%notin%c("44256","44257")),]
rownames(mat)<-mat$Gene
mat<-mat[,-1]
rankData<-rankGenes(mat)
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up[tgfb_gs_up%in%rownames(mat)],downSet = tgfb_gs_down[tgfb_gs_down%in%rownames(mat)],knownDirection = T,centerScore = F,dispersionFun = mad)

Quantile_score<-quantile(scoredf$TotalScore)
scoredf$Score_scaled_level<-"low_quantile"
scoredf[which(scoredf$TotalScore< 1.1623278 & scoredf$TotalScore>0.9571616),]$Score_scaled_level<-"median_quantile"
scoredf[which(scoredf$TotalScore>1.1623278),]$Score_scaled_level<-"high_quantile"
scoredf$Score_scaled_level<-factor(scoredf$Score_scaled_level,levels = c("low_quantile","median_quantile","high_quantile"))


scoredf<-scoredf[CD_Clinical_information$Sample_ID,]
CD_Clinical_information$TotalScore<-scoredf$TotalScore
CD_Clinical_information$Score_scaled_level<-scoredf$Score_scaled_level
CD_Clinical_information$Survival_status<-CD_Clinical_information$Survival
CD_Clinical_information$OS_MONTHS<-CD_Clinical_information$OS
mydata = CD_Clinical_information
mydata<-as.data.frame(mydata)
dfOS = mydata[which(!is.na(mydata["Survival_status"]) &
                      !is.na(mydata["OS_MONTHS"])&
                      (mydata["OS_MONTHS"] !="Unknown")),]
dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
dfOS$OS_MONTHS = as.numeric(as.character(dfOS$OS_MONTHS))
dfOS$SEX<-factor(dfOS$Sex,levels=c("Male","Female"))


fit.coxph <- coxph(Surv(OS_MONTHS,Survival_status)~Age+SEX+Score_scaled_level,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest.pdf")

fit.coxph <- coxph(Surv(OS_MONTHS,Survival_status)~Age+SEX+TotalScore,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest2.pdf")

############# GSE89749 #############################################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/GSE89749/")
library(GEOquery)
library(limma)
library(Biobase)
library(umap)
# load series and platform data from GEO
load("/Users/zouxinchen/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/GSE89749/zscore_merged.RData")
GSE89749_clinical <- read_excel("GSE89749_clinical.xlsx")
ex<-ex[,-1]
#ex<-ex[which(ex$Gene%in%c(tgfb_gs_down,tgfb_gs_up)),]
ex<-ex[,-121]
X<-as.data.frame(table(ex$Gene))
raw_TPM_1<-ex[ex$Gene%in%as.character(X[which(X$Freq==1),1]),]
raw_TPM_2<-ex[ex$Gene%in%as.character(X[which(X$Freq>1),1]),]
temperate<-matrix(NA,nrow=length(unique(raw_TPM_2$Gene)),ncol = ncol(ex)-1)
for (i in 1:length(unique(raw_TPM_2$Gene))) {
  temperate[i,]<-colMeans(raw_TPM_2[which(raw_TPM_2$Gene==unique(raw_TPM_2$Gene)[i]),-121])
}
temperate<-as.data.frame(temperate)
colnames(temperate)<-colnames(ex)[-121]
temperate$Gene<-unique(raw_TPM_2$Gene)
ex<-rbind(raw_TPM_1,temperate)
ex_gene<-ex$Gene
ex<-ex[,-121]
ex<-as.data.frame(ex)
rownames(ex)<-ex_gene
rankData<-rankGenes(ex)
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up[tgfb_gs_up%in%rownames(ex)],downSet = tgfb_gs_down[tgfb_gs_down%in%rownames(ex)],knownDirection = T,centerScore = F,dispersionFun = mad)

Quantile_score<-quantile(scoredf$TotalScore)
scoredf$Score_scaled_level<-"low_quantile"
scoredf[which(scoredf$TotalScore< 1.0683572 & scoredf$TotalScore>0.9156521),]$Score_scaled_level<-"median_quantile"
scoredf[which(scoredf$TotalScore>1.0683572),]$Score_scaled_level<-"high_quantile"
scoredf$Score_scaled_level<-factor(scoredf$Score_scaled_level,levels = c("low_quantile","median_quantile","high_quantile"))
scoredf<-scoredf[GSE89749_clinical$`Sample ID`,]

GSE89749_clinical$Survival_status<-GSE89749_clinical$`Vital state (1=Dead)`
GSE89749_clinical$OS_MONTHS<-GSE89749_clinical$`Overall survival (days)`
GSE89749_clinical$Score_scaled_level<-scoredf$Score_scaled_level
GSE89749_clinical$TotalScore<-scoredf$TotalScore


Quantile_score<-quantile(scoredf$TotalScore)
Quantile1_sample<-rownames(scoredf[which(scoredf$TotalScore<=Quantile_score[2]),])
Quantile4_sample<-rownames(scoredf[which(scoredf$TotalScore>=Quantile_score[4]),])
mydata = GSE89749_clinical
mydata$C2_1_score_group<-"Metabolic_proliferative-like (bottom 25%)"
mydata[which(mydata$`Sample ID`%in%Quantile4_sample),]$C2_1_score_group<-"Mesenchymal_immunosupressive-like (top 25%)"

mydata<-as.data.frame(mydata)
groups=c("C2_1_score_group")
for(vals in groups){
  dfOS = mydata[which( mydata["Survival_status"]!="N/A" &
                         !is.na(mydata["Survival_status"]) &
                         !is.na(mydata["OS_MONTHS"])&
                         (mydata["OS_MONTHS"] !="N/A")),]
  dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
  dfOS$OS_MONTHS = as.numeric(as.character(dfOS$OS_MONTHS))
  dfOS["group"] = dfOS[,vals]
  fitOS <- survfit(Surv(OS_MONTHS,Survival_status) ~ group,data =dfOS) 
  ggsurvplot(fitOS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="OS",
             xlab="Time(DAYS)",
             pval.size= 8,pval.coord=c(6,0.3)
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1)
    )
  name1 = paste(vals,"OS.png",sep="_")
  name2 = paste(vals,"OS.pdf",sep="_")
  ggsave(name1,width=12,height = 8)
  ggsave(name2,width=12,height = 8)
}



mydata = GSE89749_clinical
mydata<-as.data.frame(mydata)
dfOS = mydata[which( mydata["Survival_status"]!="N/A" &
                       !is.na(mydata["Survival_status"]) &
                       !is.na(mydata["OS_MONTHS"])&
                       (mydata["OS_MONTHS"] !="N/A")),]
dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
dfOS$OS_MONTHS = as.numeric(as.character(dfOS$OS_MONTHS))
dfOS$Sex<-factor(dfOS$Sex,levels=c("M","F"))
dfOS$`Anatomical subtype`<-factor(dfOS$`Anatomical subtype`,levels=c("Intrahepatic","Perihilar","Distal"))
dfOS$Stage<-factor(dfOS$Stage,levels=c("0","I","IA","II","IIA","IIB", "III","IIIA","IIIB","IV" ,"IVA","IVB"))
colnames(dfOS)<-c("Sample_ID","Sex","Age","Anatomical_subtype","TNM_staging","Stage","Vital_state","Overall_survival","Gene_expression_profiling","Survival_status","OS_MONTHS","Score_scaled_level","Totalscore")
dfOS$Age<-as.numeric(dfOS$Age)
dependent_os<-"Surv(OS_MONTHS,Survival_status)"
explanatory <- colnames(dfOS)[c(2,3,4,6,12,13)]
M<-finalfit(.data =dfOS,dependent =dependent_os,explanatory = explanatory)
colnames(M)<-c("Overall survival","", "","HR (univariable)","HR (multivariable)")
#write.table(M,file = "single_factor_survival.txt",sep="\t",quote=F)
hr_plot(.data =dfOS,dependent =dependent_os,explanatory = explanatory)

fit.coxph <- coxph(Surv(OS_MONTHS,Survival_status)~Sex+Age+Anatomical_subtype+Score_scaled_level,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest.pdf")

fit.coxph <- coxph(Surv(OS_MONTHS,Survival_status)~Sex+Age+Anatomical_subtype+Totalscore,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest2.pdf")


############################ cohort3 ####################
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/Cohort3/")
cohort_3_tpm <- read_xlsx("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/Cohort3/Cohort3_TPM.xlsx")
cohort3_protein <- read_excel("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/external_validation/Cohort3/protein_expression.xlsx")
{
  mat<-cohort_3_tpm
  mat[,-1]<- (2^cohort_3_tpm[,-1])-1
  mat<-mat[which(rowSums(mat[,-1])>0),]
  Gene_tmp<-mat$NAME
  mat<-mat[,-1]
  rownames(mat)<-Gene_tmp
  mat<-as.data.frame(mat)
  rownames(mat)<-Gene_tmp
  mat<-mat[which(!is.na(rowSums(mat))),]
  colnames(mat)<-paste0("s",colnames(mat))
}

library(ConsensusClusterPlus)
clinical_information_cohort3 <- read_excel("clinical_information.xlsx")
Merged_Cohort3<-clinical_information_cohort3
Merged_Cohort3<-Merged_Cohort3[which(Merged_Cohort3$Patient_ID%in%colnames(mat)),]
####################### protein expression
mat_2<-cohort3_protein
Gene_tmp<-mat_2$Protein_NAME
mat_2<-mat_2[,-1]
rownames(mat_2)<-Gene_tmp
mat_2 <- as.data.frame(sapply(mat_2, as.numeric))
rownames(mat_2)<-Gene_tmp
#mat<-mat[which(!is.na(rowSums(mat))),]

########################################
rankData<-rankGenes(mat)
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up[tgfb_gs_up%in%rownames(mat)],downSet = tgfb_gs_down[tgfb_gs_down%in%rownames(mat)],knownDirection = T,centerScore = F,dispersionFun = mad)
Quantile_score<-quantile(scoredf$TotalScore)
scoredf$Score_scaled_level<-"low_quantile"
Quantile_score
scoredf[which(scoredf$TotalScore< 1.0633221 & scoredf$TotalScore>0.8817689),]$Score_scaled_level<-"median_quantile"
scoredf[which(scoredf$TotalScore>1.0633221),]$Score_scaled_level<-"high_quantile"
scoredf$Score_scaled_level<-factor(scoredf$Score_scaled_level,levels = c("low_quantile","median_quantile","high_quantile"))
scoredf<-scoredf[Merged_Cohort3$Patient_ID,]
Merged_Cohort3$Survival_status<-Merged_Cohort3$OS_STATUS
Merged_Cohort3$OS_MONTHS<-Merged_Cohort3$OS_DAYS/30
Merged_Cohort3$Score_scaled_level<-scoredf$Score_scaled_level
Merged_Cohort3$TotalScore<-scoredf$TotalScore
mydata = Merged_Cohort3
mydata<-as.data.frame(mydata)
dfOS = mydata[which( mydata["Survival_status"]!="N/A" &
                       !is.na(mydata["Survival_status"]) &
                       !is.na(mydata["OS_MONTHS"])&
                       (mydata["OS_MONTHS"] !="N/A")),]
dfOS$Survival_status = as.numeric(as.character(dfOS$Survival_status))
dfOS$OS_MONTHS = as.numeric(as.character(dfOS$OS_MONTHS))
dfOS$Sex<-factor(dfOS$Sex,levels=c("Male","Female"))
dfOS$Stage<-factor(dfOS$`TNM stage`,levels=c("IA","IB","II","IIIA","IIIB","IV"))
dfOS$`Adjuvant therapies`<-factor(dfOS$`Adjuvant therapies`,levels=c("No treatment","Chemotherapy","TACE"))
colnames(dfOS)<-c("Patient","Patient_ID","Sex","Age","Intrahepatic_metastasis","Liver_fluke","HBsAg","Biliary_tract_stone_disease","Tumor_size_diameter","Vascular_invasion","Liver_cirrhosis","Regional_lymph_Node_metastasis","Distal_metastasis" , "Perineural_invasion","TB" ,"ALB","Preoperative_serum_AFP", "CA19_9","CEA","ALT","γ-GT","TNM_stage","OS_DAYS","OS_STATUS","Adjuvant_therapies","WES_seq","RNA_seq" ,"Proteome","Phosphoproteome","Microbiome","Survival_status",        "OS_MONTHS", "Score_scaled_level" ,"TotalScore","Stage")
dfOS$Age<-as.numeric(dfOS$Age)
dependent_os<-"Surv(OS_MONTHS,Survival_status)"
explanatory <- colnames(dfOS)[c(3,4,33,35)]
M<-finalfit(.data =dfOS,dependent =dependent_os,explanatory = explanatory)
colnames(M)<-c("Overall survival","", "","HR (univariable)","HR (multivariable)")
#write.table(M,file = "single_factor_survival.txt",sep="\t",quote=F)
hr_plot(.data =dfOS,dependent =dependent_os,explanatory = explanatory)

fit.coxph <- coxph(Surv(OS_MONTHS,Survival_status)~Sex+Age+Stage+Score_scaled_level,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest.pdf")

fit.coxph <- coxph(Surv(OS_MONTHS,Survival_status)~Sex+Age+Stage+TotalScore,data =dfOS)
ggforest(fit.coxph, data = dfOS)+
  theme(axis.text.x = element_text(face="bold", color="black", size=30),title =element_text(face="bold", color="black",size=30))
ggsave(filename = "./ggforest2.pdf")


groups=c("Score_scaled_level")
for(vals in groups){
  dfOS = mydata[which(!is.na(mydata[,vals]) & 
                        mydata$Score_scaled_level!="median_quantile"&
                        !is.na(mydata["OS_MONTHS"]) &
                        !is.na(mydata["OS_STATUS"])&
                        (mydata["OS_STATUS"] !="Unknown")),]
  dfOS$OS_STATUS = as.numeric(as.character(dfOS$OS_STATUS))
  dfOS$OS_MONTHS = as.numeric(as.character(dfOS$OS_MONTHS))
  dfOS["group"] = dfOS[,vals]
  fitOS <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ group,data =dfOS) 
  ggsurvplot(fitOS,pval = TRUE, conf.int = FALSE,linetype = "strata",
             risk.table = FALSE,risk.table.y.text.col = TRUE,
             palette = "ucscgb",
             legend.title="",
             ylab="OS",
             xlab="Time(Month)",
             pval.size= 8,pval.coord=c(6,0.3)
  ) +
    theme_survminer(font.main = c(15, "bold", "darkblue"),
                    font.submain = c(15, "bold", "black"),
                    font.caption = c(15, "plain", "black"),
                    font.x = c(20, "bold", "black"),
                    font.y = c(20, "bold", "black"),
                    font.tickslab = c(18, "bold", "black"),
                    base_size = 18,
                    axis.line = element_line(size=1),
                    axis.ticks = element_line(size=1)
    )
  name1 = paste(vals,"OS.png",sep="_")
  name2 = paste(vals,"OS.pdf",sep="_")
  ggsave(name1,width=12,height = 8)
  ggsave(name2,width=12,height = 8)
}


############################## 11、Additional exploration ###########
########################1, metabolic & proliferative-c2 ###################
#bile acid metabolism
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/11_additional_exploration/")
tpm_exp<-data_for_nmf[,-1]
rownames(tpm_exp)<-data_for_nmf$...1
bile_acid_genes<-c("CYP7B1","CYP27A1","CYP7A1","CYP8B1","NR1H4")
bile_acid_genes<-c("EPCAM","KRT1")
bile_acid_genes<-c("SREBF1","SCD","ACLY", "ACACA", "FASN","PKLR","PCK1","G6PC","G6PD")
bile_acid_genes<-c("CDH1","CDH2","MTOR","GPX4")
bile_acid_genes<-c("MTOR","SREBF1")

enzyme_exp<-tpm_exp[bile_acid_genes,]
rownames(enzyme_exp)<-bile_acid_genes
enzyme_exp<-t(enzyme_exp)
enzyme_exp<-enzyme_exp[Merged_data$SampleID,]
explor_info<-cbind(Merged_data,enzyme_exp)

for (i in bile_acid_genes) {
  ggplot(explor_info,aes(x=explor_info$Anno_cluster,
                         y=log2(explor_info[,i]),
                         color=Anno_cluster))+
    geom_boxplot(notch = F)+
    labs(x = "molecular subtype",
         y = paste0("log2(tpm) in ",i),
         title =paste0("log2 expression of ",i," in each molecular subtype"))+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
    stat_compare_means(label.y = 10)+
    theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10)) 
  ggsave(filename = paste0(i,".png"),width = 10,height=8)
}


#### hepatocyte
hepatocyte_genes<-c("ALB","TAT","TF","CYP3A4","KRT8","KRT18")
enzyme_exp<-tpm_exp[hepatocyte_genes,]
rownames(enzyme_exp)<-hepatocyte_genes
enzyme_exp<-t(enzyme_exp)
enzyme_exp<-enzyme_exp[Merged_data$SampleID,]
explor_info<-cbind(Merged_data,enzyme_exp)
for (i in hepatocyte_genes) {
  ggplot(explor_info,aes(x=factor(explor_info$Anno_cluster,levels = c("mesenchymal and immunosupressive","Metabolic and proliferative"))
                         ,y=log2(explor_info[,i]),
                         color=Anno_cluster))+
    geom_boxplot(notch = F)+
    labs(x = "molecular subtype",
         y = paste0("log2(tpm) in ",i),
         title =paste0("log2 expression of ",i," in each molecular subtype"))+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
    stat_compare_means(label.y = 10)+
    theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
  ggsave(filename = paste0(i,".png"),width = 10,height=8)
}

## Insulin
hepatocyte_genes<-c("AR","INS")
enzyme_exp<-tpm_exp[hepatocyte_genes,]
rownames(enzyme_exp)<-hepatocyte_genes
enzyme_exp<-t(enzyme_exp)
enzyme_exp<-enzyme_exp[Merged_data$SampleID,]
explor_info<-cbind(Merged_data,enzyme_exp)

my_comparisons = list(c("mesenchymal and immunosupressive","Metabolic and proliferative"))
label.y = c(7)
for (i in hepatocyte_genes) {
  ggplot(explor_info,aes(x=factor(explor_info$Anno_cluster,levels = c("mesenchymal and immunosupressive","Metabolic and proliferative"))
                         ,y=log2(explor_info[,i]),
                         color=Anno_cluster))+
    geom_boxplot(notch = F)+
    labs(x = "molecular subtype",
         y = paste0("log2 tpm in ",i),
         title =paste0("log2 expression of ",i," in each molecular subtype"))+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
    stat_compare_means(label.y = 10)+
    stat_compare_means(comparisons = my_comparisons, label.y = label.y, method = "t.test")+theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
  ggsave(filename = paste0(i,".png"),width = 10,height=8)
} 

#androgen 性别，stage
adrogen_tmp<-explor_info[which(explor_info$Anno_cluster=="Metabolic and proliferative"),]
ggplot(adrogen_tmp,aes(x=Sex,y=log2(adrogen_tmp[,"AR"]),
                       color=Sex))+
  geom_boxplot(notch = F)+
  labs(x = "Sex",
       y = paste0("log2(tpm) in ","AR"),
       title =paste0("log2 expression of ","AR"," each Sex"))+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
  stat_compare_means(label.y = 2)+
  theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
ggsave(filename = "Sex_androgen.png",width = 10,height=8)

################ Ferroptosis resistance ###############
Ferroptosis_related <- read_excel("ferroptosis resistance/Ferroptosis_related.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description="Ferroptosis-related")
mat<-tpm_exp
mat<-as.data.frame(mat)
rownames(mat)<-rownames(tpm_exp)
mat<-mat[which(rownames(mat)%in%gen$gene_symbol),]
mat<-mat[which(!is.na(rowSums(mat))),]
mat
gen
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}

mat=as.matrix(mat)
es <- gsva(mat, gs, method="zscore",verbose=FALSE, parallel.sz=2)
es<-t(es)
es<-as.data.frame(es)
es<-es[Merged_data$SampleID,]
explor_info<-cbind(Merged_data,es)
my_comparisons = list(c("mesenchymal and immunosupressive","Metabolic and secretory"))
label.y = c(1)

ggplot(explor_info,aes(x=factor(explor_info$Anno_cluster,levels = c("mesenchymal and immunosupressive","Metabolic and proliferative"))
                       ,y=log2(explor_info$es),
                       color=Anno_cluster))+
  geom_boxplot(notch = F)+
  labs(x = "molecular subtype",
       y = paste0("log2 enrichment score: Ferroptosis"),
       title =paste0("log2 enrichment of Ferroptosis in each molecular subtype"))+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
  stat_compare_means(label.y = 4)+
  stat_compare_means(comparisons = my_comparisons, label.y = 5, method ="t.test")+
  theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
ggsave(filename ="./ferroptosis resistance/Ferroptosis_z_score_gsva.png",width = 10,height=8)

pdf(file = paste0("./ferroptosis resistance/Ferroptosis_z_score_gsva.pdf"),width = 10,height=8)
ggboxplot(explor_info, 
          x = "Anno_cluster", y = "es", 
          color = "Anno_cluster", 
          ylab = "Ferroptosis signature Score (GSEA:z-score)", 
          xlab = "Molecular subtype")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means()
dev.off()



Merged_data<-Merged_data[order(Merged_data$Anno_cluster,decreasing=T),]
mat<-mat[,Merged_data$SampleID]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-data.frame(molecular_cluster=Merged_data$Anno_cluster)
rownames(annotation_col)<-Merged_data$SampleID
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))
###先z-score 后 scale
pdf(file="./zscore_expression_Ferroptosis.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = F,cluster_rows = F,
         main="Ferroptosis-related gene expression (z-score)",
         annotation_col = annotation_col,
         show_colnames = F,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()


############## MHC molecules, co-stimulators and co-inhibitors ################
Ferroptosis_related <- read_excel("./MHC_molecules/MHC molecules.xlsx")
gen<-data.frame(gene_symbol=Ferroptosis_related$GeneID,gs_description=Ferroptosis_related$Category)
mat<-tpm_exp
mat<-as.data.frame(mat)
rownames(mat)<-rownames(tpm_exp)
mat<-mat[which(rownames(mat)%in%gen$gene_symbol),]
mat<-mat[which(!is.na(rowSums(mat))),]

Merged_data<-Merged_data[order(Merged_data$Anno_cluster,decreasing=T),]
mat<-mat[gen$gene_symbol,Merged_data$SampleID]
mat<-t(scale(t(mat), scale=TRUE, center=TRUE))
annotation_col<-data.frame(molecular_cluster=Merged_data$Anno_cluster)
rownames(annotation_col)<-Merged_data$SampleID
annotation_row<-data.frame(Category=gen$gs_description)
rownames(annotation_row)<-gen$gene_symbol
bk<-c(seq(min(mat),0,length.out = 50),seq(0.00001,max(mat),length.out = 50))

###先z-score 后 scale
library(ComplexHeatmap)
pdf(file="./zscore_expression_molecular_MHCs_antigens.pdf",width = 12,height = 10)
pheatmap(as.matrix(mat),
         cluster_cols = F,cluster_rows = F,
         main="MHC_moleculars/co-inhibitor/Co-simulator expression (z-score)",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_colnames = F,
         row_split=annotation_row$Category,
         column_split=annotation_col$molecular_cluster,
         color =c(colorRampPalette(colors = c("Navyblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), breaks = bk)
dev.off()








########################## Mesenchymal & Immunosupressive-C1 #####################
## hypoxia
hypoxia_gene <- read_csv("hypoxia_gene.csv")
gen<-data.frame(gene_symbol=hypoxia_gene$gene,gs_description=hypoxia_gene$term)
mat<-tpm_exp
mat<-as.data.frame(mat)
rownames(mat)<-rownames(tpm_exp)
mat<-mat[which(rownames(mat)%in%gen$gene_symbol),]
mat<-mat[which(!is.na(rowSums(mat))),]
mat
gen
module=levels(as.factor(gen$gs_description))
len=length(module)
gs=list()
for(y in 1:len){
  gs[[module[y]]]<-subset(gen,gen[,2]==module[y])[,1]
}
mat=as.matrix(mat)
es <- gsva(mat, gs, method="ssgsea", verbose=FALSE, parallel.sz=2)
es<-t(es)
es<-as.data.frame(es)
es<-es[Merged_data$SampleID,]
explor_info<-cbind(Merged_data,es)
my_comparisons = list(c("mesenchymal and immunosupressive","Metabolic and secretory"))
label.y = c(1)
explor_info$Anno_cluster<-factor(explor_info$Anno_cluster,levels = c("mesenchymal and immunosupressive","Metabolic and secretory"))

for (i in c("Buffa","West","Winter")) {
  ggplot(explor_info,aes(x=factor(explor_info$Anno_cluster,levels = c("mesenchymal and immunosupressive","Metabolic and secretory"))
                         ,y=log2(explor_info[,i]),
                         color=Anno_cluster))+
    geom_boxplot(notch = F)+
    labs(x = "molecular subtype",
         y = paste0("log2 expression: ",i),
         title =paste0("log2 exp of ",i," in each molecular subtype"))+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
    stat_compare_means(label.y = 2)+
    stat_compare_means(comparisons = my_comparisons, label.y = label.y, method ="t.test")+theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10)) ggsave(filename = paste0(i,".png"),width = 10,height=8)
}


pdf(file = "./hypoxia/West_hypoxia.pdf",width = 10,height=8)
A<-ggboxplot(explor_info, 
             x = "Anno_cluster", y = "West", 
             color = "Anno_cluster", 
             ylab = "West hypoxia signature score", xlab = "Molecular subtype")+
  stat_summary(fun="mean",color="black")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means(label.y = 1)+
  ggtitle(label="West: hypoxia signature score vs molecular subtype")

print(A)
dev.off()



#### Wnt-beta catenin signaling
hepatocyte_genes<-c("SNAI1","IL1B","IL6")
enzyme_exp<-tpm_exp[hepatocyte_genes,]
rownames(enzyme_exp)<-hepatocyte_genes
enzyme_exp<-t(enzyme_exp)
enzyme_exp<-enzyme_exp[Merged_data$SampleID,]
explor_info<-cbind(Merged_data,enzyme_exp)
label.y = c(6)

for (i in hepatocyte_genes) {
  ggplot(explor_info,aes(x=factor(explor_info$Anno_cluster,levels = c("mesenchymal and immunosupressive","Metabolic and proliferative"))
                         ,y=log2(explor_info[,i]),
                         color=Anno_cluster))+
    geom_boxplot(notch = F)+
    labs(x = "molecular subtype",
         y = paste0("log2 expression: ",i),
         title =paste0("log2 exp of ",i," in each molecular subtype"))+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
    stat_compare_means(label.y = 5)+
    stat_compare_means(comparisons = my_comparisons, label.y = label.y, method ="t.test")+
    theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
  ggsave(filename = paste0(i,".png"),width = 10,height=8)
}

# EMT 
EMT_gene <- read_csv("EMT_gene.csv")
# Antigen process
EMT_gene<- read_csv("Antigen_process_machinary.csv")
EMT_gene<-c(EMT_gene$genes)
#IL6 _STAT3_inflammatory 
EMT_gene<-c("STAT3","JAK1","JAK2","IL6")
# angiogenesis
EMT_gene<-c("VEGFA","VEGFC","HIF1A")
EMT_gene<-c("SNAI1","IL1B")

es<-tpm_exp[EMT_gene,]
rownames(es)<-EMT_gene
es<-t(es)
es<-as.data.frame(es)
es<-es[Merged_data$SampleID,]
es<-as.data.frame(es)
explor_info<-cbind(Merged_data,es)

for (i in unique(EMT_gene)) {
  ggplot(explor_info,aes(x=factor(explor_info$Anno_cluster,levels = c("mesenchymal and immunosupressive","Metabolic and secretory"))
                         ,y=log2(explor_info[,i]),
                         color=Anno_cluster))+
    geom_boxplot(notch = F)+
    labs(x = "molecular subtype",
         y = paste0("log2 expression: ",i),
         title =paste0("log2 exp of ",i," in each molecular subtype"))+
    stat_summary(fun="mean",color="black")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.3)+
    stat_compare_means(label.y = 10)+
    stat_compare_means(comparisons = my_comparisons, label.y = 10, method ="t.test")+
    theme(axis.text.x =element_text(color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=14),title =element_text(face="bold", color="black", size=10))
  ggsave(filename = paste0(i,".png"),width = 10,height=8)
}


for (i in unique(EMT_gene)) {
  A<-ggboxplot(explor_info, 
               x = "Anno_cluster", y = print(i), 
               color = "Anno_cluster", 
               ylab = paste0("Expression of ",print(i)), xlab = "Molecular subtype")+
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
    stat_compare_means()+
    ggtitle(label=paste0("Expression of ",print(i)," vs molecular subtype"))
  pdf(file = paste0("./Wnt_beta//",print(i),".pdf"),width = 10,height=8)
  print(A)
  dev.off()
}


###############TIDE immune
raw_TPM_tmp<-raw_TPM[,Merged_data$SampleID]
raw_TPM_tmp<-log2(raw_TPM_tmp+1)
raw_TPM_tmp_means<-rowMeans(raw_TPM_tmp)
A<-matrix(as.numeric(rep(raw_TPM_tmp_means,164)),nrow=55880,ncol=164)
raw_TPM_tmp<-raw_TPM_tmp-A
#write.table(raw_TPM_tmp,file = "./11_additional_exploration/immunosupressive/tpm_normalized_forTIDE.txt",sep = "\t",quote=F,row.names = T)
#############去 web
TIDE_immune <- read_csv("./11_additional_exploration/immunosupressive/TIDE_immune.csv")
TIDE_immune<-TIDE_immune[which(TIDE_immune$Patient%in%Merged_data$SampleID),]
Merged_data<-merge(Merged_data,TIDE_immune,by.x="SampleID",by.y="Patient",all.x=T,all.y=T)
#TIDE,Dysfunction,Exclusion
pdf(file = paste0("./11_additional_exploration/immunosupressive/","TIDE",".pdf"),width = 10,height=8)
ggboxplot(Merged_data, 
          x = "Anno_cluster", y = "TIDE", 
          color = "Anno_cluster", 
          ylab = "TIDE Score", xlab = "Molecular subtype")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means()+
  ggtitle(label=paste0("TIDE score vs molecular subtype"))
dev.off()

pdf(file = paste0("./11_additional_exploration/immunosupressive/","T cell dysfunction",".pdf"),width = 10,height=8)
ggboxplot(Merged_data, 
          x = "Anno_cluster", y = "Dysfunction", 
          color = "Anno_cluster", 
          ylab = "T cell dysfunction Score", xlab = "Molecular subtype")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means()+
  ggtitle(label=paste0("T cell dysfunction score vs molecular subtype"))
dev.off()


pdf(file = paste0("./11_additional_exploration/immunosupressive/","T cell Exclusion",".pdf"),width = 10,height=8)
ggboxplot(Merged_data, 
          x = "Anno_cluster", y = "Exclusion", 
          color = "Anno_cluster", 
          ylab = "T cell Exclusion Score", xlab = "Molecular subtype")+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.4)+
  stat_compare_means()+
  ggtitle(label=paste0("T cell exclusion score vs molecular subtype"))
dev.off()


################################### Tree plot #####
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/1_NMF/treeplot/")
View(raw_TPM)
Sample_analysising
all_sample<-Merged_data$SampleID
USArrests<-raw_TPM[,all_sample]
dists <- dist(USArrests,method = "euclidean") 

# 进行层次聚类，method = "average"选择UPGMA聚类算法
hc <- hclust(dists, method = "ave")
hc
## 
## Call:
## hclust(d = dists, method = "ave")
## 
## Cluster method   : average 
## Distance         : euclidean 
## Number of objects: 50

# 将hclust对象转换为dendrogram对象
dend1 <- as.dendrogram(hc)
dend1
## 'dendrogram' with 2 branches and 50 members total, at height 152.314

# 绘制聚类树图，默认type = "rectangle"
plot(dend1, type = "rectangle", 
     ylab="Height",
     main="Cluster Dendrogram")







setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/10_integrated_signature/ROC_NRI/")
### ROC 曲线
library(rms)
library(timeROC)
View(scoredf)
View(Merged_data)
scoredf<-scoredf[Merged_data$SampleID,]

{
  ROC_data<-data.frame(time=Merged_data$OS,status=Merged_data$Survival_status,riskScore=scoredf$TotalScore,SampleID=Merged_data$SampleID,Death_reason=Merged_data$Death_reason,staging=Merged_data$Stage,anatomical_location=Merged_data$histological_type)
  ROC_data$staging_category<-0
  ROC_data[which(ROC_data$staging%in%c("I","IA","IB")),]$staging_category<-1
  ROC_data[which(ROC_data$staging%in%c("II","IIA","IIB")),]$staging_category<-2
  ROC_data[which(ROC_data$staging%in%c("IIIA","IIIB","IIIC")),]$staging_category<-3
  ROC_data[which(ROC_data$staging%in%c("IV","IVA","IVB")),]$staging_category<-4
  
  #intrahepatic #hilar #extrahepatic
  ROC_data<-ROC_data[which(ROC_data$anatomical_location=="hilar"),]
  ROC_data<-ROC_data[which(!is.na(ROC_data$time) & ROC_data$time!="Unknown" & ROC_data$SampleID%notin% ROC_data[which(ROC_data$Death_reason==" 0"),]$SampleID),]
}

ROC_data2<-data.frame(time=Merged_Cohort3$OS_DAYS,status=Merged_Cohort3$OS_STATUS,riskScore=Merged_Cohort3$TotalScore,SampleID=Merged_Cohort3$Patient_ID)
ROC_data2<-ROC_data2[which(!is.na(ROC_data2$time) &!is.na(ROC_data2$status)),]

ROC_data3<-data.frame(time=scoredf_GSE89749$`Overall survival (days)`,status=scoredf_GSE89749$Survival_status,riskScore=scoredf_GSE89749$TotalScore,SampleID=scoredf_GSE89749$`Sample ID`)
ROC_data3<-ROC_data3[which(ROC_data3$time!="N/A" &ROC_data3$status!="N/A"),]

#ROC_data_combined<-rbind(ROC_data,ROC_data2,ROC_data3)

ROC <- timeROC(T= (as.numeric(ROC_data$time)/(365.5/12)), #将结局时间转换为月
               delta = as.numeric(ROC_data$status), #生存结局
               marker = ROC_data$staging_category, #预测变量 
               cause = 1, #阳性结局赋值
               weighting = "marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times = c(12,24,36,48),
               iid = T) #只有marginal可以计算置信区间 
ROC


pdf(file = "./ROC_TNM_pCCA.pdf")
{
  plot(ROC, time=12*1, col = "blue", add = F, title = F)
  plot(ROC, time=12*2, col = "black", add = T)
  plot(ROC, time=12*3, col = "green", add = T)
  #plot(ROC, time=12*4, col = "red", add = T)
  legend("bottomright", lty = 1, cex = 1.0,
         col = c("blue", "black", "green","red"),
         legend = c("1 year AUC:52.66", 
                    "2 year AUC:55.52",
                    "3 year AUC:56.58 ",
                    "4 year AUC:75.92"))
}
dev.off()



#treeplot
setwd("~/Desktop/GenePlus/project/胆管癌/2022_12_28/1_NMF/treeplot/")
Merged_data[which(Merged_data$histological_type=="Distal"),]$histological_type<-"dCCA"
Merged_data[which(Merged_data$histological_type=="perihilar"),]$histological_type<-"pCCA"
Merged_data[which(Merged_data$histological_type=="intrahepatic"),]$histological_type<-"iCCA"
USArrests<-raw_TPM[Protein_coding_genes,Merged_data$SampleID]
colnames(USArrests)<-Merged_data$histological_type
USArrests<-USArrests[which(rowSums(USArrests)>0),]
USArrests<-t(USArrests)
dists <- dist(USArrests,method = "euclidean")
head(dists)
# 进行层次聚类，method = "average"选择UPGMA聚类算法
hc <- hclust(dists, method = "ave")
ggdendrogram(hc)+ geom_text(data=label(ddata_x),
                            aes(label=labs$group,x=x, y=0, colour=labs$group),position ="identity")

hc<-as.dendrogram(hc)
ddata_x <- dendro_data(hc)
p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
labs <- label(ddata_x)
labs$group <- NA
labs[which(labs$label=="dCCA"),]$group<-"D"
labs[which(labs$label=="iCCA"),]$group<-"I"
labs[which(labs$label=="pCCA"),]$group<-"P"
p2 + geom_text(data=label(ddata_x),
               aes(label=label,x=x, y=0, colour=labs$group))+theme_bw()


#NRI
library(nricens)

{
  ROC_data<-data.frame(time=Merged_data$OS,status=Merged_data$Survival_status,riskScore=scoredf$TotalScore,SampleID=Merged_data$SampleID,Death_reason=Merged_data$Death_reason,staging=Merged_data$Stage,anatomical_location=Merged_data$histological_type,sex=Merged_data$Sex,age=Merged_data$Age)
  ROC_data$staging_category<-0
  ROC_data[which(ROC_data$staging%in%c("I","IA","IB")),]$staging_category<-1
  ROC_data[which(ROC_data$staging%in%c("II","IIA","IIB")),]$staging_category<-2
  ROC_data[which(ROC_data$staging%in%c("IIIA","IIIB","IIIC")),]$staging_category<-3
  ROC_data[which(ROC_data$staging%in%c("IV","IVA","IVB")),]$staging_category<-4
  
  #intrahepatic #hilar #extrahepatic
  ROC_data<-ROC_data[which(ROC_data$anatomical_location=="intrahepatic"),]
  ROC_data<-ROC_data[which(!is.na(ROC_data$time) & ROC_data$time!="Unknown" & ROC_data$SampleID%notin% ROC_data[which(ROC_data$Death_reason==" 0"),]$SampleID),]
}
ROC_data$riskScore_category<-2
ROC_data[which(ROC_data$riskScore>as.numeric(quantile(ROC_data$riskScore)[4])),]$riskScore_category<-3
ROC_data[which(ROC_data$riskScore<as.numeric(quantile(ROC_data$riskScore)[2])),]$riskScore_category<-1
ROC_data$age<-as.numeric(ROC_data$age)
time= as.numeric(ROC_data$time)
event= ROC_data$status
z.std= as.matrix(subset(ROC_data, select = c(age,staging_category)))
#z.new=as.matrix(subset(ROC_data, select = c(age,riskScore)))
z.new= as.matrix(subset(ROC_data, select = c(age,staging_category,riskScore)))
mstd= coxph(Surv(time,event) ~ ., data.frame(time,event,z.std), x=TRUE)
mnew= coxph(Surv(time,event) ~ ., data.frame(time,event,z.new), x=TRUE)
p.std= get.risk.coxph(mstd, t0=500)
p.new= get.risk.coxph(mnew, t0=500)

##KM-estimate 
A<-nricens(mdl.std= mstd, mdl.new = mnew, t0 = 500, cut = 0.05,
           niter = 1000, updown = 'diff',alpha = 0.05)
A
#dCCA
#新模型较旧模型重分类正确的比例提高24.83%，
#换句话说增加一个预测变量的新模型预测的准确度比旧模型更好。

B<-nricens(mdl.std= mstd, mdl.new = mnew, t0 = 500, cut = 0.05,
           niter = 1000, updown = 'diff',alpha = 0.05)
B
#pCCA
#新模型较旧模型重分类正确的比例提高19.57%，
#换句话说增加一个预测变量的新模型预测的准确度比旧模型更好。

C<-nricens(mdl.std= mstd, mdl.new = mnew, t0 = 500, cut = 0.05,
           niter = 1000, updown = 'diff',alpha = 0.05)
C
#iCCA
#新模型较旧模型重分类正确的比例提高69.01%，
#换句话说增加一个预测变量的新模型预测的准确度比旧模型更好。


#####桑基图###
setwd("/GeneCloud003/genecloud/Org_terminal/org_139/terminal/zuolj_18846149734/2022_12_28_CHOL_FZQ/0515_result/01_妗戝熀鍥綺438sample")
.libPaths("/GeneCloud003/genecloud/Org_terminal/org_139/terminal/zuolj_18846149734/Package/R")

if (!require("remotes")) install.packages("remotes")
remotes::install_github("wch/webshot")
devtools::install_github("fbreitwieser/sankeyD3")
if(!is_phantomjs_installed()){install_phantomjs()}
webshot::install_phantomjs()

# 
library(sankeyD3)
library(webshot)


# ==== Nodes ====
nodes = data.frame("name" = c("All Samples",
                              "Tissue contamination <= 25 %", "Tissue contamination > 25 %",
                              "FDR >= 0.1",  "FDR < 0.1",
                              "Purified Cohort","Verification Cohort"),
                   "group" = c("C1","C2","C3","C4","C5","C6","C7"))

nodes = data.frame("name" = c("All Samples",
                              "AIS", "MIA", "IAC",
                              "HAS-group","LAS-group"),
                   "group" = c("C1","C2","C3","C4","C5","C6"))
links = as.data.frame(matrix(c(
  0, 1, 3, 1,9,
  0, 2, 3, 2,9,
  0, 3, 3, 3,9,
  1, 4, 1, 4,1,
  1, 5, 2, 5,1,
  2, 4, 1, 4,2,
  2, 5, 2, 5,2,
  3, 4, 3, 4,3),
  byrow = TRUE, ncol = 5))

# ==== Links ====
links = as.data.frame(matrix(c(
  0, 1, 301, 1,9,
  0, 2, 137, 2,9,
  1, 3, 164, 3,1,
  1, 4, 137, 4,1,
  2, 3, 10,  3,2,
  2, 4, 127, 4,2,
  3, 5, 164, 5,3,
  3, 6, 10,  6,3,
  4, 6, 264, 6,4),
  byrow = TRUE, ncol = 5))
names(links) = c("source", "target", "value","group","group2")
my_color <- 'd3.scaleOrdinal() 
.domain(["C1","C2","C3","C4","C5","C6","C7","1","2","3","4","5","6"]) 
.range(["#8dd3c7","#33a02c","#b2df8a", "#1f78b4","#a6cee3","#fb9a99","#ff7f00",
"#33a02c","#b2df8a", "#1f78b4","#a6cee3","#fb9a99","#ff7f00"])'


## ==== Plot ====
sn <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "source", Target = "target",Value = "value", NodeID = "name",
                    fontSize= 12, nodeWidth = 30, orderByPath = TRUE,
                    nodePadding = 15,
                    NodeGroup = "group",LinkGroup = "group",
                    nodeShadow = TRUE,
                    colourScale = my_color)

sn
saveNetwork(sn, "sn_pink_scRNA.html")

webshot("sn_pink_scRNA.html" , "sn_pink_scRNA.pdf",
        #vwidth = 992,
        #vheight = 300,
        debug = TRUE)


