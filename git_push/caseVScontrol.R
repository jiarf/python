

setwd("/home/devdata/jiarongf/nyy/tabfile")

# 1.数据合并 ------------------------------------------------------------------

# case_133076_P1  case_192005_P1  case_192010_P2  case_192020_P1  control_060312_P1  control_192008_P2  control_192013_P1  control_192022_P1
rm(list=ls())
case_133076_P1_gene_abund<- read.delim("/home/devdata/jiarongf/nyy/tabfile/case_133076_P1_gene_abund.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)#60676
###colnames(P4_case_gene_abund)
gene_matrix <- case_133076_P1_gene_abund[,c(1,2,9)]
colnames(gene_matrix) <- c('Gene.ID','Gene.Name','case_133076_P1')

res = c( 'case_192005_P1' , 'case_192010_P2',  'case_192020_P1' , 'control_060312_P1' , 'control_192008_P2',  'control_192013_P1'  ,'control_192022_P1')
for (i in res) {
  data <- read.delim(file = paste("/home/devdata/jiarongf/nyy/tabfile/",i,"_gene_abund.tab",sep=''))
  data <- data[,c(1,9)]
  colnames(data) <- c('Gene.ID',i)
  gene_matrix <- merge(x=gene_matrix,y=data,by="Gene.ID",all =TRUE) ##all =TRUE
}

gene_matrix<-gene_matrix[gene_matrix$Gene.Name !='-',]#name
# gene_matrix<-gene_matrix[gene_matrix$Gene.Name !='.',]
exp_data <- aggregate(gene_matrix[,c(-1,-2)],by = list(gene_matrix$Gene.Name), FUN = mean)
rownames(exp_data)<-exp_data[,1]
exp_data<-exp_data[,-1]
exp_data<-exp_data[rowSums(exp_data)>0,]
write.table(exp_data, file = "/home/devdata/jiarongf/nyy/tabfile/gene_matrix.tab", quote = FALSE,sep="\t",row.names = TRUE)

# 2.可视化样本间的相似性分析 ---------------------------------------------------------------
rm(list=ls())

gene_TPM_matrix <- read.delim("/home/devdata/jiarongf/nyy/tabfile/gene_matrix.tab",stringsAsFactors = F)


#cor_psy <- cor(gene_TPM_matrix)
# gene_TPM_matrix <- as.numeric(unlist(gene_TPM_matrix))
#cortest_psy <- cor.test(gene_TPM_matrix,gene_TPM_matrix,method = "pearson")
# gene_TPM_matrix<-as.matrix(gene_TPM_matrix)
data <- round(cor(gene_TPM_matrix),6)
# cor
library(corrplot)
pdf(file = "./cor.pdf",width = 10,height = 10)
corrplot(data,method = "circle",order = "original",type = "upper",tl.pos = "lt",mar = c(0,0,0,0),tl.srt = 45,
         tl.col = "blue",tl.cex = 1.2,col = colorRampPalette(c("blue","white","red"))(10),cl.pos = "r", is.corr = FALSE)
corrplot(data,number.digits = 5,number.cex = 1.2,add = TRUE,method = "number",
         type = "lower",order = "AOE",diag = FALSE,tl.pos = "n",cl.pos = "n")
dev.off()
# 3.PCA ----------------------------------------------------------------
###2.1 PCA.R
rm(list=ls())
library(ropls)
source("PCA.R")
exp_data <- read.delim("/home/devdata/jiarongf/nyy/tabfile/gene_matrix.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
sampledata<-read.delim("/home/devdata/jiarongf/nyy/tabfile/metedata.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)
sample1 <- as.data.frame(sampledata[,3])
sample1$name <- colnames(exp_data)

colnames(sample1) <- c("Group","SampleID")
color=list(case = "#fc8d59",control = "#99d594")
# dir.create("/home/devdata/jiarongf/nyy/tabfile/result/")
PCA_analysis(exp_data,pData=sample1,color=color, "", group_name='Group',filepath="./result/")



###2.2 PCAtools
rm(list=ls())
library(PCAtools)
exp_data <- read.delim("/home/devdata/jiarongf/nyy/tabfile/gene_matrix.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
metadata <- read.delim("/home/devdata/jiarongf/nyy/tabfile/metedata.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)
rownames(metadata) = colnames(exp_data)
p <- pca(exp_data, metadata = metadata)
pdf('./result/PCA_PCAtool.pdf', width = 10, height = 8)
biplot(p,colby = 'group',title = 'PCA',pointSize = 8,labSize = 6.0,labhjust = 1.0,labvjust = 1.5,legendPosition = 'right')
dev.off()



# 4.差异 --------------------------------------------------------------------
# p<0.05
----------------------------------------------------------------------
  rm(list = ls())
# t.test
edata <- read.delim("/home/devdata/jiarongf/nyy/tabfile/gene_matrix.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
pvalue<-c();tstat<-c();meansControl<-c();meansCase<-c();FC<-c();log2FC<-c()
Control_list<-c('control_060312_P1','control_192008_P2','control_192013_P1','control_192022_P1')
Case_list<-c('case_133076_P1','case_192005_P1','case_192010_P2','case_192020_P1')
edata<-edata+0.01


for (i in 1:nrow(edata)){
  #t.test
  Control_value<-edata[i,Control_list]
  Case_value<-edata[i,Case_list]
  result<-t.test(as.numeric(Case_value), as.numeric(Control_value), paired=FALSE);
  pvalue[i]<-result[["p.value"]]
  tstat[i]<-result[["statistic"]][["t"]]
  meansControl[i]<-mean(as.numeric(Control_value))
  meansCase[i]<-mean(as.numeric(Case_value))
  FC[i]<-mean(as.numeric(Case_value))/mean(as.numeric(Control_value))
  log2FC[i]<-log2(FC[i])
}
p_bf = p.adjust(pvalue, method = "bonferroni")
dif_data = data.frame(rownames(edata),pvalue,tstat,meansControl,meansCase,FC,log2FC,p_bf,stringsAsFactors=FALSE)

colnames(dif_data)[1] = 'SYMBOL'
rownames(dif_data) = dif_data$SYMBOL

#将表达矩阵同t.test结果合并
edata$SYMBOL<-rownames(edata)
dif_data_all = merge(edata,dif_data,by="SYMBOL",all =TRUE) ##all =TRUE
logFC_cutoff = log2(1)
dif_data_all$change =factor(ifelse(dif_data_all$pvalue < 0.05 & abs(dif_data_all$log2FC) >= logFC_cutoff,
                                   ifelse(dif_data_all$log2FC >= logFC_cutoff , 'UP', 'DOWN' ), 'NOT'),levels=c('UP', 'DOWN', 'NOT'))
write.table(dif_data_all, file = "./dif_data.tab", quote = FALSE,sep="\t",row.names = FALSE) ##总的差异的基因

dif_up_data <- dif_data_all[dif_data_all$change == "UP",]######up基因
dif_down_data <- dif_data_all[dif_data_all$change == "DOWN",]######DOWN基因
write.table(dif_data_all, file = "./result/dif_data_all.csv",sep = ",",row.names=FALSE)#######输出总的差异基因
write.table(dif_up_data, file = "./result/dif_up_data.csv",sep = ",",row.names=FALSE)#######输出上调的差异基因
write.table(dif_down_data, file = "./result/dif_down_data.csv",sep = ",",row.names=FALSE)######输出下调的差异基因
#####################################the first part
edata <- read.delim("/home/devdata/jiarongf/nyy/tabfile/data/dif_data.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
first_part_up_data <- edata[edata$change == "UP",]
first_part_down_data <- edata[edata$change == "DOWN",]
write.table(first_part_up_data, file = "./data/first_part_up_data.csv",sep = ",",row.names=FALSE)#######输出上调的差异基因
write.table(first_part_down_data, file = "./data/first_part_down_data.csv",sep = ",",row.names=FALSE)#######输出下调的差异基因
# 4.1volcano --------------------------------------------------------------
#  volcanomap_up down second
rm(list=ls())
library("ggplot2")
library("ggrepel")
nrDEG = read.table('./dif_data.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
# nrDEG$change_V1 = factor(nrDEG$change,levels = c('UP', 'DOWN','NOT'))
nrDEG = nrDEG[nrDEG$meansControl > 0,]
nrDEG = nrDEG[nrDEG$meansCase > 0,]

this_tile <- paste0( 'The number of up gene is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),"\n",
                     'The number of down gene is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ))

upgene = nrDEG[nrDEG$change == 'UP',];upgene <- upgene[order(upgene$pvalue),]
rownames(upgene) = upgene$SYMBOL;top5_up = rownames(upgene)[1:5]

downgene = nrDEG[nrDEG$change == 'DOWN',];downgene <- downgene[order(downgene$pvalue),]
rownames(downgene) = downgene$SYMBOL;top5_down = rownames(downgene)[1:5]

volcano = 
  ggplot(data = nrDEG, aes( x = log2FC, y = -log10(pvalue), color = change)) +
  scale_x_continuous(limits = c(-4,4),breaks = seq(-4,4,1))+
  scale_y_continuous(limits = c(0,5),breaks = seq(0,6,1))+
  geom_point( alpha = 0.4, size = 1.75) +  
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+  #y轴分界线
  # geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+ #x轴分界线
  theme_bw() + theme(legend.position = "none") +
  xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
  ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
  scale_colour_manual( values = c('blue','black','red')) +
  geom_text_repel(data=subset(nrDEG, SYMBOL %in% top5_up | SYMBOL %in% top5_down ), aes(label=SYMBOL),col="black",alpha = 1)
# geom_label_repel(data=subset(nrDEG, change_V2 == 'UP' | change_V2 == 'DOWN'),aes(label=SYMBOL,fill=factor(change_V2)), col="black",alpha = 1)
print(volcano)
plotfile=paste('./result/','CASEvsCONTROL_second','_Volcano.png',sep='')
plotfile=paste('./result/','CASEvsCONTROL_second','_Volcano.pdf',sep='')
ggsave(plotfile, plot=volcano, dpi = 600,width = 10, height = 8)
####volcano the first
rm(list=ls())
library("ggplot2")
library("ggrepel")
nrDEG = read.table('/home/devdata/jiarongf/nyy/tabfile/data/dif_data.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)

nrDEG = nrDEG[nrDEG$meansControl > 0,]
nrDEG = nrDEG[nrDEG$meansCase > 0,]

this_tile <- paste0( 'The number of up gene is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),"\n",
                     'The number of down gene is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ))

upgene = nrDEG[nrDEG$change == 'UP',];upgene <- upgene[order(upgene$pvalue),]
rownames(upgene) = upgene$SYMBOL;top5_up = rownames(upgene)[1:5]

downgene = nrDEG[nrDEG$change == 'DOWN',];downgene <- downgene[order(downgene$pvalue),]
rownames(downgene) = downgene$SYMBOL;top5_down = rownames(downgene)[1:5]

volcano =
  ggplot(data = nrDEG, aes( x = log2FC, y = -log10(pvalue), color = change)) +
  scale_x_continuous(limits = c(-5,4),breaks = seq(-4,4,1))+
  scale_y_continuous(limits = c(0,5),breaks = seq(0,6,1))+
  geom_point( alpha = 0.4, size = 1.75) +
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+  #y轴分界线
  # geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+ #x轴分界线
  theme_bw() + theme(legend.position = "none") +
  xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
  ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
  scale_colour_manual( values = c('blue','black','red')) +
  geom_text_repel(data=subset(nrDEG, SYMBOL %in% top5_up | SYMBOL %in% top5_down ), aes(label=SYMBOL),col="black",alpha = 1)
# geom_label_repel(data=subset(nrDEG, change_V2 == 'UP' | change_V2 == 'DOWN'),aes(label=SYMBOL,fill=factor(change_V2)), col="black",alpha = 1)
print(volcano)
# plotfile=paste('./result/','CASEvsCONTROL','_Volcano.png',sep='')
plotfile=paste('./result/','CASEvsCONTROL','_Volcano.pdf',sep='')
ggsave(plotfile, plot=volcano, dpi = 600,width = 10, height = 8)

# ### 4差异表达分析-火山图(有图例的) ----------------------------------------------------
# library(ggplot2)
# library(ggrepel)
# rm(list=ls())
# dif_data <- read.delim(file = './data/dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
# dif_data<-dif_data[order(dif_data$pvalue),]
# topgene<-c(dif_data$SYMBOL[dif_data$change=='UP'][1:5],dif_data$SYMBOL[dif_data$change=='DOWN'][1:5])
# sig<-rep('',nrow(dif_data))
# for (i in 1:nrow(dif_data)){
#   if (dif_data$SYMBOL[i] %in% topgene){
#     sig[i]<-'most_dif'
#   } else if (dif_data$pvalue[i]<0.01){
#     sig[i]<-'True'
#   } else{
#     sig[i]<-'False'
#   }
# }
# 
# dif_data$significant <- as.factor(sig)
# dif_data$change <- factor(dif_data$change,levels=c('UP', 'DOWN', 'NOT'))
# logFC_cutoff = log2(1)
# dif_data$change1 =factor(ifelse(dif_data$pvalue < 0.01 & abs(dif_data$log2FC) >= logFC_cutoff,
#                                 ifelse(dif_data$log2FC >= logFC_cutoff , 'UP', 'DOWN' ), 'NOT'),levels=c('UP', 'DOWN', 'NOT'))
# this_tile <- paste0(nrow(dif_data[dif_data$change =='UP', ] ),' Up-regulated genes', "; ",
#                     nrow(dif_data[dif_data$change =='DOWN', ]), ' Down-regulated genes')
# 
# p=ggplot(data=dif_data, aes(x=log2FC, y =-log10(pvalue),color=change))+
#   geom_point(alpha=0.8, size=1.2)+
#   scale_x_continuous(limits = c(-9,9),breaks = seq(-9,9,1))+
#   scale_color_manual(values =c('red',"blue","#595959"))+
#   labs(title="", x="log2 (fold change)",y="-log10 (p-value)")+
#   ggtitle( this_tile )+
#   theme( plot.title = element_text(size =15, hjust = 0.5)) +
#   geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
#   geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
#   theme_set( theme_set(theme_bw( base_size = 12)))+
#   geom_text_repel(data=subset(dif_data, significant =='most_dif'), aes(label=SYMBOL),col="black",size=2)
# plotfile=paste('./result/','dif_data','_Volcano.png',sep='')
# ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)
# dif_data_most<-dif_data[dif_data$significant=='most_dif',]
# write.table(dif_data_most, file = "./data/dif_data_mosttop20.tab", quote = FALSE,sep="\t",row.names = FALSE)


# 5.kegg_GO(-cancer) ---------------------------------------------------------------
rm(list=ls())
library( "clusterProfiler")
load(file = '/home/devdata/zhuzn/enrichment_annotation/crab_eating_enrichment.RData')
dif_data <- read.delim(file = './dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
gene_up = dif_data[dif_data$change=='UP','SYMBOL']#175
gene_down = dif_data[dif_data$change=='DOWN','SYMBOL']#145
# gene_up = dif_data[dif_data$pvalue<0.01 & dif_data$FC>1, 'SYMBOL']#255
# gene_down = dif_data[dif_data$pvalue<0.01 & dif_data$FC<1, 'SYMBOL']#233
# dif_gene<-c(gene_up,gene_down)



###GOBP
term2gene <- gobp[,c('GO','SYMBOL')]
term2gene<-term2gene[!is.na(term2gene$SYMBOL),]
term2name <- gobp[,c('GO','Name')]

ego_BP <- enricher(
  gene = gene_up,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_UP_result<-as.data.frame(ego_BP@result)



ego_BP <- enricher(
  gene = gene_down,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_DOWN_result<-as.data.frame(ego_BP@result)
# write.table(ego_BP_DOWN_result, file = "./result/ego_BP_DOWN_result.tab", quote = FALSE,sep="\t",row.names = FALSE)
# write.table(ego_BP_UP_result, file = "./result/ego_BP_UP_result.tab", quote = FALSE,sep="\t",row.names = FALSE)


###BP_RESULT
BP_UP_top20 <- ego_BP_UP_result[order(ego_BP_UP_result$pvalue),][1:20,]
BP_DOWN_top20 <- ego_BP_DOWN_result[order(ego_BP_DOWN_result$pvalue),][1:20,]
overlap <- intersect(BP_UP_top20$Description,BP_DOWN_top20$Description)
path <- c(BP_UP_top20$Description,BP_DOWN_top20$Description)
BP_UP <- ego_BP_UP_result[ego_BP_UP_result$Description %in% path,]
BP_UP <- BP_UP[,c(2,5,9)]
BP_UP$change <- rep("UP")
BP_DOWN <- ego_BP_DOWN_result[ego_BP_DOWN_result$Description %in% path,]
BP_DOWN <- BP_DOWN[,c(2,5,9)]
BP_DOWN$change <- rep("DOWN")
BP_path <- rbind(BP_UP,BP_DOWN)


## MF
term2gene <- gomf[,c(1,5)]
term2name <- gomf[,c(1,2)]
### UP
ego_MF_UP <- enricher(gene = gene_up,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
ego_MF_UP_result<-as.data.frame(ego_MF_UP@result) # 148


### DOWN
ego_MF_DOWN <- enricher(gene = gene_down,
                        pvalueCutoff = 0.05,pAdjustMethod = "BH",
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                        TERM2GENE = term2gene,TERM2NAME = term2name)
ego_MF_DOWN_result<-as.data.frame(ego_MF_DOWN@result) # 118

###MF_RESULT
MF_UP_top20 <- ego_MF_UP_result[order(ego_MF_UP_result$pvalue),][1:20,]
MF_DOWN_top20 <- ego_MF_DOWN_result[order(ego_MF_DOWN_result$pvalue),][1:20,]
overlap <- intersect(MF_UP_top20$Description,MF_DOWN_top20$Description)
path <- c(MF_UP_top20$Description,MF_DOWN_top20$Description)
MF_UP <- ego_MF_UP_result[ego_MF_UP_result$Description %in% path,]
MF_UP <- MF_UP[,c(2,5,9)]
MF_UP$change <- rep("UP")
MF_DOWN <- ego_MF_DOWN_result[ego_MF_DOWN_result$Description %in% path,]
MF_DOWN <- MF_DOWN[,c(2,5,9)]
MF_DOWN$change <- rep("DOWN")
MF_path <- rbind(MF_UP,MF_DOWN)

## CC
term2gene <- gocc[,c(1,5)]
term2name <- gocc[,c(1,2)]
### UP
ego_CC_UP <- enricher(gene = gene_up,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
ego_CC_UP_result<-as.data.frame(ego_CC_UP@result) # 95
### DOWN
ego_CC_DOWN <- enricher(gene = gene_down,
                        pvalueCutoff = 0.05,pAdjustMethod = "BH",
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                        TERM2GENE = term2gene,TERM2NAME = term2name)
ego_CC_DOWN_result<-as.data.frame(ego_CC_DOWN@result) # 100
###CC_RESULT
CC_UP_top20 <- ego_CC_UP_result[order(ego_CC_UP_result$pvalue),][1:20,]
CC_DOWN_top20 <- ego_CC_DOWN_result[order(ego_CC_DOWN_result$pvalue),][1:20,]
overlap <- intersect(CC_UP_top20$Description,CC_DOWN_top20$Description)
path <- c(CC_UP_top20$Description,CC_DOWN_top20$Description)
CC_UP <- ego_CC_UP_result[ego_CC_UP_result$Description %in% path,]
CC_UP <- CC_UP[,c(2,5,9)]
CC_UP$change <- rep("UP")
CC_DOWN <- ego_CC_DOWN_result[ego_CC_DOWN_result$Description %in% path,]
CC_DOWN <- CC_DOWN[,c(2,5,9)]
CC_DOWN$change <- rep("DOWN")
CC_path <- rbind(CC_UP,CC_DOWN)


grep("cancer",gocc$Name)
## KEGG
term2gene <- kegg[,c(1,5)]
term2name <- kegg[,c(1,2)]
###UP
ekegg_UP <- enricher(gene = gene_up,
                     pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                     TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_UP<-as.data.frame(ekegg_UP@result) # 148
###DOWN
ekegg_DOWN <- enricher(gene = gene_down,
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",
                       minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                       TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_DOWN<-as.data.frame(ekegg_DOWN@result) # 118
###kegg_RESULT
ekegg_UP_top20 <- ekegg_UP[order(ekegg_UP$pvalue),][1:20,]
ekegg_DOWN_top20 <- ekegg_DOWN[order(ekegg_DOWN$pvalue),][1:20,]
overlap <- intersect(ekegg_UP_top20$Description,ekegg_DOWN_top20$Description)
path <- c(ekegg_UP_top20$Description,ekegg_DOWN_top20$Description)
ekegg_UP <- ekegg_UP[ekegg_UP$Description %in% path,]
ekegg_UP <- ekegg_UP[,c(2,5,9)]
ekegg_UP$change <- rep("UP")
ekegg_DOWN <- ekegg_DOWN[ekegg_DOWN$Description %in% path,]
ekegg_DOWN <- ekegg_DOWN[,c(2,5,9)]
ekegg_DOWN$change <- rep("DOWN")
ekegg_path <- rbind(ekegg_UP,ekegg_DOWN)

write.table(ekegg_path, file = "./enrichment/CASEvsCONTROL_kegg_enrich_dif.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(BP_path, file = "./enrichment/CASEvsCONTRO_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(MF_path, file = "./enrichment/CASEvsCONTRO_MF_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(CC_path, file = "./enrichment/CASEvsCONTRO_CC_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)
##############kegg富集时削去包含cancer的pathway
class(kegg)
path<-as.data.frame(c(kegg[grep("cancer",kegg$Name),]))#2171
no_cancer_kegg_tab_result <- ekegg_path[!(ekegg_path$Description %in% path$Name),]
write.table(no_cancer_kegg_tab_result, file = "./enrichment/CASEvsCONTROL_kegg_enrich_no_cancer_kegg_tab_result.tab", quote = FALSE,sep="\t",row.names = FALSE)


#气泡图
rm(list=ls())
#install.packages("ggnewscale")
library(ggnewscale)
res = c('CASEvsCONTROL_kegg_enrich_no_cancer_kegg_tab_result','CASEvsCONTRO_BP_enrich','CASEvsCONTRO_MF_enrich','CASEvsCONTRO_CC_enrich')
i = 'CASEvsCONTROL_kegg_enrich_dif'   #limits:10/12
i = 'CASEvsCONTRO_BP_enrich'         #limits:10/6
i = 'CASEvsCONTROL_MF_enrich'         #limits:15/3
i = 'CASEvsCONTROL_CC_enrich'         #limits:32/8
########直线图
for (i in res) {
  path<-read.delim(file = paste("./enrichment/",i,'.tab',sep=''),header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  p <- ggplot(path,aes(x = change,y=Description,size=Count,ylab=''))+scale_x_discrete(limits=c("UP","DOWN"))+
    scale_y_discrete(limits=rev(unique(path$Description)))+
    geom_point(data = subset(path,change=='UP'),
               mapping = aes(change,Description,colour =-log10(pvalue)),
    )+
    scale_colour_gradient(low="#FFB6C1",high="red",limits = c(0,32)) +
    new_scale_color() +
    geom_point(data = subset(path,change=='DOWN'),
               mapping = aes(change,Description,colour =-log10(pvalue)),
    )+
    scale_colour_gradient(low="#ADD8E6",high="blue",limits = c(0,8)) +
   
    scale_size_area(name="genecounts")+theme_bw()+
    scale_size_continuous(range=c(1,8))+
    labs(y='',x='',title = i)+
    geom_vline(xintercept = 0.05,linetype =2,colour = 'black')+
    theme(axis.text=element_text(size=12),
         
          axis.title=element_text(size=18),
          panel.grid.major=element_line(colour=NA),
          plot.title = element_text(hjust = 0.5,size = 30),
          legend.text=element_text(size=12),legend.title = element_text(size=12))
  plotfile = paste('./enrichment/',i,'_res.pdf')
  ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 10)
}

# ###5.1 第一批数据进行kegg和gobp分析(1)有cancer --------------------------------------------------
#  volcanomap_up down
rm(list=ls())
library("ggplot2")
library("ggrepel")
nrDEG = read.delim('/home/devdata/jiarongf/nyy/tabfile/the first data/dif_data.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
# nrDEG$change_V1 = factor(nrDEG$change,levels = c('UP', 'DOWN','NOT'))
nrDEG = nrDEG[nrDEG$meansControl > 0,]
nrDEG = nrDEG[nrDEG$meansCase > 0,]

this_tile <- paste0( 'The number of up gene is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),"\n",
                     'The number of down gene is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ))

upgene = nrDEG[nrDEG$change == 'UP',];upgene <- upgene[order(upgene$pvalue),]
rownames(upgene) = upgene$SYMBOL;top5_up = rownames(upgene)[1:5]

downgene = nrDEG[nrDEG$change == 'DOWN',];downgene <- downgene[order(downgene$pvalue),]
rownames(downgene) = downgene$SYMBOL;top5_down = rownames(downgene)[1:5]

volcano = 
  ggplot(data = nrDEG, aes( x = log2FC, y = -log10(pvalue), color = change)) +
  scale_x_continuous(limits = c(-4,4),breaks = seq(-4,4,1))+
  scale_y_continuous(limits = c(0,5),breaks = seq(0,6,1))+
  geom_point( alpha = 0.4, size = 1.75) +  
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+  #y轴分界线
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+ #x轴分界线
  theme_bw() + theme(legend.position = "none") +
  xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
  ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
  scale_colour_manual( values = c('red','blue','black')) +
  geom_text_repel(data=subset(nrDEG, SYMBOL %in% top5_up | SYMBOL %in% top5_down ), aes(label=SYMBOL),col="black",alpha = 1)
# geom_label_repel(data=subset(nrDEG, change_V2 == 'UP' | change_V2 == 'DOWN'),aes(label=SYMBOL,fill=factor(change_V2)), col="black",alpha = 1)
print(volcano)
plotfile=paste('./result/','CASEvsCONTROL','_Volcano.png',sep='')
plotfile=paste('./result/','CASEvsCONTROL','_Volcano.pdf',sep='')
ggsave(plotfile, plot=volcano, dpi = 600,width = 10, height = 8)


#kegg and go
rm(list=ls())
library( "clusterProfiler")
load(file = '/home/devdata/zhuzn/enrichment_annotation/crab_eating_enrichment.RData')
dif_data <- read.delim(file = '/home/devdata/jiarongf/nyy/tabfile/the first data/dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
gene_up = dif_data[dif_data$change=='UP','SYMBOL']#374
gene_down = dif_data[dif_data$change=='DOWN','SYMBOL']#131
# gene_up = dif_data[dif_data$pvalue<0.01 & dif_data$FC>1, 'SYMBOL']#255
# gene_down = dif_data[dif_data$pvalue<0.01 & dif_data$FC<1, 'SYMBOL']#233
# dif_gene<-c(gene_up,gene_down)



###GOBP
term2gene <- gobp[,c('GO','SYMBOL')]
term2gene<-term2gene[!is.na(term2gene$SYMBOL),]
term2name <- gobp[,c('GO','Name')]

ego_BP <- enricher(
  gene = gene_up,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_UP_result<-as.data.frame(ego_BP@result)



ego_BP <- enricher(
  gene = gene_down,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_DOWN_result<-as.data.frame(ego_BP@result)
# write.table(ego_BP_DOWN_result, file = "./result/ego_BP_DOWN_result_up_down(1).tab", quote = FALSE,sep="\t",row.names = FALSE)
# write.table(ego_BP_UP_result, file = "./result/ego_BP_UP_result_up_down(1).tab", quote = FALSE,sep="\t",row.names = FALSE)


###BP_RESULT
BP_UP_top20 <- ego_BP_UP_result[order(ego_BP_UP_result$pvalue),][1:20,]
BP_DOWN_top20 <- ego_BP_DOWN_result[order(ego_BP_DOWN_result$pvalue),][1:20,]
overlap <- intersect(BP_UP_top20$Description,BP_DOWN_top20$Description)
path <- c(BP_UP_top20$Description,BP_DOWN_top20$Description)
BP_UP <- ego_BP_UP_result[ego_BP_UP_result$Description %in% path,]
BP_UP <- BP_UP[,c(2,5,9)]
BP_UP$change <- rep("UP")
BP_DOWN <- ego_BP_DOWN_result[ego_BP_DOWN_result$Description %in% path,]
BP_DOWN <- BP_DOWN[,c(2,5,9)]
BP_DOWN$change <- rep("DOWN")
BP_path <- rbind(BP_UP,BP_DOWN)



## MF
term2gene <- gomf[,c(1,5)]
term2name <- gomf[,c(1,2)]
### UP
ego_MF_UP <- enricher(gene = gene_up,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
ego_MF_UP_result<-as.data.frame(ego_MF_UP@result) # 148


### DOWN
ego_MF_DOWN <- enricher(gene = gene_down,
                        pvalueCutoff = 0.05,pAdjustMethod = "BH",
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                        TERM2GENE = term2gene,TERM2NAME = term2name)
ego_MF_DOWN_result<-as.data.frame(ego_MF_DOWN@result) # 118

###MF_RESULT
MF_UP_top20 <- ego_MF_UP_result[order(ego_MF_UP_result$pvalue),][1:20,]
MF_DOWN_top20 <- ego_MF_DOWN_result[order(ego_MF_DOWN_result$pvalue),][1:20,]
overlap <- intersect(MF_UP_top20$Description,MF_DOWN_top20$Description)
path <- c(MF_UP_top20$Description,MF_DOWN_top20$Description)
MF_UP <- ego_MF_UP_result[ego_MF_UP_result$Description %in% path,]
MF_UP <- MF_UP[,c(2,5,9)]
MF_UP$change <- rep("UP")
MF_DOWN <- ego_MF_DOWN_result[ego_MF_DOWN_result$Description %in% path,]
MF_DOWN <- MF_DOWN[,c(2,5,9)]
MF_DOWN$change <- rep("DOWN")
MF_path <- rbind(MF_UP,MF_DOWN)

## CC
term2gene <- gocc[,c(1,5)]
term2name <- gocc[,c(1,2)]
### UP
ego_CC_UP <- enricher(gene = gene_up,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
ego_CC_UP_result<-as.data.frame(ego_CC_UP@result) # 95
### DOWN
ego_CC_DOWN <- enricher(gene = gene_down,
                        pvalueCutoff = 0.05,pAdjustMethod = "BH",
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                        TERM2GENE = term2gene,TERM2NAME = term2name)
ego_CC_DOWN_result<-as.data.frame(ego_CC_DOWN@result) # 100
###CC_RESULT
CC_UP_top20 <- ego_CC_UP_result[order(ego_CC_UP_result$pvalue),][1:20,]
CC_DOWN_top20 <- ego_CC_DOWN_result[order(ego_CC_DOWN_result$pvalue),][1:20,]
overlap <- intersect(CC_UP_top20$Description,CC_DOWN_top20$Description)
path <- c(CC_UP_top20$Description,CC_DOWN_top20$Description)
CC_UP <- ego_CC_UP_result[ego_CC_UP_result$Description %in% path,]
CC_UP <- CC_UP[,c(2,5,9)]
CC_UP$change <- rep("UP")
CC_DOWN <- ego_CC_DOWN_result[ego_CC_DOWN_result$Description %in% path,]
CC_DOWN <- CC_DOWN[,c(2,5,9)]
CC_DOWN$change <- rep("DOWN")
CC_path <- rbind(CC_UP,CC_DOWN)

## KEGG
term2gene <- kegg[,c(1,5)]
term2name <- kegg[,c(1,2)]
###UP
ekegg_UP <- enricher(gene = gene_up,
                     pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                     TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_UP<-as.data.frame(ekegg_UP@result) # 148
###DOWN
ekegg_DOWN <- enricher(gene = gene_down,
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",
                       minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                       TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_DOWN<-as.data.frame(ekegg_DOWN@result) # 118
###kegg_RESULT
ekegg_UP_top20 <- ekegg_UP[order(ekegg_UP$pvalue),][1:20,]
ekegg_DOWN_top20 <- ekegg_DOWN[order(ekegg_DOWN$pvalue),][1:20,]
overlap <- intersect(ekegg_UP_top20$Description,ekegg_DOWN_top20$Description)
path <- c(ekegg_UP_top20$Description,ekegg_DOWN_top20$Description)
ekegg_UP <- ekegg_UP[ekegg_UP$Description %in% path,]
ekegg_UP <- ekegg_UP[,c(2,5,9)]
ekegg_UP$change <- rep("UP")
ekegg_DOWN <- ekegg_DOWN[ekegg_DOWN$Description %in% path,]
ekegg_DOWN <- ekegg_DOWN[,c(2,5,9)]
ekegg_DOWN$change <- rep("DOWN")
ekegg_path <- rbind(ekegg_UP,ekegg_DOWN)

write.table(ekegg_path, file = "./the first data/CASEvsCONTROL_kegg_enrich_dif1.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(BP_path, file = "./the first data/CASEvsCONTRO_BP_enrich1.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(MF_path, file = "./the first data/CASEvsCONTRO_MF_enrich1.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(CC_path, file = "./the first data/CASEvsCONTRO_CC_enrich1.tab", quote = FALSE,sep="\t",row.names = FALSE)



#气泡图
rm(list=ls())
#install.packages("ggnewscale")
library(ggnewscale)
res = c('CASEvsCONTROL_kegg_enrich_dif1','CASEvsCONTRO_BP_enrich1','CASEvsCONTRO_MF_enrich1','CASEvsCONTRO_CC_enrich1')
i = 'CASEvsCONTROL_kegg_enrich_dif'   #limits:10/12
i = 'CASEvsCONTRO_BP_enrich'         #limits:10/6
i = 'CASEvsCONTROL_MF_enrich'         #limits:15/3
i = 'CASEvsCONTROL_CC_enrich'         #limits:32/8
########直线图
for (i in res) {
  path<-read.delim(file = paste("./the first data/",i,'.tab',sep=''),header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  p <- ggplot(path,aes(x = change,y=Description,size=Count,ylab=''))+scale_x_discrete(limits=c("UP","DOWN"))+
    scale_y_discrete(limits=rev(unique(path$Description)))+
    geom_point(data = subset(path,change=='UP'),
               mapping = aes(change,Description,colour =-log10(pvalue)),
    )+
    scale_colour_gradient(low="#FFB6C1",high="red",limits = c(0,32)) +
    new_scale_color() +
    geom_point(data = subset(path,change=='DOWN'),
               mapping = aes(change,Description,colour =-log10(pvalue)),
    )+
    scale_colour_gradient(low="#ADD8E6",high="blue",limits = c(0,8)) +
    #scale_colour_discrete(c("blue","red"))+
    #scale_colour_gradient(low="blue",high="red") + 
    #scale_shape_manual(values = c(16,10))+
    scale_size_area(name="genecounts")+theme_bw()+
    scale_size_continuous(range=c(1,8))+
    labs(y='',x='',title = i)+
    geom_vline(xintercept = 0.05,linetype =2,colour = 'black')+
    theme(axis.text=element_text(size=12),
          #axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title=element_text(size=18),
          panel.grid.major=element_line(colour=NA),
          plot.title = element_text(hjust = 0.5,size = 30),
          legend.text=element_text(size=12),legend.title = element_text(size=12))
  plotfile = paste('./the first data/',i,'_res.pdf')
  ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 10)
}



# ###5.2 第一批数据进行kegg和gobp分析(1.5)有cancer---------------------------------------------------------------------
rm(list=ls())
library( "clusterProfiler")
load(file = '/home/devdata/zhuzn/enrichment_annotation/crab_eating_enrichment.RData')
dif_data <- read.delim(file = './up_down(1.5).xls',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
gene_up = dif_data[dif_data$change=='UP','SYMBOL']#374
gene_down = dif_data[dif_data$change=='DOWN','SYMBOL']#131
# gene_up = dif_data[dif_data$pvalue<0.01 & dif_data$FC>1, 'SYMBOL']#255
# gene_down = dif_data[dif_data$pvalue<0.01 & dif_data$FC<1, 'SYMBOL']#233
# dif_gene<-c(gene_up,gene_down)



###GOBP
term2gene <- gobp[,c('GO','SYMBOL')]
term2gene<-term2gene[!is.na(term2gene$SYMBOL),]
term2name <- gobp[,c('GO','Name')]

ego_BP <- enricher(
  gene = gene_up,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_UP_result<-as.data.frame(ego_BP@result)



ego_BP <- enricher(
  gene = gene_down,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_DOWN_result<-as.data.frame(ego_BP@result)
# write.table(ego_BP_DOWN_result, file = "./result/ego_BP_DOWN_result_up_down(1).tab", quote = FALSE,sep="\t",row.names = FALSE)
# write.table(ego_BP_UP_result, file = "./result/ego_BP_UP_result_up_down(1).tab", quote = FALSE,sep="\t",row.names = FALSE)


###BP_RESULT
BP_UP_top20 <- ego_BP_UP_result[order(ego_BP_UP_result$pvalue),][1:20,]
BP_DOWN_top20 <- ego_BP_DOWN_result[order(ego_BP_DOWN_result$pvalue),][1:20,]
overlap <- intersect(BP_UP_top20$Description,BP_DOWN_top20$Description)
path <- c(BP_UP_top20$Description,BP_DOWN_top20$Description)
BP_UP <- ego_BP_UP_result[ego_BP_UP_result$Description %in% path,]
BP_UP <- BP_UP[,c(2,5,9)]
BP_UP$change <- rep("UP")
BP_DOWN <- ego_BP_DOWN_result[ego_BP_DOWN_result$Description %in% path,]
BP_DOWN <- BP_DOWN[,c(2,5,9)]
BP_DOWN$change <- rep("DOWN")
BP_path <- rbind(BP_UP,BP_DOWN)



## MF
term2gene <- gomf[,c(1,5)]
term2name <- gomf[,c(1,2)]
### UP
ego_MF_UP <- enricher(gene = gene_up,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
ego_MF_UP_result<-as.data.frame(ego_MF_UP@result) # 148


### DOWN
ego_MF_DOWN <- enricher(gene = gene_down,
                        pvalueCutoff = 0.05,pAdjustMethod = "BH",
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                        TERM2GENE = term2gene,TERM2NAME = term2name)
ego_MF_DOWN_result<-as.data.frame(ego_MF_DOWN@result) # 118

###MF_RESULT
MF_UP_top20 <- ego_MF_UP_result[order(ego_MF_UP_result$pvalue),][1:20,]
MF_DOWN_top20 <- ego_MF_DOWN_result[order(ego_MF_DOWN_result$pvalue),][1:20,]
overlap <- intersect(MF_UP_top20$Description,MF_DOWN_top20$Description)
path <- c(MF_UP_top20$Description,MF_DOWN_top20$Description)
MF_UP <- ego_MF_UP_result[ego_MF_UP_result$Description %in% path,]
MF_UP <- MF_UP[,c(2,5,9)]
MF_UP$change <- rep("UP")
MF_DOWN <- ego_MF_DOWN_result[ego_MF_DOWN_result$Description %in% path,]
MF_DOWN <- MF_DOWN[,c(2,5,9)]
MF_DOWN$change <- rep("DOWN")
MF_path <- rbind(MF_UP,MF_DOWN)

## CC
term2gene <- gocc[,c(1,5)]
term2name <- gocc[,c(1,2)]
### UP
ego_CC_UP <- enricher(gene = gene_up,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                      TERM2GENE = term2gene,TERM2NAME = term2name)
ego_CC_UP_result<-as.data.frame(ego_CC_UP@result) # 95
### DOWN
ego_CC_DOWN <- enricher(gene = gene_down,
                        pvalueCutoff = 0.05,pAdjustMethod = "BH",
                        minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                        TERM2GENE = term2gene,TERM2NAME = term2name)
ego_CC_DOWN_result<-as.data.frame(ego_CC_DOWN@result) # 100
###CC_RESULT
CC_UP_top20 <- ego_CC_UP_result[order(ego_CC_UP_result$pvalue),][1:20,]
CC_DOWN_top20 <- ego_CC_DOWN_result[order(ego_CC_DOWN_result$pvalue),][1:20,]
overlap <- intersect(CC_UP_top20$Description,CC_DOWN_top20$Description)
path <- c(CC_UP_top20$Description,CC_DOWN_top20$Description)
CC_UP <- ego_CC_UP_result[ego_CC_UP_result$Description %in% path,]
CC_UP <- CC_UP[,c(2,5,9)]
CC_UP$change <- rep("UP")
CC_DOWN <- ego_CC_DOWN_result[ego_CC_DOWN_result$Description %in% path,]
CC_DOWN <- CC_DOWN[,c(2,5,9)]
CC_DOWN$change <- rep("DOWN")
CC_path <- rbind(CC_UP,CC_DOWN)

## KEGG
term2gene <- kegg[,c(1,5)]
term2name <- kegg[,c(1,2)]
###UP
ekegg_UP <- enricher(gene = gene_up,
                     pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                     TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_UP<-as.data.frame(ekegg_UP@result) # 148
###DOWN
ekegg_DOWN <- enricher(gene = gene_down,
                       pvalueCutoff = 0.05,pAdjustMethod = "BH",
                       minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                       TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_DOWN<-as.data.frame(ekegg_DOWN@result) # 118
###kegg_RESULT
ekegg_UP_top20 <- ekegg_UP[order(ekegg_UP$pvalue),][1:20,]
ekegg_DOWN_top20 <- ekegg_DOWN[order(ekegg_DOWN$pvalue),][1:20,]
overlap <- intersect(ekegg_UP_top20$Description,ekegg_DOWN_top20$Description)
path <- c(ekegg_UP_top20$Description,ekegg_DOWN_top20$Description)
ekegg_UP <- ekegg_UP[ekegg_UP$Description %in% path,]
ekegg_UP <- ekegg_UP[,c(2,5,9)]
ekegg_UP$change <- rep("UP")
ekegg_DOWN <- ekegg_DOWN[ekegg_DOWN$Description %in% path,]
ekegg_DOWN <- ekegg_DOWN[,c(2,5,9)]
ekegg_DOWN$change <- rep("DOWN")
ekegg_path <- rbind(ekegg_UP,ekegg_DOWN)

write.table(ekegg_path, file = "./enrichment/CASEvsCONTROL_kegg_enrich_dif(1.5).tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(BP_path, file = "./enrichment/CASEvsCONTRO_BP_enrich(1.5).tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(MF_path, file = "./enrichment/CASEvsCONTRO_MF_enrich(1.5).tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(CC_path, file = "./enrichment/CASEvsCONTRO_CC_enrich(1.5).tab", quote = FALSE,sep="\t",row.names = FALSE)



#气泡图
rm(list=ls())
#install.packages("ggnewscale")
library(ggnewscale)
res = c('CASEvsCONTROL_kegg_enrich_dif(1.5)','CASEvsCONTRO_BP_enrich(1.5)','CASEvsCONTRO_MF_enrich(1.5)','CASEvsCONTRO_CC_enrich(1.5)')
i = 'CASEvsCONTROL_kegg_enrich_dif'   #limits:10/12
i = 'CASEvsCONTRO_BP_enrich'         #limits:10/6
i = 'CASEvsCONTROL_MF_enrich'         #limits:15/3
i = 'CASEvsCONTROL_CC_enrich'         #limits:32/8
########直线图
for (i in res) {
  path<-read.delim(file = paste("./enrichment/",i,'.tab',sep=''),header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  p <- ggplot(path,aes(x = change,y=Description,size=Count,ylab=''))+scale_x_discrete(limits=c("UP","DOWN"))+
    scale_y_discrete(limits=rev(unique(path$Description)))+
    geom_point(data = subset(path,change=='UP'),
               mapping = aes(change,Description,colour =-log10(pvalue)),
    )+
    scale_colour_gradient(low="#FFB6C1",high="red",limits = c(0,32)) +
    new_scale_color() +
    geom_point(data = subset(path,change=='DOWN'),
               mapping = aes(change,Description,colour =-log10(pvalue)),
    )+
    scale_colour_gradient(low="#ADD8E6",high="blue",limits = c(0,8)) +
    #scale_colour_discrete(c("blue","red"))+
    #scale_colour_gradient(low="blue",high="red") + 
    #scale_shape_manual(values = c(16,10))+
    scale_size_area(name="genecounts")+theme_bw()+
    scale_size_continuous(range=c(1,8))+
    labs(y='',x='',title = i)+
    geom_vline(xintercept = 0.05,linetype =2,colour = 'black')+
    theme(axis.text=element_text(size=12),
          #axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title=element_text(size=18),
          panel.grid.major=element_line(colour=NA),
          plot.title = element_text(hjust = 0.5,size = 30),
          legend.text=element_text(size=12),legend.title = element_text(size=12))
  plotfile = paste('./enrichment/',i,'_res.pdf')
  ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 10)
}

# 6.heatmap ----------------------------------------------------------------------
#########heatmap_up and down
rm(list=ls())
library(pheatmap)
dif_data <- read.delim(file = './dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dim(dif_data)
dif_data<-dif_data[order(dif_data$pvalue),]
topgene<-c(dif_data$SYMBOL[dif_data$change=='UP'][1:20],dif_data$SYMBOL[dif_data$change=='DOWN'][1:20])
rownames(dif_data)<-dif_data$SYMBOL
gene_norm<-dif_data[dif_data$SYMBOL %in% topgene,1:9]
Case<-grep("case",colnames(gene_norm))
Control<-grep("control",colnames(gene_norm))
gene_norm<-gene_norm[,c(Case,Control)]
Z_score <- (gene_norm - apply(gene_norm, 1, mean)) / apply(gene_norm, 1, sd);

annotation_col <- data.frame(
  condition = c(rep("Case",length(Case)),rep("Control",length(Control)))
)
rownames(annotation_col) <- colnames(gene_norm);
ancols = c("#9BBB59","#7F7F7F")
names(ancols) = c("Case","Control")
ann_colors <- list(condition=ancols)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         clustering_distance_cols = "correlation",show_rownames=T,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=9,filename = "./result/gene_heatmap_v1.png")

########v2
rm(list=ls())
library(pheatmap)
dif_data <- read.delim(file = './dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data<-dif_data[order(dif_data$pvalue),]
rownames(dif_data)<-dif_data$SYMBOL
# gene_norm<-dif_data[dif_data$SYMBOL %in% topgene,1:7]
Case<-grep("case",colnames(dif_data))
Control<-grep("control",colnames(dif_data))
gene_norm<-dif_data[,c(Case,Control)]
Z_score <- (gene_norm - apply(gene_norm, 1, mean)) / apply(gene_norm, 1, sd);

annotation_col <- data.frame(
  condition = c(rep("Case",length(Case)),rep("Control",length(Control)))
)
rownames(annotation_col) <- colnames(gene_norm);
ancols = c("#9BBB59","#7F7F7F")
names(ancols) = c("Case","Control")
ann_colors <- list(condition=ancols)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         clustering_distance_cols = "euclidean",show_rownames=F,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=9,filename = "./result/gene_heatmap_v2.png")


# 7.差异p<0.01 --------------------------------------------------------------

rm(list = ls())
# t.test
edata <- read.delim("/home/devdata/jiarongf/nyy/tabfile/gene_matrix.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
pvalue<-c();tstat<-c();meansControl<-c();meansCase<-c();FC<-c();log2FC<-c()
Control_list<-c('control_060312_P1','control_192008_P2','control_192013_P1','control_192022_P1')
Case_list<-c('case_133076_P1','case_192005_P1','case_192010_P2','case_192020_P1')
edata<-edata+0.01

for (i in 1:nrow(edata)){
  #t.test
  Control_value<-edata[i,Control_list]
  Case_value<-edata[i,Case_list]
  result<-t.test(as.numeric(Case_value), as.numeric(Control_value), paired=FALSE);
  pvalue[i]<-result[["p.value"]]
  tstat[i]<-result[["statistic"]][["t"]]
  meansControl[i]<-mean(as.numeric(Control_value))
  meansCase[i]<-mean(as.numeric(Case_value))
  FC[i]<-mean(as.numeric(Case_value))/mean(as.numeric(Control_value))
  log2FC[i]<-log2(FC[i])
}
p_bf = p.adjust(pvalue, method = "bonferroni")
dif_data = data.frame(rownames(edata),pvalue,tstat,meansControl,meansCase,FC,log2FC,p_bf,stringsAsFactors=FALSE)

colnames(dif_data)[1] = 'SYMBOL'
rownames(dif_data) = dif_data$SYMBOL

#将表达矩阵同t.test结果合并
edata$SYMBOL<-rownames(edata)
dif_data_all = merge(edata,dif_data,by="SYMBOL",all =TRUE) ##all =TRUE
logFC_cutoff = log2(1)
dif_data_all$change =factor(ifelse(dif_data_all$pvalue < 0.01 & abs(dif_data_all$log2FC) >= logFC_cutoff,
                                   ifelse(dif_data_all$log2FC >= logFC_cutoff , 'UP', 'DOWN' ), 'NOT'),levels=c('UP', 'DOWN', 'NOT'))##总的差异的基因
write.table(diff_data_all, file = "./dif_data001.tab", quote = FALSE,sep="\t",row.names = FALSE) 



# 8.富集BPp<0.01 --------------------------------------------------------------------
rm(list=ls())
library( "clusterProfiler")
load(file = '/home/devdata/zhuzn/enrichment_annotation/crab_eating_enrichment.RData')
dif_data <- read.delim(file = './dif_data001.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
gene_up = dif_data[dif_data$change=='UP','SYMBOL']#12
gene_down = dif_data[dif_data$change=='DOWN','SYMBOL']#13
# gene_up = dif_data[dif_data$pvalue<0.01 & dif_data$FC>1, 'SYMBOL']#255
# gene_down = dif_data[dif_data$pvalue<0.01 & dif_data$FC<1, 'SYMBOL']#233
# dif_gene<-c(gene_up,gene_down)
###GOBP
term2gene <- gobp[,c('GO','SYMBOL')]
term2gene<-term2gene[!is.na(term2gene$SYMBOL),]
term2name <- gobp[,c('GO','Name')]

ego_BP <- enricher(
  gene = gene_up,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_UP_result<-as.data.frame(ego_BP@result)


ego_BP <- enricher(
  gene = gene_down,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
ego_BP_DOWN_result<-as.data.frame(ego_BP@result)



  library(ggplot2)
  res.filter1 <- ego_BP_UP_result[ego_BP_UP_result$pvalue<0.05,]
  res.filter1 <- res.filter1[order(res.filter1$pvalue,decreasing = F),]
  res.filter1 <- res.filter1[1:20,]
  res.filter1$Change<-'Up'
  
  res.filter2 <- ego_BP_DOWN_result[ego_BP_DOWN_result$pvalue<0.05,]
  res.filter2 <- res.filter2[order(res.filter2$pvalue,decreasing = F),]
  res.filter2 <- res.filter2[1:20,]
  res.filter2$Change<-'Down'
  
  dt <- rbind(res.filter1,res.filter2)
  dt$log.p <- -log10(dt$pvalue)
  dt$Change <- factor(dt$Change,levels = c('Up','Down'))
  dt$Count<-as.numeric(gsub('/[0-9]+','',dt$GeneRatio))
  
  p<- ggplot(dt,aes(log.p,y=reorder(Description,-pvalue)))+ #以pvalue和pathway为X轴和Y轴
    geom_point(aes(size=Count,color=Change))+ #点图大小和颜色数据
    labs(y='',x='-log10(Pvalue)')+
    scale_size_area(name="Gene count")+
    theme_bw()+
    geom_vline(xintercept = -log10(0.05),linetype =2,colour = 'black')+
    theme(axis.text=element_text(size=15,face = "bold", color="gray50"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=16),legend.title = element_text(size=15))
  
  ggsave("./result/GO_BP001.buble.png", plot=p, dpi = 600,width = 15, height = 10)
