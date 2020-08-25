setwd('/home/devdata/zr/proj/Transcriptome/nyy')

# 1.数据合并-------------------------------------------------------------------------
rm(list=ls())
P4_case_gene_abund <- read.delim("/home/devdata/zr/proj/Transcriptome/nyy/1.hisat2/tabfile/P4_case_gene_abund.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)#60676
###colnames(P4_case_gene_abund)
gene_matrix <- P4_case_gene_abund[,c(1,2,9)]
colnames(gene_matrix) <- c('Gene.ID','Gene.Name','P4_case')

res = c('P4_control','P7_case','P7_control','P10_case','P10_control')
for (i in res) {
  data <- read.delim(file = paste("/home/devdata/zr/proj/Transcriptome/nyy/1.hisat2/tabfile/",i,"_gene_abund.tab",sep=''))
  data <- data[,c(1,9)]
  colnames(data) <- c('Gene.ID',i)
  gene_matrix <- merge(x=gene_matrix,y=data,by="Gene.ID",all =TRUE) ##all =TRUE
}

gene_matrix<-gene_matrix[gene_matrix$Gene.Name !='-',]
gene_matrix<-gene_matrix[gene_matrix$Gene.Name !='.',]
exp_data <- aggregate(gene_matrix[,c(-1,-2)],by = list(gene_matrix$Gene.Name), FUN = mean)
rownames(exp_data)<-exp_data[,1]
exp_data<-exp_data[,-1]
exp_data<-exp_data[rowSums(exp_data)>0,]
write.table(exp_data, file = "./data/exp_data.tab", quote = FALSE,sep="\t",row.names = TRUE)

# 2.PCA -------------------------------------------------------------------
###2.1 PCA.R
rm(list=ls())
library(ropls)
source("PCA.R")
exp_data <- read.delim("./data/exp_data.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
sampledata<-read.delim("./data/metadata.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)
colnames(sampledata) <- c("SampleID","Group")
color=list(case = "#fc8d59",control = "#99d594")
dir.create("./result/")
PCA_analysis(exp_data,pData=sampledata,color=color, "casevscontrol", group_name='Group',filepath="./result/")

###2.2 PCAtools
rm(list=ls())
library(PCAtools)
exp_data <- read.delim("./data/exp_data.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
sampledata<-read.delim("./data/metadata.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE)
rownames(sampledata) = sampledata[,1]
p <- pca(exp_data, metadata = sampledata)
pdf('./result/PCA_PCAtool.pdf', width = 10, height = 8)
biplot(p,colby = 'Group',title = 'PCA',pointSize = 8,labSize = 6.0,labhjust = 1.0,labvjust = 1.5,legendPosition = 'right')
dev.off()


# 3. DEGs -----------------------------------------------------------------
# t.test
edata <- read.delim("./data/exp_data.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
pvalue<-c();tstat<-c();meansControl<-c();meansCase<-c();FC<-c();log2FC<-c()
Control_list<-c('P10_control','P4_control','P7_control')
Case_list<-c('P10_case','P4_case','P7_case')
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
write.table(dif_data_all, file = "./data/dif_data.tab", quote = FALSE,sep="\t",row.names = FALSE)

# 4差异表达分析-火山图 ----------------------------------------------------
library(ggplot2)
library(ggrepel)
rm(list=ls())
dif_data <- read.delim(file = './data/dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data<-dif_data[order(dif_data$pvalue),]
topgene<-c(dif_data$SYMBOL[dif_data$change=='UP'][1:20],dif_data$SYMBOL[dif_data$change=='DOWN'][1:20])
sig<-rep('',nrow(dif_data))
for (i in 1:nrow(dif_data)){
  if (dif_data$SYMBOL[i] %in% topgene){
    sig[i]<-'most_dif'
  } else if (dif_data$pvalue[i]<0.01){
    sig[i]<-'True'
  } else{
    sig[i]<-'False'
  }
}

dif_data$significant <- as.factor(sig)
dif_data$change <- factor(dif_data$change,levels=c('UP', 'DOWN', 'NOT'))
logFC_cutoff = log2(1)
dif_data$change1 =factor(ifelse(dif_data$pvalue < 0.01 & abs(dif_data$log2FC) >= logFC_cutoff,
                               ifelse(dif_data$log2FC >= logFC_cutoff , 'UP', 'DOWN' ), 'NOT'),levels=c('UP', 'DOWN', 'NOT'))
this_tile <- paste0(nrow(dif_data[dif_data$change =='UP', ] ),' Up-regulated genes', "; ",
                    nrow(dif_data[dif_data$change =='DOWN', ]), ' Down-regulated genes')

p=ggplot(data=dif_data, aes(x=log2FC, y =-log10(pvalue),color=change))+
   geom_point(alpha=0.8, size=1.2)+
  scale_x_continuous(limits = c(-9,9),breaks = seq(-9,9,1))+
  scale_color_manual(values =c('red',"blue","#595959"))+
  labs(title="", x="log2 (fold change)",y="-log10 (p-value)")+
  ggtitle( this_tile )+
  theme( plot.title = element_text(size =15, hjust = 0.5)) +
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  theme_set( theme_set(theme_bw( base_size = 12)))+
  geom_text_repel(data=subset(dif_data, significant =='most_dif'), aes(label=SYMBOL),col="black",size=2)
plotfile=paste('./result/','dif_data','_Volcano.png',sep='')
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)
dif_data_most<-dif_data[dif_data$significant=='most_dif',]
write.table(dif_data_most, file = "./data/dif_data_mosttop20.tab", quote = FALSE,sep="\t",row.names = FALSE)


# 4 .heatmap --------------------------------------------------------------

#########heatmap_v2
rm(list=ls())
library(pheatmap)
dif_data <- read.delim(file = './data/dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dim(dif_data)
dif_data<-dif_data[order(dif_data$pvalue),]
topgene<-c(dif_data$SYMBOL[dif_data$change=='UP'][1:20],dif_data$SYMBOL[dif_data$change=='DOWN'][1:20])
rownames(dif_data)<-dif_data$SYMBOL
gene_norm<-dif_data[dif_data$SYMBOL %in% topgene,1:7]
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
######v2
rm(list=ls())
library(pheatmap)
dif_data <- read.delim(file = './data/dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
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

# 5.kegg_GO ---------------------------------------------------------------

setwd('/home/devdata/zr/proj/Transcriptome/nyy')
# 5.kegg_GO ---------------------------------------------------------------
rm(list=ls())
library( "clusterProfiler")
load(file = '/home/devdata/zhuzn/enrichment_annotation/crab_eating_enrichment.RData')
dif_data <- read.delim(file = './data/dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
gene_up = dif_data[dif_data$change=='UP','SYMBOL']#1187
gene_down = dif_data[dif_data$change=='DOWN','SYMBOL']#882
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
ego_BP_up_result<-as.data.frame(ego_BP@result)

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
ego_BP_down_result<-as.data.frame(ego_BP@result)

write.table(ego_BP_down_result, file = "./result/ego_BP_down_result.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_up_result, file = "./result/ego_BP_up_result.tab", quote = FALSE,sep="\t",row.names = FALSE)


####plot GO
library(ggplot2)
res.filter1 <- ego_BP_up_result[ego_BP_up_result$pvalue<0.05,]
res.filter1 <- res.filter1[order(res.filter1$pvalue,decreasing = F),]
res.filter1 <- res.filter1[1:20,]
res.filter1$Change<-'Up'

res.filter2 <- ego_BP_down_result[ego_BP_down_result$pvalue<0.05,]
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

ggsave("./result/GO_BP.buble.png", plot=p, dpi = 600,width = 15, height = 10)
ggsave("./result/GO_BP.buble.pdf", plot=p, dpi = 600,width = 15, height = 10)


# ###kegg -----------------------------------------------------------------

term2gene <- kegg[,c(2,6)]
term2gene<-term2gene[term2gene$SYMBOL !='',]
term2name <- kegg[,c(2,3)]

kegg_en <- enricher(
  gene = gene_up,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2gene
)
kegg_result<-as.data.frame(kegg_en@result)


############## p4_vs_p7_vs_p10 ################
# 1. DEGs -----------------------------------------------------------------
# t.test
rm(list=ls())
edata <- read.delim("./data/exp_data.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
p4vsp7<-c();p7vsp10<-c();p4vsp10<-c();log2p4vsp7<-c();log2p7vsp10<-c();log2p4vsp10<-c();control_not<-c();P4vsP7vsP10<-c();
P4_list <- 'P4_case'
P7_list <- 'P7_case'
P10_list <- 'P10_case'
contorl_log2FC <- 0
n_log2FC <- 0
edata<-edata+0.01

for (i in 1:nrow(edata)){
  P4 <- edata[i,P4_list]
  P7 <- edata[i,P7_list]
  P10 <- edata[i,P10_list]
  P4_control <- edata[i,"P4_control"]
  P7_control <- edata[i,"P7_control"]
  P10_control <- edata[i,"P10_control"]
  p4vsp7[i] <- as.numeric(P4)/as.numeric(P7)
  p7vsp10[i] <- as.numeric(P7)/as.numeric(P10)
  p4vsp10[i] <- as.numeric(P4)/as.numeric(P10)
  log2p4vsp7[i] <- log2(p4vsp7[i])
  log2p7vsp10[i] <- log2(p7vsp10[i])
  log2p4vsp10[i] <- log2(p4vsp10[i])
  control_not[i] = ifelse((abs(P4_control/P7_control) <= contorl_log2FC) & (abs(P7_control/P10_control) <= contorl_log2FC) & (abs(P4_control/P10_control) <= contorl_log2FC), 'NOT', 'DIFF' )
  #if(control_not[i] == "NOT"){
    if(log2p4vsp7[i] > n_log2FC & log2p7vsp10[i] > n_log2FC){
    #if(log2p4vsp7[i] > n_log2FC & log2p4vsp10[i] > n_log2FC){
    P4vsP7vsP10[i] = "UP_UP"} else if(log2p4vsp7[i] > n_log2FC & abs(log2p7vsp10[i]) <= n_log2FC) {
      P4vsP7vsP10[i] = "UP_NOT"} else if(abs(log2p4vsp7[i]) <= n_log2FC & log2p7vsp10[i] >n_log2FC){
        P4vsP7vsP10[i] = "NOT_UP"} else if(log2p4vsp7[i] > n_log2FC & log2p7vsp10[i] < -(n_log2FC)){
          #P4vsP7vsP10[i] = "UP_DOWN"} else if(log2p4vsp7[i] < -(n_log2FC) & log2p7vsp10[i] < -(n_log2FC)){
          P4vsP7vsP10[i] = "UP_DOWN"} else if(log2p4vsp7[i] < -(n_log2FC) & log2p4vsp10[i] < -(n_log2FC)){
            P4vsP7vsP10[i] = "DOWN_DOWN"} else if(log2p4vsp7[i] < -(n_log2FC) & abs(log2p7vsp10[i]) <= n_log2FC){
              P4vsP7vsP10[i] = "DOWN_NOT"} else if(abs(log2p4vsp7[i]) <= n_log2FC & log2p7vsp10[i] < -(n_log2FC)){
                P4vsP7vsP10[i] = "NOT_DOWN"} else if(log2p4vsp7[i] < -(n_log2FC) & log2p7vsp10[i] > n_log2FC){
                  P4vsP7vsP10[i] = "DOWN_UP"} else {P4vsP7vsP10[i] = "NA"}
  #}
  #else{P4vsP7vsP10[i] = "NA"}
}
#p_bf = p.adjust(pvalue, method = "bonferroni")
dif_data = data.frame(rownames(edata),p4vsp7,p7vsp10,p4vsp10,log2p4vsp7,log2p7vsp10,log2p4vsp10,control_not,P4vsP7vsP10,stringsAsFactors=FALSE)

colnames(dif_data)[1] = 'SYMBOL'
rownames(dif_data) = dif_data$SYMBOL

write.table(dif_data,file = "result/all_UP_DOWN.xls",sep = "\t",quote = F,row.names = F)
table(dif_data$control_not == "NOT")#1504
table(dif_data$P4vsP7vsP10)

################ heatmap #################
library("pheatmap")
library("ggplot2")
library("dplyr")
rm(list=ls())
edata <- read.delim("./data/exp_data.tab",sep='\t',header=TRUE,stringsAsFactors=FALSE)
all_UP_DOWN <- read.delim("/home/devdata/zr/proj/Transcriptome/nyy/result/all_UP_DOWN.xls")
dif_data <- read.delim(file = './data/dif_data.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
df <- all_UP_DOWN[which(all_UP_DOWN$P4vsP7vsP10 == "DOWN_DOWN"),]
genelist1 <- df$SYMBOL
genelist2 <- dif_data$SYMBOL[which(dif_data$change == "UP" | dif_data$change == "DOWN")]
genelist <- intersect(genelist1,genelist2)
table(dif_data$change)
# DOWN   NOT    UP 
# 882 13660  1187 
table(all_UP_DOWN$P4vsP7vsP10)
# DOWN_DOWN  DOWN_NOT   DOWN_UP        NA  NOT_DOWN    NOT_UP   UP_DOWN    UP_NOT     UP_UP 
# 3        84       170     15146        85        88        49        94        10
edata <- edata[genelist,]
zsorse <- (edata - apply(edata, 1, mean)) / apply(edata, 1, sd)
#zsorse <- scale(log10(edata+0.00000001))
#zsorse[zsorse >1] =2
#zsorse[zsorse < -1] = -2
zsorse <- zsorse[,c(2,4,6,1,3,5)]
edata['HLA-F',]
zsorse['HLA-F',]
p <- pheatmap(zsorse,border_color="#C5C5C5",color = colorRampPalette(c("#0F4FA8", "white", "#eb4d4b"))(50),
         cluster_cols = FALSE,
         clustering_distance_cols = "euclidean",show_rownames=T,
         #annotation_col=annotation_col,
         annotation_colors=ann_colors,
         #filename = "./result/gene_heatmap_v2.png",
         fontsize=9
         )
ggsave("case_up_down.heatmap.pdf",plot = p,width = 8,height = 12)
#########################################
library( "clusterProfiler")
load(file = '/home/devdata/zhuzn/enrichment_annotation/crab_eating_enrichment.RData')
#query_data = dif_data$SYMBOL[dif_data$P4vsP7vsP10 != "NA" & dif_data$P4vsP7vsP10 != "UP_DOWN" & dif_data$P4vsP7vsP10 != "DOWN_UP"]

term2gene <- gobp[,c('GO','SYMBOL')]
term2gene<-term2gene[!is.na(term2gene$SYMBOL),]
term2name <- gobp[,c('GO','Name')]

ego_BP <- enricher(
  gene = query_data,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
query_BP_result<-as.data.frame(ego_BP@result)
write.table(query_BP_result, file = "./result/query_BP_result.tab", quote = FALSE,sep="\t",row.names = FALSE)

res.filter1 <- query_BP_result[query_BP_result$pvalue<0.05,]
res.filter1 <- res.filter1[order(res.filter1$pvalue,decreasing = F),]
res.filter1 <- res.filter1[1:20,]

dt <-res.filter1
dt$log.p <- -log10(dt$pvalue)
dt$Count<-as.numeric(gsub('/[0-9]+','',dt$GeneRatio))


p<- ggplot(dt,aes(log.p,y=reorder(Description,-pvalue)))+ #以pvalue和pathway为X轴和Y轴
  geom_point(aes(size=Count,color=Count))+ #点图大小
  labs(y='',x='-log10(Pvalue)')+
  scale_size_area(name="Gene count")+
  theme_bw()+
  geom_vline(xintercept = -log10(0.05),linetype =2,colour = 'black')+
  #scale_color_brewer(palette = 1)+
  theme(axis.text=element_text(size=15,face = "bold", color="gray50"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),legend.title = element_text(size=15))

ggsave("./result/query_GO_BP.buble.png", plot=p, dpi = 600,width = 15, height = 10)
ggsave("./result/query_GO_BP.buble.pdf", plot=p, dpi = 600,width = 15, height = 10)