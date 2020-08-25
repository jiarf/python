setwd('/home/devdata/jiarongf/ATAC-library-110225-AK8856')
################################part1#########################
rm(list = ls())
library(SnapATAC);

file.name = "Intestine_ATAC_ALLv2.snap"
x.sp = createSnap(
  file=file.name,
  sample="ATAC-test",
  num.cores=1
);

write.table(x.sp@barcode, file = "../barcode_list.txt", quote = FALSE,sep="\t",row.names = F)


str(x.sp@barcode)
showBinSizes(file.name);
# [1]5000 10000  50000 100000
x.sp = addBmatToSnap(x.sp, bin.size=100000, num.cores=1);
x.sp = makeBinary(x.sp, mat="bmat");

tmp<-x.sp@bmat
tmp<-tmp[rowSums(tmp)>0,]
barcode<-tmp@Dimnames[[1]]

x.sp = x.sp[which(x.sp@barcode %in% barcode), 1:dim(tmp)[2], mat = c("bmat", "pmat", "gmat"),
            drop = "missing"]


###########修改bed文件，gene id to gene name
#read bed table
library('clusterProfiler')
library('org.Mmu.eg.db')
data<-read.table("Macaca_mulatta.Mmul_10.100.bed",header=T,sep="\t",stringsAsFactors = F)#data.csv自定义包含输入ID的数据集
columns(org.Mmu.eg.db)
#  [1] "ACCNUM"       "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [10] "GO"           "GOALL"        "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PMID"         "REFSEQ"       "SYMBOL"       "UNIPROT" 
#class(data)

#data_ENSEMBL <- data[,4]
#write.table(data_ENSEMBL, file = "./data_ENSEMBL.txt", quote = FALSE,sep="\t",row.names = F)#DAVID中进行了转化
#data_list<- read.table("data_ENSEMBL.txt",header=F,sep="\t",stringsAsFactors = F)

##function bitr to make the id to name,but no fully
gene_name <- bitr(data[,4], fromType = "ENSEMBL", toType="SYMBOL",OrgDb = org.Mmu.eg.db)#Warning message:In bitr(data[, 4], fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mmu.eg.db) :57.5% of input gene IDs are fail to map...



#load(file = '/home/devdata/zhuzn/enrichment_annotation/crab_eating_enrichment.RData')

##merge name and ensemble
colnames(data)<-c("CHR","START","END","ENSEMBL")
gene_data_bed<-merge(data,gene_name,by="ENSEMBL")
gene_bed <- gene_data_bed[,-1]
# bed <- as.data.frame(gene_data_bed$SYMBOL)   
# bed$start <- gene_data_bed$start
# bed$end <- gene_data_bed$end

write.table(gene_bed, file = "./Macaca_mulatta.Mmul_10.100_genename.bed", quote = FALSE,sep="\t",row.names = F)



#使用diffusion maps的方法进行数据降维
#tmp <- x.sp@bmat
#tmp<-tmp[rowSums(tmp)>0,]

#x.sp@bmat <- tmp

x.sp = runDimReduct(
  obj=x.sp,
  pc.num = 50,
  input.mat="bmat"
);

#6. Determine significant components
#接下来，我们基于数据降维的结果确定用于下游分析的维数。我们绘制不同维数之间的配对散点图，并选择散点图开始看起来零散的维数。在下面的示例中，我们选择了前20个维度
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
);
# Step 7. Graph-based clustering ------------------------------------------
#选择好有效的降维维度后，我们基于它们来构造一个K近邻(KNN)的聚类图(K =15)。在该图中，每个点代表一个细胞，并根据欧氏距离确定每个细胞的k近邻个点
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
);

x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
);

x.sp@metaData$cluster = x.sp@cluster;
# Step 8. Visualization ---------------------------------------------------
#SnapATAC可以使用tSNE(FI-tsne)或UMAP的方法对降维聚类后的结果进行可视化展示。在此例中，我们计算t-SNE embedding，使用tSNE方法进行可视化展示。同时，我们还将测序深度或其他偏差投射到t-SNE embedding上
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="Rtsne",
  seed.use=10
);

par(mfrow = c(2, 2));
plotViz(
  obj=x.sp,
  method="tsne", 
  main="10X Brain Cluster",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);

plotViz(
  obj=x.sp,
  method="tsne", 
  main="10X Brain TN",
  point.color=x.sp@metaData$TN, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);
plotViz(
  obj=x.sp,
  method="tsne", 
  main="10X Brain um",
  point.color=x.sp@metaData$UM, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);
plotViz(
  obj=x.sp,
  method="tsne", 
  main="10X Brain pp",
  point.color=x.sp@metaData$PP, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);




# Step 9. Gene based annotation -------------------------------------------

#为了帮助注释聚类后分群的细胞簇，SnapATAC接下来将创建cell-by-gene计数矩阵，
#并可视化marker基因的富集情况，根据marker基因的表达情况进行cluster的注释

genes = read.table("Macaca_mulatta.Mmul_10.bed",header = F);
genes.gr = GRanges(genes[,1], 
                     IRanges(genes[,2], genes[,3]), name=genes[,4]
);
marker.genes = c(
  "CCNF", "TEDC2", "NAGA",
  "PTPRC", "EMC7", "KATNBL1", 
  "ADSS", "MS4A15", "MS4A10"
);

genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];

# re-add the cell-by-bin matrix to the snap object;
x.sp = addBmatToSnap(x.sp);
x.sp = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.gr,
  do.par=TRUE,
  num.cores=10
);



# smooth the cell-by-gene matrix
x.sp = runMagic(
  obj=x.sp,
  input.mat="gmat",
  step.size=3
);

par(mfrow = c(3, 3));
for(i in 1:9){
  plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@gmat[, marker.genes[i]],
    method="tsne", 
    main=marker.genes[i],
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0, 1)
  )};

####################cell by gene
tmp<-x.sp@gmat
tmp<-tmp[rowSums(tmp)>0,]
tmp@Dimnames[[1]] <- barcode

x.sp = x.sp[which(x.sp@barcode %in% barcode), 1:dim(tmp)[2], mat = c("bmat", "pmat", "gmat"),
            drop = "missing"]
tmp <- as.data.frame(tmp)
tmp_t <- as.data.frame(t(tmp))

tmp_t<-tmp_t[rowSums(tmp_t)>0,]###38 obs 185 variables
tmp_t$sum <- rowSums(tmp_t)
#write.table(tmp, file = "./x.sp_gmat.txt", quote = FALSE,sep="\t",row.names = F) 
tmp_t_order <- tmp_t[order(tmp_t$sum,decreasing = T),]
# tmp_t_order$sum
# [1] 30.7570380  6.2040239  3.3310750  2.5315693  1.8676676  1.8136436  1.8123016  1.8087290  1.8086460  1.8085857
# [11]  1.7641256  1.5665697  1.5325065  1.5313160  1.4985924  1.4370979  1.4028365  1.3692089  1.3378944  1.2747906
# [21]  1.2745018  1.2075544  1.1352784  1.1323443  1.1319765  0.6362299  0.6036635  0.6035875  0.6028820  0.6001295
# [31]  0.5676392  0.5656764  0.5654663  0.5324709  0.5324709  0.5324709  0.5324709  0.5324709
top10 <- rownames(tmp_t_order[1:10,])
top10
#[1] "TRPC3"    "FOLR2"    "NALCN"    "ADD3"     "PRKCB"    "WDHD1"    "ANGPTL5"  "ZNF222"   "SLC25A51" "FSTL5" 
marker_gene <- top10


# Step 10. Heretical clustering -------------------------------------------

#接下来，我们将属于同一cluster的细胞汇集到一起，用以创建每个cluster的聚合信号（aggregate signal）。
# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})

# cluster using 1-cor as distance  
# data <- 1 - cor(t(do.call(rbind, ensemble.ls)))#matrix
# data
# #1        2        3        4  5
# #1  0       NA       NA       NA NA
# #2 NA 0.000000 1.491652 1.178451 NA
# #3 NA 1.491652 0.000000 1.550227 NA
# #4 NA 1.178451 1.550227 0.000000 NA
# #5 NA       NA       NA       NA  0
# ##检测缺失值
# library(VIM)
# aggr(data,prop = T,numbers=T)
# ##处理缺失值，
# data_cut_NA <- data[-1,]
# data_cut_NA <- data_cut_NA[-4,]
# data_cut_NA <- data_cut_NA[,-1]
# data_cut_NA <- data_cut_NA[,-4]
# data_cut_NA <- as.dist(data_cut_NA)
hc_2 = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))) ,method="ward.D2");
# hc_2 = hclust(data_cut_NA ,method="ward.D2");
plotViz(
  obj=x.sp,
  method="tsne", 
  main="10X Brain Cluster",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);
plot(hc_2, hang=-1, xlab="");
#Step 11. Identify peaks -------------------------------------------------
  #接下来，我们将每个cluster群的细胞信息聚合起来，创建一个用于peak calling和可视化的集成track。
  #在该步骤中，将生成一个.narrowPeak的文件，其中包含识别出的所有peak的信息，和一个.bedGraph文件，可以用于可视化展示。
  #为了获得最robust的结果，我们不建议对细胞数目小于100的cluster执行此步骤。
  # 查看snaptools的安装路径
  system("which snaptools")
#/home/devdata/Software/soft/anaconda3/bin/snaptools
# 查看macs2的安装路径
system("which macs2")
#/home/devdata/Software/soft/anaconda3/bin/macs2

# 调用macs2进行peak callling
runMACS(
  obj=x.sp[which(x.sp@cluster==1),],  
  output.prefix="Intestine_ATAC_ALLv2.1",
  path.to.snaptools="/home/devdata/Software/soft/anaconda3/bin/snaptools",
  path.to.macs="/home/devdata/Software/soft/anaconda3/bin/macs2",
  gsize="hs",
  buffer.size=500,
  num.cores=1,
  macs.options="--nomodel --shift 75 --ext 150 --qval 1e-2 -B --SPMR --call-summits",
  tmp.folder=tempdir(),
  keep.minimal = TRUE
);
#Build Peak Model... 
#INFO  @ Tue, 04 Aug 2020 16:53:35: #2 Skipped... 
#但是对于双端测序而言，它本身测的就是文库的两端，因此建立模型没有必要，偏倚也没有必要
#NAME_peaks.narrowPeak NAME_peaks.broadPeak 类似。后面4列表示为， integer score for display， fold-change，-log10pvalue，-log10qvalue，relative summit position to peak start。内容和NAME_peaks.xls基本一致，适合用于导入R进行分析。


# call peaks for all cluster with more than 100 cells
clusters.sel = names(table(x.sp@cluster));
class(clusters.sel)
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("Intestine_ATAC_ALLv2.1.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/home/devdata/Software/soft/anaconda3/bin/snaptools",
    path.to.macs="/home/devdata/Software/soft/anaconda3/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=1,
    macs.options="--nomodel --shift 37 --ext 73 --qval 5e-2 -B --SPMR",
    tmp.folder=devdata()
  );
}, mc.cores=5);

# assuming all .narrowPeak files in the current folder are generated from the clusters
peaks.names = system("ls | grep narrowPeak", intern=TRUE);
class(peaks.names)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})



