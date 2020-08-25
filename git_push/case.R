setwd("/home/devdata/jiarongf/nyy/tabfile")
rm(list=ls())
# 1.数据合并-------------------------------------------------------------------------
rm(list=ls())
case_192010-P2_gene_abund <- read.table(file = "/home/devdata/jiarongf/nyy/tabfile/case_192010-P2_gene_abund.tab",fill = TRUE)#60676
# test <- NO_9_IPC_gene_abund[which(NO_9_IPC_gene_abund$Gene.ID == "ENSG00000228252"),]###重复行
# NO_9_IPC_gene_abund <- unique(NO_9_IPC_gene_abund)
control_060312-P1_gene_abund <- control_060312-P1_gene_abund[!duplicated(control_060312-P1_gene_abund$Gene.ID), ]#60675
gene_matrix <- NO_9_IPC_gene_abund[,c(1,2,9)]

res = c('NO_9_MSC','NO_10IPC','NO_10_MSC','NO_11IPC','NO_11_MSC')
for (i in res) {
  data <- read.delim(file = paste("/data/siwei/hisat2/stringtie/",i,"_gene_abund.tab",sep=''))
  data1 <- data[,c(1,9)]
  gene_matrix <- merge(x=gene_matrix,y=data1,by="Gene.ID",all =TRUE)
  gene_matrix <- gene_matrix
  # colnames(data1) <- paste(i,"_",colnames(data1),sep='')
  # data2 <- cbind(gene_matrix,data1)
  # gene_matrix <- data2
}
colnames(gene_matrix) <- c("Gene.ID","Gene.Name","NO_9_IPC",'NO_9_MSC','NO_10IPC','NO_10_MSC','NO_11IPC','NO_11_MSC')
rownames(gene_matrix) <- gene_matrix[,1]

gene_matrix1 <- gene_matrix[,-c(1,2)]
# gene_matrix1 <- gene_matrix1[rowSums(gene_matrix1) > 0, ]#29624
gene_matrix1 <- gene_matrix1[rowSums(gene_matrix1 > 0) != 0 ,]#29624

gene_matrix2 <- gene_matrix[rownames(gene_matrix1),]
gene_matrix_rmdup <- aggregate(gene_matrix2[,-c(1,2)],by = list(gene_matrix2$Gene.Name), FUN = mean)#29221
rownames(gene_matrix_rmdup) <- gene_matrix_rmdup[,1]
gene_matrix_rmdup <- gene_matrix_rmdup[,-1]
write.table(gene_matrix_rmdup, file = "./prepDE/gene_TPM_matrix.tab", quote = FALSE,sep="\t",row.names = TRUE)