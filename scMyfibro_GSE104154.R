
# GSE104154_d0_d21_sma_tm_Expr_raw from GSE104154
# mart_export from https://www.ensembl.org/biomart/martview/ed6c0e929a46a633a562d74c11afdc49

library(Seurat)
library(openxlsx)
GSE104154_SC<-read.csv('/Users/yihuawang/Downloads/GSE104154_d0_d21_sma_tm_Expr_raw.csv',sep = ',', row.names = 1) #reads in counts data from scRNA-sqe experiment
GSE104154_SC<-GSE104154_SC[,c(GSE104154_d0 $Barcode,GSE104154_d21$Barcode)]

biomart<-read.table('/Users/yihuawang/Downloads/mart_export (1).txt',sep = '\t',header = T)
GSE104154_SC$Gene.stable.ID<-rownames(GSE104154_SC)

GSE104154_SC<-merge(GSE104154_SC,biomart,by='Gene.stable.ID')
GSE104154_SC<-aggregate(GSE104154_SC[,2:5361],by=list(GSE104154_SC$Gene.name),FUN=mean)

colnames(GSE104154_SC)[c(1,5359,5360)]
rownames(GSE104154_SC)<-GSE104154_SC$Group.1
GSE104154_SC<-GSE104154_SC[-1,-1]


GSE104154_SC_col<-as.data.frame((colnames(GSE104154_SC)))

GSE104154_SC_col$order<-1:5360
colnames(GSE104154_SC_col)[1]<-'Barcode'
GSE104154_SC[1:10,1:10]
GSE104154_SC_col<-merge(GSE104154_SC_col,GSE104154_d0,by='Barcode',all = T)
GSE104154_SC_col<-merge(GSE104154_SC_col,GSE104154_d21,by='Barcode',all = T)
for (i in 1:5360){
  GSE104154_SC_col$Celltype[i]<- paste0(na.omit(GSE104154_SC_col$defined.x[i]),na.omit(GSE104154_SC_col$defined.y[i]))  
}



GSE104154_SC_col<-GSE104154_SC_col[order(GSE104154_SC_col$order),]
colnames(GSE104154_SC)<-GSE104154_SC_col$Celltype
GSE104154_SC[10000,1:5]
write.table(GSE104154_SC,'GSE104154_SC.txt',sep = '\t')
save(GSE104154_SC,file='scMyfiro_GSE104154.RData')

