##### DropViz data set http://dropviz.org/ #####
### orignal study: Saunders et al. bioRxiv. 2018
### aggregated gene expression profiles for 565 cell population was downloaded from http://dropviz.org/
### Note that we did not use individual cell level expression value
### Two levels of cell type annotation was used
### level 1: "class" column in the annotation file
### level 2: "subcluster" clumn in the annotation file
### Note that we kept cell types from different brain regions separate
### e.g. Neurons from Cerebellum and Frontal Cortex were kept separate as CB.neurons and FC.Neurons
### For level 2, we used the column "full_name" but also concatenated "tissue" and "subcluster" to make 565 unique labels
### We also created per brain region cell type expression data set for both level 1 and level 2
### In total (1 (all together) + 9 (brain regions)) x 2 (levels) = 20 outputs

library(data.table)
library(Hmisc)

##### metadata #####
samples <- fread(input="gzip -cd DropViz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv.gz", data.table=F, header=T)
samples$class <- capitalize(tolower(samples$class))
samples$subcluster <- sub("-", "_", samples$subcluster)
samples$full_name <- gsub(" ", "_", gsub(",", "", gsub("\\.", "_", samples$full_name)))

##### expression #####
umi <- fread(input="gzip -cd DropViz/metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv.gz", data.table=F)
genes <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- genes
### convert to cpm since the downloadable umi is the sum of umi counts in each subcluster
lib_size <- apply(umi, 2, sum)/1000000
umi <- t(apply(umi, 1, function(x){x/lib_size}))

### map genes to hs ENSG ID
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(umi), stringsAsFactors = F)
genes$mm.ensg <- ncbi_mouse$ensg[match(genes$symbol, ncbi_mouse$Symbol)]
ncbi_mouse$Synonyms <- paste0("|", ncbi_mouse$Synonyms, "|")
genes$mm.ensg[is.na(genes$mm.ensg)] <- sapply(genes$symbol[is.na(genes$mm.ensg)], function(x){
  n <- grep(paste0("|",x,"|"), ncbi_mouse$Synonyms)
  if(length(n)==1){ncbi_mouse$ensg[n]}
  else{NA}
})
genes$hs.ensg <- mm2hs$hs.ensg[match(genes$mm.ensg, mm2hs$mm.ensg)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
umi <- umi[genes$symbol[!is.na(genes$hs.ensg)],]

### log transform
log_umi <- log2(umi+1)

###### level 1 in all regions #####
cell_label <- unique(paste(gsub(" ","_",samples$tissue), gsub(" ", "_",samples$class), sep="."))
ts_cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
colnames(ts_cell_log_umi) <- cell_label
ts <- unique(samples$tissue)
for(t in ts){
  ct <- unique(samples$class[samples$tissue==t])
  for(c in ct){
    n <- which(samples$tissue==t & samples$class==c)
    if(length(n)>1){
      ts_cell_log_umi[,paste(gsub(" ","_",t), gsub(" ", "_",c), sep=".")] <- apply(log_umi[,samples$tissue_subcluster[n]], 1, mean)
    }else{
      ts_cell_log_umi[,paste(gsub(" ","_",t), gsub(" ", "_",c), sep=".")] <- log_umi[,samples$tissue_subcluster[n]]
    }
  }
}
ts_cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], ts_cell_log_umi)
ts_cell_log_umi$Average <- apply(ts_cell_log_umi[,2:ncol(ts_cell_log_umi)], 1, mean)
write.table(ts_cell_log_umi, "processed_data/DropViz_all_level1.txt", quote=F, row.names=F, sep="\t")

##### level 2 in all regions #####
ts_cell_log_umi <- data.frame(log_umi)
colnames(ts_cell_log_umi) <- paste(samples$tissue, samples$full_name, samples$subcluster, sep=".")
ts_cell_log_umi <- cbind(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], ts_cell_log_umi)
ts_cell_log_umi$Average <- apply(ts_cell_log_umi[,2:ncol(ts_cell_log_umi)], 1, mean)
write.table(ts_cell_log_umi, "processed_data/DropViz_all_level2.txt", quote=F, row.names=F, sep="\t")

##### level 1 per brain region #####
for(t in ts){
  tmp_umi <- log_umi[, samples$tissue_subcluster[samples$tissue==t]]
  ct <- unique(samples$class[samples$tissue==t])
  cell_log_umi <- matrix(nrow=nrow(tmp_umi), ncol=length(ct))
  colnames(cell_log_umi) <- ct
  for(c in ct){
    n <- which(samples$tissue==t & samples$class==c)
    if(length(n)>1){
      cell_log_umi[,c] <- apply(tmp_umi[,samples$tissue_subcluster[n]], 1, mean)
    }else{
      cell_log_umi[,c] <- tmp_umi[,samples$tissue_subcluster[n]]
    }
  }
  cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(tmp_umi), genes$symbol)], cell_log_umi)
  cell_log_umi$Average <- apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
  write.table(cell_log_umi, paste0("processed_data/DropViz_", gsub(" ", "_", t), "_level1.txt"), quote=F, row.names=F, sep="\t")
}

##### level 2 per brain region #####
for(t in ts){
  cell_log_umi <- data.frame(log_umi[,samples$tissue_subcluster[samples$tissue==t]])
  ct <- unique(paste(samples$full_name[samples$tissue==t], samples$subcluster[samples$tissue==t], sep="."))
  colnames(cell_log_umi) <- ct
  cell_log_umi <- cbind(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
  cell_log_umi$Average <- apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
  write.table(cell_log_umi, paste0("processed_data/DropViz_", gsub(" ", "_", t), "_level2.txt"), quote=F, row.names=F, sep="\t")
}
