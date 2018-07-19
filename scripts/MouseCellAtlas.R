##### Mouse Cell Atlas data set http://bis.zju.edu.cn/MCA/ #####
### original study: Ham et al. Cell. 2018 (PMID: 29474909)
### Both expression and cell labels were downloaded from from the Mouse Cell Atlas website at http://bis.zju.edu.cn/MCA/contact.html
### Used batch removed data sets used for Figure 2 in the original study
### Data is also available at GSE108097
### We created multiple data sets from MCA
### 1) all together (731 unique cell types)
### 2) only adult mouse samples (437 cell types from 18 tissue)
### 3) only embryo samples (including fetal tissues, 137 cell types)
### 4) only neonatal samples (108 cell types)
### 5) per tissue per sample type (adult, emmbryo, neonatal and cell line, 33 combinatioin in total)
### In total 1 (all) + 1 (adult) + 1 (embryo) + 1 (neonatal) + 33 (per tissue per type) = 37 data sets

library(data.table)

##### samples and expression #####
samples <- fread("MouseCellAtlas/cell-ann.csv", data.table=F)
umi <- fread(input="gzip -cd MouseCellAtlas/Figure2-batch-removed.txt.gz", data.table=F, header=T, fill=T)
cols <- colnames(umi)[1:(ncol(umi)-1)]
colnames(umi)[1] <- "Gene"
colnames(umi)[2:ncol(umi)] <- cols
g <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- g
samples <- samples[samples$Cell.name %in% colnames(umi),]
umi <- umi[,samples$Cell.name]

##### tissue cell type sort #####
samples$Tissue <- gsub("\\.", "_", gsub("-", "_", samples$Tissue))
samples$Tissue <- sub("MammaryGland", "Mammary_Gland", samples$Tissue)
samples$ts <- samples$Tissue
samples$ts[grepl("^Bone_Marrow", samples$Tissue)] <- "Bone_Marrow"
samples$ts[grepl("^Mammary_Gland", samples$Tissue)] <- "Mammary_Gland"
samples$cell_type <- sub("\\(.+\\)", "", samples$Cell.Anno)
samples$type <- samples$`Mouse-Sex-Age`
samples$type[grepl("week", samples$type)] <- "Adult"
samples$type[grepl("pups", samples$type)] <- "Neonatal"
samples$type[grepl("embryo", samples$type)] <- "Embryo"
samples$type[grepl("cell line", samples$type)] <- "Cell_line"

##### Map genes to hs ENSG ID #####
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
notes <- c(notes, paste0("Number of genes with unique hs ENSG: ", nrow(umi)))

##### process #####
log_umi <- log2(umi+1)

### all tissues
tmp <- unique(samples[c("type", "ts", "cell_type")])
tmp$type <- factor(tmp$type, levels=c("Adult", "Embryo", "Neonatal", "Cell_line"))
tmp <- tmp[with(tmp, order(type, ts, cell_type)),]
cell_label <- unique(paste(tmp$ts, gsub(" ","_",tmp$cell_type), sep=":"))
ts_cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
colnames(ts_cell_log_umi) <- cell_label
ts <- unique(samples$ts)
for(t in ts){
  ct <- unique(samples$cell_type[samples$ts==t])
  for(c in ct){
    if(length(which(samples$ts==t & samples$cell_type==c))==1){
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]]
    }else{
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = apply(log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]], 1, mean)
    }
  }
}
ts_cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], ts_cell_log_umi)
ts_cell_log_umi$Average <- apply(ts_cell_log_umi[,2:ncol(ts_cell_log_umi)], 1, mean)

write.table(ts_cell_log_umi, "processed_data/MouseCellAtlas_all.txt", quote=F, row.names=F, sep="\t")

### only adult samples
cell_label <- sort(unique(paste(samples$ts[samples$type=="Adult"], gsub(" ","_",samples$cell_type[samples$type=="Adult"]), sep=":")))
ts_cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
colnames(ts_cell_log_umi) <- cell_label
ts <- unique(samples$ts[samples$type=="Adult"])
for(t in ts){
  ct <- unique(samples$cell_type[samples$ts==t])
  for(c in ct){
    if(length(which(samples$ts==t & samples$cell_type==c))==1){
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]]
    }else{
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = apply(log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]], 1, mean)
    }
  }
}
ts_cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], ts_cell_log_umi)
ts_cell_log_umi$Average <- apply(ts_cell_log_umi[,2:ncol(ts_cell_log_umi)], 1, mean)

write.table(ts_cell_log_umi, "processed_data/MouseCellAtlas_Adult_all.txt", quote=F, row.names=F, sep="\t")

### only embryo samples
cell_label <- sort(unique(paste(samples$ts[samples$type=="Embryo"], gsub(" ","_",samples$cell_type[samples$type=="Embryo"]), sep=":")))
ts_cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
colnames(ts_cell_log_umi) <- cell_label
ts <- unique(samples$ts[samples$type=="Embryo"])
for(t in ts){
  ct <- unique(samples$cell_type[samples$ts==t])
  for(c in ct){
    if(length(which(samples$ts==t & samples$cell_type==c))==1){
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]]
    }else{
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = apply(log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]], 1, mean)
    }
  }
}
ts_cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], ts_cell_log_umi)
ts_cell_log_umi$Average <- apply(ts_cell_log_umi[,2:ncol(ts_cell_log_umi)], 1, mean)

write.table(ts_cell_log_umi, "processed_data/MouseCellAtlas_Embryo_all.txt", quote=F, row.names=F, sep="\t")

### only neonatal samples
cell_label <- sort(unique(paste(samples$ts[samples$type=="Neonatal"], gsub(" ","_",samples$cell_type[samples$type=="Neonatal"]), sep=":")))
ts_cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
colnames(ts_cell_log_umi) <- cell_label
ts <- unique(samples$ts[samples$type=="Neonatal"])
for(t in ts){
  ct <- unique(samples$cell_type[samples$ts==t])
  for(c in ct){
    if(length(which(samples$ts==t & samples$cell_type==c))==1){
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]]
    }else{
      ts_cell_log_umi[,paste(t,gsub(" ","_",c), sep=":")] = apply(log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]], 1, mean)
    }
  }
}
ts_cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], ts_cell_log_umi)
ts_cell_log_umi$Average <- apply(ts_cell_log_umi[,2:ncol(ts_cell_log_umi)], 1, mean)

write.table(ts_cell_log_umi, "processed_data/MouseCellAtlas_Neonatal_all.txt", quote=F, row.names=F, sep="\t")

### per tissue
ts <- unique(samples$ts)
for(t in ts){
  ts2 <- unique(samples$Tissue[samples$ts==t])
  if(length(ts2)==1){
    ct <- sort(unique(samples$cell_type[samples$ts==t]))
    if(length(ct)==1){next}
    cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
    colnames(cell_log_umi) <- gsub(" ", "_", ct)
    for(c in ct){
      if(length(which(samples$ts==t & samples$cell_type==c))==1){
        cell_log_umi[,gsub(" ","_",c)] = log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]]
      }else{
        cell_log_umi[,gsub(" ", "_", c)] = apply(log_umi[,samples$Cell.name[samples$ts==t & samples$cell_type==c]], 1, mean)
      }
    }
  }else{
    cell_label <- sort(unique(paste(samples$Tissue[samples$ts==t], gsub(" ", "_",samples$cell_type[samples$ts==t]), sep=":")))
    cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
    colnames(cell_log_umi) <- cell_label
    for(t2 in ts2){
      ct <- sort(unique(samples$cell_type[samples$Tissue==t2]))
      for(c in ct){
        if(length(which(samples$Tissue==t2 & samples$cell_type==c))==1){
          cell_log_umi[,paste(t2,gsub(" ","_",c), sep=":")] = log_umi[,samples$Cell.name[samples$Tissue==t2 & samples$cell_type==c]]
        }else{
          cell_log_umi[,paste(t2,gsub(" ","_",c), sep=":")] = apply(log_umi[,samples$Cell.name[samples$Tissue==t2 & samples$cell_type==c]], 1, mean)
        }
      }
    }
  }
  cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
  cell_log_umi$Average <- apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)

  write.table(cell_log_umi, paste0("processed_data/MouseCellAtlas_", t, ".txt"), quote=F, row.names=F, sep="\t")
}
