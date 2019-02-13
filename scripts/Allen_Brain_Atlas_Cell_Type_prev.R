##### Allen Brain Atlas Cell Types data #####
### data was downloaded from Broadinstitute Single Cell Portal at https://portals.broadinstitute.org/single_cell
### (since cell type label was not available from the Allen Brain Atlas website)
### This script refers to previously released 3 mouse scRNA-seq datasets

library(data.table)

##### VISp (GSE71585) #####
### samples (cell info)
samples <- fread("Allen/VISp/clusters.txt", data.table=F)
samples <- samples[2:nrow(samples),]
colnames(samples) <- c("cell_id", "level1", "level2")

### expresion
rpkm <- fread(input="Allen/VISp/expression_matrix.txt", data.table=F)
g <- rpkm[,1]
rpkm <- as.matrix(rpkm[,2:ncol(rpkm)])
rownames(rpkm) <- g

### map genes to hs ENSG
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(rpkm), stringsAsFactors = F)
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
rpkm <- rpkm[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_rpkm <- log2(rpkm+1)
## level1
ct <- sort(unique(samples$level1))
cell_log_rpkm <- matrix(nrow=nrow(log_rpkm), ncol=length(ct))
colnames(cell_log_rpkm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_rpkm[,gsub(" ", "_", c)] <- apply(log_rpkm[,samples$cell_id[samples$level1==c]], 1, mean)
}
cell_log_rpkm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_rpkm), genes$symbol)], cell_log_rpkm)
cell_log_rpkm$Average = apply(cell_log_rpkm[,2:ncol(cell_log_rpkm)], 1, mean)
write.table(cell_log_rpkm, paste0("processed_data/Allen_Mouse_VISp_level1.txt"), quote=F, row.names=F, sep="\t")

## level 2
ct <- sort(unique(samples$level2))
cell_log_rpkm <- matrix(nrow=nrow(log_rpkm), ncol=length(ct))
colnames(cell_log_rpkm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_rpkm[,gsub(" ", "_", c)] <- apply(log_rpkm[,samples$cell_id[samples$level2==c]], 1, mean)
}
cell_log_rpkm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_rpkm), genes$symbol)], cell_log_rpkm)
cell_log_rpkm$Average = apply(cell_log_rpkm[,2:ncol(cell_log_rpkm)], 1, mean)
write.table(cell_log_rpkm, paste0("processed_data/Allen_Mouse_VISp_level2.txt"), quote=F, row.names=F, sep="\t")

##### ALM #####
### samples
samples <- fread("Allen/ALM/clusters.txt", data.table=F)
samples <- samples[2:nrow(samples),]
colnames(samples) <- c("cell_id", "level1", "level2")
samples$level2 <- sub("  ", " ", samples$level2)

### expression
rpkm <- fread(input="Allen/ALM/expression.txt", data.table=F)
g <- rpkm[,1]
rpkm <- as.matrix(rpkm[,2:ncol(rpkm)])
rownames(rpkm) <- g

### map genes to hs ENSG
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(rpkm), stringsAsFactors = F)
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
rpkm <- rpkm[genes$symbol[!is.na(genes$hs.ensg)],]
cell_zero_lt_80 <- cell_zero_lt_80[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_rpkm <- log2(rpkm+1)
## level1
ct <- sort(unique(samples$level1))
cell_log_rpkm <- matrix(nrow=nrow(log_rpkm), ncol=length(ct))
colnames(cell_log_rpkm) <- gsub("-","_",gsub(" ", "_", ct))
for(c in ct){
  if(length(which(samples$level1==c))==1){
    cell_log_rpkm[,gsub("-","_",gsub(" ", "_", c))] <- log_rpkm[,samples$cell_id[samples$level1==c]]
  }else{
    cell_log_rpkm[,gsub("-","_",gsub(" ", "_", c))] <- apply(log_rpkm[,samples$cell_id[samples$level1==c]], 1, mean)
  }
}
cell_log_rpkm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_rpkm), genes$symbol)], cell_log_rpkm)
cell_log_rpkm$Average = apply(cell_log_rpkm[,2:ncol(cell_log_rpkm)], 1, mean)
write.table(cell_log_rpkm, paste0("processed_data/Allen_Mouse_ALM_level1.txt"), quote=F, row.names=F, sep="\t")

## level 2
ct <- sort(unique(samples$level2))
cell_log_rpkm <- matrix(nrow=nrow(log_rpkm), ncol=length(ct))
colnames(cell_log_rpkm) <- gsub("-","_",gsub(" ", "_", ct))
for(c in ct){
  if(length(which(samples$level2==c))==1){
    cell_log_rpkm[,gsub("-","_",gsub(" ", "_", c))] <- log_rpkm[,samples$cell_id[samples$level2==c]]
  }else{
    cell_log_rpkm[,gsub("-","_",gsub(" ", "_", c))] <- apply(log_rpkm[,samples$cell_id[samples$level2==c]], 1, mean)
  }
}
cell_log_rpkm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_rpkm), genes$symbol)], cell_log_rpkm)
cell_log_rpkm$Average = apply(cell_log_rpkm[,2:ncol(cell_log_rpkm)], 1, mean)
write.table(cell_log_rpkm, paste0("processed_data/Allen_Mouse_ALM_level2.txt"), quote=F, row.names=F, sep="\t")

##### LGd #####
### cells with label "Outlier" were excluded
### samples
samples <- fread("Allen/LGd/clusters.txt", data.table=F)
samples <- samples[2:nrow(samples),]
colnames(samples) <- c("cell_id", "level1", "level2")
samples <- samples[samples$level1 !="Outlier",]

### expression
rpkm <- fread(input="Allen/LGd/expression_matrix.txt", data.table=F)
g <- rpkm[,1]
rpkm <- as.matrix(rpkm[,2:ncol(rpkm)])
rownames(rpkm) <- g

### map genes to hs ENSG
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(rpkm), stringsAsFactors = F)
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
rpkm <- rpkm[genes$symbol[!is.na(genes$hs.ensg)],]
cell_zero_lt_80 <- cell_zero_lt_80[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_rpkm <- log2(rpkm+1)
## level1
ct <- sort(unique(samples$level1))
cell_log_rpkm <- matrix(nrow=nrow(log_rpkm), ncol=length(ct))
colnames(cell_log_rpkm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_rpkm[,gsub(" ", "_", c)] <- apply(log_rpkm[,samples$cell_id[samples$level1==c]], 1, mean)
}
cell_log_rpkm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_rpkm), genes$symbol)], cell_log_rpkm)
cell_log_rpkm$Average = apply(cell_log_rpkm[,2:ncol(cell_log_rpkm)], 1, mean)
write.table(cell_log_rpkm, paste0("processed_data/Allen_Mouse_LGd_level1.txt"), quote=F, row.names=F, sep="\t")

## level 2
ct <- sort(unique(samples$level2))
cell_log_rpkm <- matrix(nrow=nrow(log_rpkm), ncol=length(ct))
colnames(cell_log_rpkm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_rpkm[,gsub(" ", "_", c)] <- apply(log_rpkm[,samples$cell_id[samples$level2==c]], 1, mean)
}
cell_log_rpkm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_rpkm), genes$symbol)], cell_log_rpkm)
cell_log_rpkm$Average = apply(cell_log_rpkm[,2:ncol(cell_log_rpkm)], 1, mean)
write.table(cell_log_rpkm, paste0("processed_data/Allen_Mouse_LGd_level2.txt"), quote=F, row.names=F, sep="\t")
