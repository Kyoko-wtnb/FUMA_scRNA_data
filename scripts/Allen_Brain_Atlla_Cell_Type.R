##### Allen Brain Atlas Cell Types data #####
### data was downloaded directly from Allen Brain Atlas web-site
### (http://celltypes.brain-map.org/rnaseq)
### This script contains 2 human (MTG and LGN) and
### 3 mouse (VIPs, ALM and LGd) datasets.
### Gene expression was computed as the sum of exonic and intronic expression per gene.

library(data.table)
library(Hmisc)

##### Human MTG #####
### Human Middle temporal gyrus
### Two levels of cell type annotations were used
### level 1: "class" column
### level 2: "cluster" column
### Cells with label "no class" were excluded.

### expression and sample data
samples <- fread("Allen/Human/MTG/human_MTG_2018-06-14_samples-columns.csv", data.table=F)
genes <- fread("Allen/Human/MTG/human_MTG_2018-06-14_genes-rows.csv", data.table=F)
exon <- fread("Allen/Human/MTG/human_MTG_2018-06-14_exon-matrix.csv", data.table=F)
exon <- as.matrix(exon[,2:ncol(exon)])
intron <- fread("Allen/Human/MTG/human_MTG_2018-06-14_intron-matrix.csv", data.table=F)
intron <- as.matrix(intron[,2:ncol(intron)])
count <- exon+intron
rm(exon, intron)
rownames(count) <- genes$gene
samples <- samples[samples$cluster!="no class",]
count <- count[,samples$sample_name]

### map genes to ENSG
ENSG <- fread("BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(count), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
genes$hs.ensg[is.na(genes$hs.ensg)] <- sapply(genes$symbol[is.na(genes$hs.ensg)], function(x){
  n <- grepl(paste0("|",x,"|"), ENSG$alias)
  if(length(n)==1){ENSG$ensembl_gene_id[n]}
  else{NA}
})
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
count <- count[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

### convert to cpm
lib_size <- apply(count, 2, sum)/1000000
cpm <- t(apply(count, 1, function(x){x/lib_size})) # high RAM
log_cpm <- log2(cpm+1)
rm(count, cpm)

### level 1
ct <- sort(unique(samples$class))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$class==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Human_MTG_level1.txt"), quote=F, row.names=F, sep="\t")

### level 2
ct <- sort(unique(samples$cluster))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$cluster==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Human_MTG_level2.txt"), quote=F, row.names=F, sep="\t")

##### Human LGN #####
### Human Lateral geniculate
### Two levels of cell type annotations were used
### level 1: "class" column
### level 2: "cluster" column
### Cells with label "no class" were excluded.

### expression and sample data
samples <- fread("Allen/Human/LGN/human_LGN_2018-06-14_samples-columns.csv", data.table=F)
genes <- fread("Allen/Human/LGN/human_LGN_2018-06-14_genes-rows.csv", data.table=F)
exon <- fread("Allen/Human/LGN/human_LGN_2018-06-14_exon-matrix.csv", data.table=F)
exon <- as.matrix(exon[,2:ncol(exon)])
intron <- fread("Allen/Human/LGN/human_LGN_2018-06-14_intron-matrix.csv", data.table=F)
intron <- as.matrix(intron[,2:ncol(intron)])
count <- exon+intron
rm(exon, intron)
rownames(count) <- genes$gene
samples <- samples[samples$cluster!="no class",]
count <- count[,samples$sample_name]

### map genes to ENSG
ENSG <- fread("BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(count), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
genes$hs.ensg[is.na(genes$hs.ensg)] <- sapply(genes$symbol[is.na(genes$hs.ensg)], function(x){
  n <- grepl(paste0("|",x,"|"), ENSG$alias)
  if(length(n)==1){ENSG$ensembl_gene_id[n]}
  else{NA}
})
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
count <- count[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

### convert to cpm
lib_size <- apply(count, 2, sum)/1000000
cpm <- t(apply(count, 1, function(x){x/lib_size}))
log_cpm <- log2(cpm+1)
rm(count, cpm)

### level1
ct <- sort(unique(samples$class))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$class==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Human_LGN_level1.txt"), quote=F, row.names=F, sep="\t")

### level 2
ct <- sort(unique(samples$cluster))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$cluster==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Human_LGN_level2.txt"), quote=F, row.names=F, sep="\t")

##### Mouse VIPs v2 #####
### Mouse Visual cortex
### Three levels of cell type annotations were used
### level 1: "class" column
### level 2: "subclass" column
### level 3: "cluster" column
### Cells with label "No Class" or "Low Quality" were excluded.

### expression and sample data
samples <- fread("Allen/Mouse/VISp2/Mouse_VISp_2018-06-14_samples-columns.csv", data.table=F)
genes <- fread("Allen/Mouse/VISp2/Mouse_VISp_2018-06-14_genes-rows.csv", data.table=F)
exon <- fread("Allen/Mouse/VISp2/Mouse_VISp_2018-06-14_exon-matrix.csv", data.table=F)
exon <- as.matrix(exon[,2:ncol(exon)])
intron <- fread("Allen/Mouse/VISp2/Mouse_VISp_2018-06-14_intron-matrix.csv", data.table=F)
intron <- as.matrix(intron[,2:ncol(intron)])
count <- exon+intron
rm(exon, intron)
rownames(count) <- genes$gene_symbol
samples <- samples[!samples$class%in%c("No Class", "Low Quality"),]
count <- count[,samples$sample_name]

### map genes to hs ENSG ID
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(count), stringsAsFactors = F)
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
count <- count[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

### convert to cpm
lib_size <- apply(count, 2, sum)/1000000
cpm <- t(apply(count, 1, function(x){x/lib_size})) # high RAM
log_cpm <- log2(cpm+1)
rm(count, cpm)

### level1
ct <- sort(unique(samples$class))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$class==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_VISp2_level1.txt"), quote=F, row.names=F, sep="\t")

### level 2
ct <- sort(unique(samples$subclass))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$subclass==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_VISp2_level2.txt"), quote=F, row.names=F, sep="\t")

### level 3
ct <- sort(unique(samples$cluster))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  if(length(which(samples$cluster==c))==1){
    cell_log_cpm[,gsub(" ", "_", c)] <- log_cpm[,samples$sample_name[samples$cluster==c]]
  }else{
    cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$cluster==c]], 1, mean)
  }
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_VISp2_level3.txt"), quote=F, row.names=F, sep="\t")

##### Mouse ALM v2 #####
### Mouse Anterior lateral motor cortex
### level 1: "class" column
### level 2: "subclass" column
### level 3: "cluster" column
### Cells with label "No Class" or "Low Quality" were excluded.

### expression and sample data
samples <- fread("Allen/Mouse/ALM2/Mouse_ALM_2018-06-14_samples-columns.csv", data.table=F)
genes <- fread("Allen/Mouse/ALM2/Mouse_ALM_2018-06-14_genes-rows.csv", data.table=F)
exon <- fread("Allen/Mouse/ALM2/Mouse_ALM_2018-06-14_exon-matrix.csv", data.table=F)
exon <- as.matrix(exon[,2:ncol(exon)])
intron <- fread("Allen/Mouse/ALM2/Mouse_ALM_2018-06-14_intron-matrix.csv", data.table=F)
intron <- as.matrix(intron[,2:ncol(intron)])
count <- exon+intron
rm(exon, intron)
rownames(count) <- genes$gene_symbol
samples <- samples[!samples$class%in%c("No Class", "Low Quality"),]
count <- count[,samples$sample_name]

### map genes to hs ENSG ID
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(count), stringsAsFactors = F)
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
count <- count[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

### convert to cpm
lib_size <- apply(count, 2, sum)/1000000
cpm <- t(apply(count, 1, function(x){x/lib_size})) # high RAM
log_cpm <- log2(cpm+1)
rm(count, cpm)

### level1
ct <- sort(unique(samples$class))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$class==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_ALM2_level1.txt"), quote=F, row.names=F, sep="\t")

### level 2
ct <- sort(unique(samples$subclass))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$subclass==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_ALM2_level2.txt"), quote=F, row.names=F, sep="\t")

### level 3
ct <- sort(unique(samples$cluster))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  if(length(which(samples$cluster==c))==1){
    cell_log_cpm[,gsub(" ", "_", c)] <- log_cpm[,samples$sample_name[samples$cluster==c]]
  }else{
    cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$cluster==c]], 1, mean)
  }
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_ALM2_level3.txt"), quote=F, row.names=F, sep="\t")

##### Mouse LGd v2 #####
### Mouse Dorsolateral geniculate
### level 1: "class" column
### level 2: "subclass" column
### level 3: "cluster" column
### Cells with label "No Class" or "Low Quality" were excluded.

### expression and sample data
samples <- fread("Allen/Mouse/LGd2/Mouse_LGd_2018-06-14_samples-columns.csv", data.table=F)
genes <- fread("Allen/Mouse/LGd2/Mouse_LGd_2018-06-14_genes-rows.csv", data.table=F)
exon <- fread("Allen/Mouse/LGd2/Mouse_LGd_2018-06-14_exon-matrix.csv", data.table=F)
exon <- as.matrix(exon[,2:ncol(exon)])
intron <- fread("Allen/Mouse/LGd2/Mouse_LGd_2018-06-14_intron-matrix.csv", data.table=F)
intron <- as.matrix(intron[,2:ncol(intron)])
count <- exon+intron
rm(exon, intron)
rownames(count) <- genes$gene_symbol
samples <- samples[samples$class!="Outlier",]
count <- count[,samples$sample_name]

### map genes to hs ENSG ID
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(count), stringsAsFactors = F)
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
count <- count[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

### convert to cpm
lib_size <- apply(count, 2, sum)/1000000
cpm <- t(apply(count, 1, function(x){x/lib_size}))
log_cpm <- log2(cpm+1)
rm(count, cpm)

### level1
ct <- sort(unique(samples$class))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$class==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_LGd2_level1.txt"), quote=F, row.names=F, sep="\t")

### level 2
ct <- sort(unique(samples$subclass))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$subclass==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_LGd2_level2.txt"), quote=F, row.names=F, sep="\t")

### level 3
ct <- sort(unique(samples$cluster))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub(" ", "_", ct)
for(c in ct){
  if(length(which(samples$cluster==c))==1){
    cell_log_cpm[,gsub(" ", "_", c)] <- log_cpm[,samples$sample_name[samples$cluster==c]]
  }else{
    cell_log_cpm[,gsub(" ", "_", c)] <- apply(log_cpm[,samples$sample_name[samples$cluster==c]], 1, mean)
  }
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Allen_Mouse_LGd2_level3.txt"), quote=F, row.names=F, sep="\t")
