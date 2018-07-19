##### Other data sets from GEO #####
### Unless specified, data was downloaded from corresponding GEO repository
### Detailed were provided for each data set separately

library(data.table)

##### GSE87544 #####
### original study: chen et al. Cell Rep. 2017 (PMID: 28355573)
### expression data was downloaded from GSE87544
### Only cells with "Normal" consition were used
### To be consistent with teh original study, cells were filtered on those with >2000 expressed genes (non-zero expression)
### Cells with label "zothers" were excluded
gse <- "GSE87544"

### samples
samples <- fread(input="gzip -cd GEO/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv.gz", data.table=F, header=T)
samples <- samples[,2:3]
colnames(samples) <- c("cell_id", "cell_type")
notes <- paste0("Number of samples: ", nrow(samples))
samples <- samples[grepl("Normal", samples$cell_id),]
samples <- samples[samples$cell_type!="zothers",]

### expression data
umi <- fread(input="gzip -cd GEO/GSE87544_Merged_17samples_14437cells_count.txt.gz", data.table=F)
genes <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- genes
umi <- umi[,samples$cell_id]
n_genes <- apply(umi, 2, function(x){length(which(x>0))})
umi <- umi[,n_genes>2000]
samples <- samples[samples$cell_id %in% colnames(umi),]

### map symbols to hs ENSG
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

### process
log_umi <- log(umi+1)
ct <- sort(unique(samples$cell_type))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  if(length(which(samples$cell_type==c))==1){
    cell_log_umi[,c] <- log_umi[,samples$cell_id[samples$cell_type==c]]
  }else{
    cell_log_umi[,c] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
  }
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/", gse, "_Mouse_Hypothalamus.txt"), quote=F, row.names=F, sep="\t")

##### GSE98816 #####
### original study: Vanlandewijck et al. Nature 2018 (PMID: 29443965)
### exxpression data was downloaded from GSE98816 (used raw count matrix)
### Cell type label was shaed by the authors
gse <- "GSE98816"

### samples
samples <- fread("GEO/GSE98816_celltype.txt", data.table=F)

### expresion
count <- fread(input="gzip -cd GEO/GSE98816_Brain_samples_raw_read_counts_matrix.txt.gz", data.table=F)
g <- count[,1]
count <- as.matrix(count[,2:ncol(count)])
rownames(count) <- g

### map symbols to hs ENSG ID
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

### process
lib_size <- apply(count, 2, sum)/1000000
log_cpm <- log2(t(apply(count, 1, function(x){x/lib_size}))+1)
ct <- sort(unique(samples$cell_type))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/", gse, "_Mouse_Brain_Vascular.txt"), quote=F, row.names=F, sep="\t")

##### GSE99235 #####
### original study: Vanlandewijck et al. Nature 2018 (PMID: 29443965)
### exxpression data was downloaded from GSE99235 (used raw count matrix)
### Cell type label was shaed by the authors
gse <- "GSE99235"

### samples
samples <- fread("GEO/GSE99235_celltype.txt", data.table=F)

### expresion
count <- fread(input="gzip -cd GEO/GSE99235_Lung_samples_raw_reads_counts_matrix.txt.gz", data.table=F)
g <- count[,1]
count <- as.matrix(count[,2:ncol(count)])
rownames(count) <- g

### map symbols to hs ENSG ID
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

### process
lib_size <- apply(count, 2, sum)/1000000
log_cpm <- log2(t(apply(count, 1, function(x){x/lib_size}))+1)
ct <- sort(unique(samples$cell_type))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/", gse, "_Mouse_Lung_Vascular.txt"), quote=F, row.names=F, sep="\t")

##### GSE81547 #####
### original study: Enge et al. Cell 2017 (PMID: 28965763)
### expression data was downloaded from GSE81547
### metadata was extracted from family soft file
### raw expression per individual cell were concatanated prior to this process on R
### see README for details
gse <- "GSE81547"

### samples
samples <- fread("GEO/GSE81547_metadata.txt", data.table=F)
## reprace unsure to "PP" to be consistent with the original study
samples$cell_type[samples$cell_type=="unsure"] <- "PP"

### expression
count <- fread(input="gzip -cd GEO/GSE81547_count.txt.gz", data.table=F)
g <- count[,1]
count <- as.matrix(count[,2:ncol(count)])
rownames(count) <- g

### map genes to ENSG
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(count), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
count <- count[genes$symbol[!is.na(genes$hs.ensg)],]

### process
lib_size <- apply(count, 2, sum)/1000000
log_cpm <- log2(t(apply(count, 1, function(x){x/lib_size}))+1)
ct <- sort(unique(samples$cell_type))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/", gse, "_Human_Pancreas.txt"), quote=F, row.names=F, sep="\t")

##### GSE104276 #####
### original study: Zhong et al. Nature 2018 (PMID: 29539641)
### expression data was downloaded from GSE104276 (TPM_NOERCC file)
### Cell type was extracted from excel sheet (GSE104276_readme_sample_barcode.xlsx, sheet named "SampleInfo")
### We created the followoing two data sets,
### 1) per cell type average expression across developmental stage
### 2) per cell type per developmental stage average expression
gse <- "GSE104276"

### samples
samples <- fread("GEO/GSE104276_cell_type.txt", data.table=F)
colnames(samples)[3] <- "age"

### expreission
tpm <- fread(input="gzip -cd GEO/GSE104276_all_pfc_2394_UMI_TPM_NOERCC.csv.gz", data.table=F)
g <- tpm[,1]
colnames(tpm)[2:ncol(tpm)] <- colnames(tpm)[1:(ncol(tpm)-1)]
tpm <- as.matrix(tpm[2:ncol(tpm)])
rownames(tpm) <- g

### map genes to ENSG
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(tpm), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
tpm <- tpm[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_tpm <- log2(tpm+1)
## age mixed
ct <- sort(unique(samples$cell_type))
cell_log_tpm <- matrix(nrow=nrow(log_tpm), ncol=length(ct))
colnames(cell_log_tpm) <- gsub(" ", "_", ct)
for(c in ct){
  cell_log_tpm[,gsub(" ", "_", c)] <- apply(log_tpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_tpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_tpm), genes$symbol)], cell_log_tpm)
cell_log_tpm$Average = apply(cell_log_tpm[,2:ncol(cell_log_tpm)], 1, mean)
write.table(cell_log_tpm, paste0("processed_data/", gse, "_Human_Prefrontal_cortex_all_ages.txt"), quote=F, row.names=F, sep="\t")
## per age
ct <- unique(paste(samples$age, gsub(" ", "_", samples$cell_type), sep="_"))
cell_log_tpm <- matrix(nrow=nrow(log_tpm), ncol=length(ct))
colnames(cell_log_tpm) <- ct
for(c in ct){
  if(length(which(paste(samples$age, gsub(" ", "_", samples$cell_type), sep="_")==c))==1){
    cell_log_tpm[,c] <- log_tpm[,samples$cell_id[paste(samples$age, gsub(" ", "_", samples$cell_type), sep="_")==c]]
  }else{
    cell_log_tpm[,c] <- apply(log_tpm[,samples$cell_id[paste(samples$age, gsub(" ", "_", samples$cell_type), sep="_")==c]], 1, mean)
  }
}
cell_log_tpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_tpm), genes$symbol)], cell_log_tpm)
cell_log_tpm$Average = apply(cell_log_tpm[,2:ncol(cell_log_tpm)], 1, mean)
write.table(cell_log_tpm, paste0("processed_data/", gse, "_Human_Prefrontal_cortex_per_ages.txt"), quote=F, row.names=F, sep="\t")

##### GSE82187 #####
### original study: Gokce et al. Cell Repo. 2016 (PMID: 27425622)
### expression data was downloaded from GSE82187
### we only used microfluid data since FACS data was limited to neurons
gse <- "GSE82187"

### expression and metadata
cpm <- fread("gzip -cd GEO/GSE82187_cast_all_forGEO.csv.gz", data.table=F, sep=",", header=T)
samples <- cpm[,2:5]
colnames(samples)[1:2] <- c("cell_id", "cell_type")
cpm <- t(cpm[,6:ncol(cpm)])
colnames(cpm) <- samples$cell_id

#### Mic scRNA
samples <- samples[samples$protocol=="Mic-scRNA-Seq",]
log_cpm <- cpm[,samples$cell_id[samples$protocol=="Mic-scRNA-Seq"]]

### map symbols to hs ENSG ID
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(log_cpm), stringsAsFactors = F)
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
log_cpm <- log_cpm[genes$symbol[!is.na(genes$hs.ensg)],]

### process
ct <- sort(unique(samples$cell_type))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- gsub("-", "_", ct)
for(c in ct){
  cell_log_cpm[,gsub("-", "_", c)] <- apply(log_cpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/", gse, "_Mouse_Striatum.txt"), quote=F, row.names=F, sep="\t")


##### GSE89232 #####
### oeiginal study: Breton et al. J. Exp. Med. 2016 (PMID: 27864467)
### expression data was downloaded from GSE89232
### metadata was extracted from family soft file
### see README for detils
### Cell type labels were manually aligned with the orignal study
gse <- "GSE89232"

### samples
samples <- fread("GEO/GSE89232_celltype.txt", data.table=F)
samples$cell_id <- sub("Blood pre-cDC", "blood_predc", samples$cell_id)
samples$cell_id <- sub("CB pre-cDC", "cb_predc", samples$cell_id)
samples$cell_id <- sub("CD1c\\+ cDC", "blood_predc", samples$cell_id)
samples$cell_id <- sub("CD141\\+ cDC", "blood_predc", samples$cell_id)

### expression
tpm <- fread(input="gzip -cd GEO/GSE89232_expMatrix.txt.gz", data.table=F, header=T, fill=T)
colnames(tpm)[2:ncol(tpm)] <- colnames(tpm)[1:(ncol(tpm)-1)]
g <- tpm[,1]
tpm <- as.matrix(tpm[,2:ncol(tpm)])
rownames(tpm) <- g

### map genes to ENSG
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(tpm), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
tpm <- tpm[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_tpm <- tpm
ct <- sort(unique(samples$cell_type))
cell_log_tpm <- matrix(nrow=nrow(log_tpm), ncol=length(ct))
colnames(cell_log_tpm) <- gsub(" ", "_", sub("\\+", "", sub("-", "_", ct)))
for(c in ct){
  cell_log_tpm[,gsub(" ", "_", sub("\\+", "", sub("-", "_", c)))] <- apply(log_tpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_tpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_tpm), genes$symbol)], cell_log_tpm)
cell_log_tpm$Average = apply(cell_log_tpm[,2:ncol(cell_log_tpm)], 1, mean)
write.table(cell_log_tpm, paste0("processed_data/", gse, "_Human_Blood.txt"), quote=F, row.names=F, sep="\t")

##### GSE100597 #####
### original study: Mohammed, et al. Cell Rep. 2017 (PMID: 28768204)
### expression data was downloaded from GSE100597
### only developmental stage was used as the cell label
### E6.75 was replaced with E6.5 to be consistent with the original stidu
gse <- "GSE100597"

### expression and samples
count <- fread(input="gzip -cd GEO/GSE100597_count_table_QC_filtered.txt.gz", data.table=F)
g <- count[,1]
g <- as.character(unlist(sapply(g, function(x){unlist(strsplit(x, "_"))[1]})))
count <- as.matrix(count[,2:ncol(count)])
rownames(count) <- g
samples <- data.frame(cell_id=colnames(count), cell_type=sub("(E\\d+\\.\\d+)_.+", "\\1", colnames(count)), stringsAsFactors = F)
samples$cell_type[samples$cell_type=="E6.75"] <- "E6.5"

### map symbols to hs ENSG ID
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

### process
lib_size <- apply(count, 2, sum)/1000000
log_cpm <- log2(t(apply(count, 1, function(x){x/lib_size}))+1)
ct <- sort(unique(samples$cell_type))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/", gse, "_Mouse_Embryo.txt"), quote=F, row.names=F, sep="\t")

##### GSE93374 #####
### original study: campbell et al. Nat. Neurosci. 2017 (PMID: 28166221)
### expression data was downloaded from GSE93374
### Batch corrected normalized data was used (GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt.gz)
### Level 1, level 2 and clusters for neurons were processed separately
### In total 3 data sets were created
### Cells with label "miss" were excluded
### For the 3rd data set, only neuronal cells were used
gse <- "GSE93374"

### samples
samples <- fread(input="gzip -cd GEO/GSE93374_cell_metadata.txt.gz", data.table=F)
samples[1,] <- gsub(" ", "_",sub("^\\d+\\.", "", samples[1,]))
colnames(samples) <- samples[1,]
samples <- samples[2:nrow(samples),]
samples <- samples[samples$clust_all!="miss",]
cell_label <- unique(c(samples$All_Cell_Clusters,samples$All_Cell_Subclusters, samples$Neuron_Subclusters))
cell_label <- cell_label[cell_label!="" & cell_label!="miss"]
cell_label <- data.frame(clst=sub("(.+)\\..+", "\\1", cell_label), label=sub(".+\\.(.+)", "\\1", cell_label))
cell_label$label <- gsub("/", "_", cell_label$label)

### expression
log_umi <- fread(input="gzip -cd GEO/GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt.gz", data.table=F)
g <- log_umi[,1]
log_umi <- as.matrix(log_umi[,2:ncol(log_umi)])
rownames(log_umi) <- g

### map symbols to hs ENSG
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(log_umi), stringsAsFactors = F)
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
log_umi <- log_umi[genes$symbol[!is.na(genes$hs.ensg)],]

### process
### level1
ct <- sort(unique(samples$clust_all))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$ID[samples$clust_all==c]], 1, mean)
}
colnames(cell_log_umi) <- cell_label$label[match(ct, cell_label$clst)]
cell_log_umi <- cell_log_umi[,sort(colnames(cell_log_umi))]
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/", gse, "_Mouse_Arc_ME_level1.txt"), quote=F, row.names=F, sep="\t")

### level2
ct <- sort(unique(samples$clust_all_micro))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$ID[samples$clust_all_micro==c]], 1, mean)
}
colnames(cell_log_umi) <- cell_label$label[match(ct, cell_label$clst)]
cell_log_umi <- cell_log_umi[,sort(colnames(cell_log_umi))]
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/", gse, "_Mouse_Arc_ME_level2.txt"), quote=F, row.names=F, sep="\t")

### neurons
ct <- sort(unique(samples$clust_neurons[samples$clust_neurons!="miss"]))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$ID[samples$clust_neurons==c]], 1, mean)
}
colnames(cell_log_umi) <- cell_label$label[match(ct, cell_label$clst)]
cell_log_umi <- cell_log_umi[,sort(colnames(cell_log_umi))]
cell_log_umi <- cell_log_umi[,!grepl("unassigned", colnames(cell_log_umi))]
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/", gse, "_Mouse_Arc_ME_neurons.txt"), quote=F, row.names=F, sep="\t")

##### GSE92332 #####
### original study: Haber et al. Nature 2017 (PMID: 29144463)
### expression data was downloaded from GSE923323
### GSE92332_AtlasFullLength_TPM.txt.gz for SMARTSeq data set
### GSE92332_atlas_UMIcounts.txt.gz for droplet data set
### Each data set was processed separately
gse <- "GSE92332"

#### SMARTseq data set
### expression and samples
tpm <- fread(input="gzip -cd GEO/GSE92332_AtlasFullLength_TPM.txt.gz", data.table=F, fill=T)
colnames(tpm)[2:ncol(tpm)] <- colnames(tpm)[1:(ncol(tpm)-1)]
g <- tpm[,1]
tpm <- as.matrix(tpm[,2:ncol(tpm)])
rownames(tpm) <- g
samples <- data.frame(cell_id=colnames(tpm), stringsAsFactors = F)
samples$cell_type <- sapply(samples$cell_id, function(x){tmp<-unlist(strsplit(x, "_")); tmp[length(tmp)]})

### map symbols to hs ENSG
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(tpm), stringsAsFactors = F)
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
tpm <- tpm[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_tpm <- log2(tpm+1)
ct <- sort(unique(samples$cell_type))
cell_log_tpm <- matrix(nrow=nrow(log_tpm), ncol=length(ct))
colnames(cell_log_tpm) <- gsub("-", "_", ct)
for(c in ct){
  cell_log_tpm[,gsub("-", "_", c)] <- apply(log_tpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_tpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_tpm), genes$symbol)], cell_log_tpm)
cell_log_tpm$Average = apply(cell_log_tpm[,2:ncol(cell_log_tpm)], 1, mean)
write.table(cell_log_tpm, paste0("processed_data/", gse, "_Mouse_Epithelium_SMARTseq.txt"), quote=F, row.names=F, sep="\t")

##### droplet dat aset
### expression and samples
umi <- fread(input="gzip -cd GEO/GSE92332_atlas_UMIcounts.txt.gz", data.table=F, fill=T)
colnames(umi)[2:ncol(umi)] <- colnames(umi)[1:(ncol(umi)-1)]
g <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- g
samples <- data.frame(cell_id=colnames(umi), stringsAsFactors = F)
samples$cell_type <- sapply(samples$cell_id, function(x){tmp<-unlist(strsplit(x, "_")); tmp[length(tmp)]})
samples$cell_type <- gsub("\\.", "_", samples$cell_type)

### map symbols to hs ENSG
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

### process
log_umi <- log2(umi+1)
ct <- sort(unique(samples$cell_type))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- gsub("-", "_", ct)
for(c in ct){
  cell_log_umi[,gsub("-", "_", c)] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/", gse, "_Mouse_Epithelium_droplet.txt"), quote=F, row.names=F, sep="\t")

##### GSE89164 #####
### original study: Alle et al. MBC Biol. 2017 (PMID: 28526029)
### expression data was downloaded from GSE89164 (raw data)
### Combiend count matrices for two mouse replicates
### Cell type label was manually assigned to the cluster index based on the original study
gse <- "GSE89164"

### samples
samples <- fread(input="gzip -cd GEO/GSE89164_clusters_mm.txt.gz", data.table=F)
colnames(samples) <- c("id", "clst")
samples$cell_id <- sub("(.+)_.+", "\\1", samples$id)
samples$clst <- as.character(samples$clst)
cell_label <- data.frame(clst=as.character(0:11), label=c("Excitatory_neurons", "Granule_cells", "Oligodendrocyte_progenitors", "Inhibitory_neurons", "Astrocytes",
                                                          "Cerebellar_neurons", "Oligodendrotype_myelinating", "Microglia", "Astrocyte_progenitors",
                                                          "Endothelial_progenitors", "Endothelial_cells", "SC"))

### expression
umi <- fread(input="gzip -cd GEO/GSE89164_raw/GSM2518784_mm_rep1_dge.txt.gz", data.table=F)
g <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- g
tmp <- fread(input="gzip -cd GEO/GSE89164_raw/GSM2518785_mm_rep2_dge.txt.gz", data.table=F)
g <- tmp[,1]
tmp <- as.matrix(tmp[,2:ncol(tmp)])
rownames(tmp) <- g
umi <- cbind(umi, tmp[match(rownames(umi),rownames(tmp)),])
check <- apply(umi, 1, function(x){length(which(is.na(x)))})
umi <- umi[check==0,]
umi <- umi[,samples$cell_id]

### map symbols to hs ENSG
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

### process
log_umi <- log2(umi+1)
ct <- sort(unique(samples$clst))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$cell_id[samples$clst==c]], 1, mean)
}
colnames(cell_log_umi) <- cell_label$label[match(colnames(cell_log_umi), cell_label$clst)]
cell_log_umi <- cell_log_umi[,sort(colnames(cell_log_umi))]
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/", gse, "_Mouse_Hindbrain.txt"), quote=F, row.names=F, sep="\t")

##### GSE67835 #####
### original study: Darmanis et al. Proc. Natl. Acad. Sci. USA 2015 (PMID: 26060301)
### expression data was downloaded from GSE67835 (raw data)
### metadata was extracted from family soft file
### Per cell count was concatenated into a single matrix prior to this process in R
### see README for details
### We created two data set, with and without fetal samples
gse <- "GSE67835"

### samples
samples <- fread("GEO/GSE67835_metadata.txt", data.table=F)

### expression
count <- fread(input="gzip -cd GEO/GSE67835_count.txt.gz", data.table=F)
g <- count[,1]
count <- as.matrix(count[,2:ncol(count)])
rownames(count) <- g

### map genes to ENSG
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(count), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
count <- count[genes$symbol[!is.na(genes$hs.ensg)],]

### process
lib_size <- apply(count, 2, sum)/1000000
log_cpm <- log2(t(apply(count, 1, function(x){x/lib_size}))+1)
## including fetal samples
ct <- sort(unique(samples$cell_type))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/", gse, "_Human_Cortex.txt"), quote=F, row.names=F, sep="\t")

## without fetal samples
samples <- samples[!grepl("fetal", samples$cell_type),]
ct <- sort(unique(samples$cell_type))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/", gse, "_Human_Cortex_woFetal.txt"), quote=F, row.names=F, sep="\t")

##### GSE106678 #####
### original study: Hu et al. Mol. Cell. 2017 (PMID: 29220646)
### expression data was downloaded from Broadinstitute Single Cell Portal at
### https://portals.broadinstitute.org/single_cell/study/snucdrop-seq-dissecting-cell-type-composition-and-activity-dependent-transcriptional-state-in-mammalian-brains-by-massively-parallel-single-nucleus-rna-seq#study-download
gse <- "GSE106678"

### samples
samples <- fread("GEO/GSE106678/metadata.txt", data.table=F)
samples <- samples[2:nrow(samples),]
samples$NAME <- sub("_y$", ".y", samples$NAME)
samples$NAME <- sub("_x$", ".x", samples$NAME)

### expression
umi <- fread(input="gzip -cd GEO/GSE106678/expression.txt.gz", data.table=F)
g <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- g

### map symbols to hs ENSG
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

### process
log_umi <- log2(umi+1)
ct <- sort(unique(samples$Cluster))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$NAME[samples$Cluster==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/", gse, "_Mouse_Cortex.txt"), quote=F, row.names=F, sep="\t")
