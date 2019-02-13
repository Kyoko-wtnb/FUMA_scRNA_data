##### scRNA-seq data from Linnarsson's lab #####
### Unless specified, data was downloaded from corresponding GEO repository
### Detailed were provided for each data set separately

library(data.table)

##### GSE59739 #####
### original study: Usokin et al. Nat. Neurosci. 2015 (PMID: 25420068)
### expression data was downloaded from GSE59739
### metadata was extracted from S Table 1
### Only cells with label NF, NP, PEP or TH in level 1 were used to be consistent with the original study
### We created 3 data sets for each label level 1, 2, and 3
gse <- "GSE59739"
### samples
samples <- as.data.frame(t(fread("Linnarsson/GSE59739_metadata.txt", data.table=F, header=F)), stringsAsFactors=F)
colnames(samples) <- gsub(" ", "_", samples[1,])
samples <- samples[2:nrow(samples),]
rownames(samples) <- NULL

### expression
cpm <- fread(input="gzip -cd Linnarsson/GSE59739_DataTable.txt.gz", data.table=F)
cpm <- cpm[5:nrow(cpm),]
genes <- cpm[,1]
cpm <- as.matrix(apply(cpm[,2:ncol(cpm)], 2, as.numeric))
rownames(cpm) <- genes
notes <- c(notes, paste0("Number of genes: ", nrow(cpm)))

### select 622 Neuronal cells by filtering on Level 1 label
samples <- samples[samples$Level_1 %in% c("NF", "NP", "PEP", "TH"),]
cpm <- cpm[,samples$Sample_ID]

### map symbol to hs ENSG
load("genome/mm2hs.RData")
genes <- data.frame(symbol=rownames(cpm), stringsAsFactors = F)
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
cpm <- cpm[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_cpm <- log2(cpm+1)
cell_label <- list(NF="Neurofilament_containing", NP="Non_peptidergic_nocieptor", PEP="Peptidergic_nociceptor", TH="Tyrosine_hydroxylase_containing")
## level 1
ct <- sort(unique(samples$Level_1))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$Sample_ID[samples$Level_1==c]], 1, mean)
}
colnames(cell_log_cpm) <- as.character(cell_label[ct])
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Linnarsson_", gse, "_Mouse_Dorsal_root_ganglion_level1.txt"), quote=F, row.names=F, sep="\t")

## level 2
ct <- sort(unique(samples$Level_2))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$Sample_ID[samples$Level_2==c]], 1, mean)
}
for(n in names(cell_label)){
  ct <- sub(n, cell_label[n], ct)
}
colnames(cell_log_cpm) <- ct
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Linnarsson_", gse, "_Mouse_Dorsal_root_ganglion_level2.txt"), quote=F, row.names=F, sep="\t")

## level 3
ct <- sort(unique(samples$Level_3))
cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(ct))
colnames(cell_log_cpm) <- ct
for(c in ct){
  cell_log_cpm[,c] <- apply(log_cpm[,samples$Sample_ID[samples$Level_3==c]], 1, mean)
}
for(n in names(cell_label)){
  ct <- sub(n, cell_label[n], ct)
}
colnames(cell_log_cpm) <- ct
cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(log_cpm), genes$symbol)], cell_log_cpm)
cell_log_cpm$Average = apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
write.table(cell_log_cpm, paste0("processed_data/Linnarsson_", gse, "_Mouse_Dorsal_root_ganglion_level3.txt"), quote=F, row.names=F, sep="\t")

##### GSE60361 #####
### original data: Zeisel et al. Science 2015 (PMID: 25700174)
### data was downloaded from http://linnarssonlab.org/cortex/
### We created 2 data sets for each label level 1 and 2
### For level 2, cells with label "none" were excluded (189 cells)
gse <- "GSE60361"
### expression and samples
umi <- fread(input="gzip -cd Linnarsson/GSE60361_expression_mRNA_17-Aug-2014.txt.gz", data.table=F, header=F)
samples <- data.frame(t(umi[1:10,3:ncol(umi)]))
samples <- data.frame(apply(samples, 2, as.character), stringsAsFactors = F)
colnames(samples) <- gsub(" ", "_", umi[1:10, 2])
rownames(samples) <- NULL
genes <- umi[12:nrow(umi),1]
umi <- as.matrix(apply(umi[12:nrow(umi),3:ncol(umi)], 2, as.numeric))
colnames(umi) <- samples$cell_id
rownames(umi) <- genes

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
## level 1
ct <- sort(unique(samples$level1class))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- gsub(" ", "_", gsub("-","_",ct))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-","_",c))] <- apply(log_umi[,samples$cell_id[samples$level1class==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Cortex_Hippocampus_level1.txt"), quote=F, row.names=F, sep="\t")

## level 2
ct <- sort(sub("\\(none\\)", "none", gsub(" ", "_", gsub("-","_", unique(samples$level2class)))))
ct <- ct[ct != "none"]
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-","_", c))] <- apply(log_umi[,samples$cell_id[samples$level2class==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Cortex_Hippocampus_level2.txt"), quote=F, row.names=F, sep="\t")

##### GSE67602 #####
### original study: Joost et al. Cell Syst. 2016 (PMID: 27641957)
### expression data was downloaded from GSE67602
### cell type label was extracted from family soft file
### see README for details
gse <- "GSE67602"

### samples
samples <- fread("Linnarsson/GSE67602_celltype.txt", data.table=F)

### expression
umi <- fread(input="gzip -cd Linnarsson/GSE67602_Joost_et_al_expression.txt.gz", data.table=F)
umi <- umi[!grepl("^ERCC", umi[,1]),]
genes <- umi[,1]
umi <- as.matrix(apply(umi[,2:ncol(umi)], 2, as.numeric))
rownames(umi) <- genes

### map gsenes
load("genome/mm2hs.RData")
## map symbol to ENSG
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
cell_zero_lt_80 <- cell_zero_lt_80[genes$symbol[!is.na(genes$hs.ensg)],]
notes <- c(notes, paste0("Number of genes with unique hs ENSG: ", nrow(umi)))
save(list=c("umi", "cell_zero_lt_80"), file=paste0("Linnarsson/", gse, "_umi.RData"))

### process
log_umi <- log2(umi+1)
ct <- sort(unique(samples$cell_type))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- gsub(" ", "_", gsub("-","_",ct))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-","_",c))] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Skin_Epidermis.txt"), quote=F, row.names=F, sep="\t")

##### GSE74672 #####
### original study: Romanov et al. NAt. Neurosci. 2017 (PMID: 27991900)
### expression data was downloaded from GSE74672
### We created data sets for level 1 and 2 separately
### Level 2 label is only available for neurons
### For level 2, unclassified cells with label "uc" were excluded
gse <- "GSE74672"
### expression and samples
umi <- fread("gzip -cd Linnarsson/GSE74672_expressed_mols_with_classes.csv.gz", data.table=F)
samples <- data.frame(cellID = colnames(umi)[2:ncol(umi)], t(umi[1:8, 2:ncol(umi)]), stringsAsFactors = F)
colnames(samples)[2:ncol(samples)] <- gsub(" ", "_", umi[1:8,1])
rownames(samples) <- NULL
genes <- umi[12:nrow(umi),1]
umi <- as.matrix(apply(umi[12:nrow(umi), 2:ncol(umi)], 2, as.numeric))
rownames(umi) <- genes

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

## level 1
ct <- sort(unique(samples$level1_class))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$cellID[samples$level1_class==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Hypothalamus_level1.txt"), quote=F, row.names=F, sep="\t")

## level 2 neuron
ct <- sort(unique(samples$`level2_class_(neurons_only)`[samples$level1_class=="neurons" & samples$`level2_class_(neurons_only)`!="uc"]))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- gsub(" ", "_", gsub(",", "", gsub("\\(", "", gsub(")", "", ct))))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub(",", "", gsub("\\(", "", gsub(")", "", c))))] <- apply(log_umi[,samples$cellID[samples$`level2_class_(neurons_only)`==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
colnames(cell_log_umi) <- gsub("\\.", "", colnames(cell_log_umi))
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Hypothalamus_Neurons_level2.txt"), quote=F, row.names=F, sep="\t")

##### GSE75330 #####
### original study: Marques et al. Sience 2016 (PMID: 27284195)
### expression data was downloaded from GSE75330
### metadata was extracted from family soft file
### see README for details
gse <- "GSE75330"

### samples
samples <- fread("Linnarsson/GSE75330_metadata.txt", data.table=F)
samples$cell_id <- sub("T2-", "-", samples$cell_id)

### expression
umi <- fread(input="gzip -cd Linnarsson/GSE75330_Marques_et_al_mol_counts2.tab.gz", data.table=F)
genes <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- genes

### map symbols to hs ENSG
load("genome/mm2hs.RData")
## map symbol to ENSG
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
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Oligodendrocytes.txt"), quote=F, row.names=F, sep="\t")

##### GSE76381 #####
### original study: La Manno et al. Cell 2016 (PMID: 27716510)
### expression data was downloaded from GSE76381
### Only two data sets for mouse embryo and human embryo were used
### Cells with label "Unk" (unknown) were excluded
gse <- "GSE76381"

#### mouse data
umi <- fread(input="gzip -cd Linnarsson/GSE76381_MouseEmbryoMoleculeCounts.cef.txt.gz", data.table=F)
samples <- data.frame(t(umi[3:6,3:ncol(umi)]), stringsAsFactors = F)
colnames(samples) <- umi[3:6,2]
genes <- umi[8:nrow(umi),1]
umi <- as.matrix(apply(umi[8:nrow(umi), 3:ncol(umi)], 2, as.numeric))
colnames(umi) <- samples$Cell_ID
rownames(umi) <- genes
samples <- samples[samples$Cell_type!="mUnk",]
umi <- umi[,samples$Cell_ID]

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
ct <- sort(unique(samples$Cell_type))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$Cell_ID[samples$Cell_type==c]], 1, mean)
}
ct <- sub("^m", "", ct)
colnames(cell_log_umi) <- ct
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Midbrain.txt"), quote=F, row.names=F, sep="\t")

#### human data
umi <- fread(input="gzip -cd Linnarsson/GSE76381_EmbryoMoleculeCounts.cef.txt.gz", data.table=F)
samples <- data.frame(t(umi[2:4,3:ncol(umi)]), stringsAsFactors = F)
colnames(samples) <- umi[2:4,2]
genes <- umi[6:nrow(umi),1]
umi <- as.matrix(apply(umi[6:nrow(umi), 3:ncol(umi)], 2, as.numeric))
colnames(umi) <- samples$Cell_ID
rownames(umi) <- genes
samples <- samples[samples$Cell_type!="Unk",]
umi <- umi[,samples$Cell_ID]

### map genes to ENSG
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(umi), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
umi <- umi[genes$symbol[!is.na(genes$hs.ensg)],]
cell_zero_lt_80 <- cell_zero_lt_80[genes$symbol[!is.na(genes$hs.ensg)],]
notes <- c(notes, paste0("Number of genes with unique hs ENSG: ", nrow(umi)))
save(list=c("umi", "cell_zero_lt_80", "samples"), file=paste0("Linnarsson/", gse, "_human_umi.RData"))

### process
log_umi <- log2(umi+1)
ct <- sort(unique(samples$Cell_type))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$Cell_ID[samples$Cell_type==c]], 1, mean)
}
ct <- sub("^h", "", ct)
colnames(cell_log_umi) <- ct
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Human_Midbra.txt"), quote=F, row.names=F, sep="\t")

##### GSE78845 #####
### original study: Furlan et al. Nat. Neurosci. 2016 (PMIT: 27571008)
### expression data was downloaded from gSE78845
### metadata was extracted from family soft file
### see README for details
### Cells with label "Unclassified" were excluded

gse <- "GSE78845"
### samples
samples <- fread("Linnarsson/GSE78845_metadata.txt", data.table=F)
samples <- samples[samples$cell_type!="Unclassified",]

### expression
umi <- fread(input="gzip -cd Linnarsson/GSE78845_Furlan_et_al_Expression.tab.gz", data.table=F)
genes <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- genes
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
ct <- sort(unique(samples$cell_type))
cell_label <- list(Ach="Cholinergic_neurons_",
                   Glutamatergic="Glutamatergic_neurons",
                   Noradrenergic="Noradregergic_neurons")
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
for(n in names(cell_label)){
  ct <- sub(n, paste0(cell_label[n],"_"), ct)
}
colnames(cell_log_umi) <- ct
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Ganglia.txt"), quote=F, row.names=F, sep="\t")

##### GSE95315 #####
### original study: Hochgerner et al. Nat. Neurosci. 2018 (PMID: 29335606)
### expression data was downloaded from GSE95315
### metadata was extracted from family soft file
### see README for details
### cell types were edited to create 22 cell types as described in the original study
gse <- "GSE95315"

### samples
samples <- fread("Linnarsson/GSE95315_metadata.txt", data.table=F)
# merge cell types into 22 as presented in the paper
samples$cell_type <- sub(" glia-like", "", samples$cell_type)
samples$cell_type[!grepl("Neuroblast", samples$cell_type)] <- sub("\\d+$", "", samples$cell_type[!grepl("Neuroblast", samples$cell_type)])

### expression
umi <- fread(input="gzip -cd Linnarsson/GSE95315_10X_expression_data.tab.gz", data.table=F)
genes <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- genes

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
colnames(cell_log_umi) <- gsub(" ", "_", gsub("-", "_", ct))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-", "_", c))] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Dentate_gyrus.txt"), quote=F, row.names=F, sep="\t")

##### GSE95752 #####
### original study: original study: Hochgerner et al. Nat. Neurosci. 2018 (PMID: 29335606)
### expression data was downloaded from GSE95315
### metadata was extracted from family soft file
### see README for details
gse <- "GSE95752"

### samples
samples <- fread("Linnarsson/GSE95752_metadata.txt", data.table=F)

### expression
umi <- fread(input="gzip -cd Linnarsson/GSE95752_C1_expression_data.tab.gz", data.table=F)
genes <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- genes

### map symbol to hs ENSG
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
colnames(cell_log_umi) <- gsub(" ", "_", gsub("-", "_", ct))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-", "_", c))] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Dentate_gyrus.txt"), quote=F, row.names=F, sep="\t")

##### GSE104323 #####
### original study: Hochgerner et al. Nat. Neurosci. 2018 (PMID: 29335606)
### expression data was downloaded from GSE95315
gse <- "GSE104323"

### samples
samples <- fread(input="gzip -cd Linnarsson/GSE104323_metadata_barcodes_24185cells.txt.gz", data.table=F)
colnames(samples)[1] <- "cell_id"
colnames(samples)[7] <- "cell_type"
notes <- paste0("Number of samples: ", nrow(samples))
samples <- samples[samples$cell_type!="",]

### expression
umi <- fread(input="gzip -cd Linnarsson/GSE104323_10X_expression_data_V2.tab.gz", data.table=F)
genes <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- genes
umi <- umi[,samples$cell_id]

### map symbols to hs ENSG ID
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
colnames(cell_log_umi) <- gsub(" ", "_", gsub("-", "_", ct))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-", "_", c))] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Dentate_gyrus.txt"), quote=F, row.names=F, sep="\t")

##### GSE101601 #####
### original study: Hochgerner et al. Sci. Rep. 2017 (PMID: 29180631)
### expression data was obtained from GSE101601
### Human and mouse data sets were processed separately
gse <- "GSE101601"

#### Mouse data
umi <- fread(input="gzip -cd Linnarsson/GSE101601_molecules_mouse.csv.gz", data.table=F)
samples <- data.frame(cell_id=colnames(umi)[2:ncol(umi)], cell_type=as.character(umi[1,2:ncol(umi)]), stringsAsFactors = F)
genes <- umi[2:nrow(umi),1]
umi <- as.matrix(apply(umi[2:nrow(umi), 2:ncol(umi)], 2, as.numeric))
rownames(umi) <- genes

### map symbol to hs ENSG
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
cell_label <- list(astro="Astrocyte", cop="Differentiation_commited_oligodendrocyte_precursors",
                   endo="Endothelial", int="Interneuron", L="Pyramidal_cells_layer_",
                   mfol="Myelin_forming_oligodendrocyte", mgl="Microglia",
                   mol="Mature_oligodendrocyte", nfol="Myelin-forming_oligodendrotye",
                   opc="Olygodendrocyte_precursor_cell", peric="Pericytes",
                   pvm="Perivascular_macrophages", rbc="Red_blood_cell", vsm="Vascular_smooth_muscle")
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- gsub(" ", "_", gsub("-", "_", ct))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-", "_", c))] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
for(n in names(cell_label)){
  ct <- sub(n, cell_label[n], ct)
}
colnames(cell_log_umi) <- ct
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Somatosensory_cortex.txt"), quote=F, row.names=F, sep="\t")

#### Human data
umi <- fread(input="gzip -cd Linnarsson/GSE101601_molecules_human.csv.gz", data.table=F)
samples <- data.frame(cell_id=colnames(umi)[2:ncol(umi)], cell_type=as.character(umi[1,2:ncol(umi)]), stringsAsFactors = F)
genes <- umi[2:nrow(umi),1]
umi <- as.matrix(apply(umi[2:nrow(umi), 2:ncol(umi)], 2, as.numeric))
rownames(umi) <- genes

### map to ENSG
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(umi), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
umi <- umi[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_umi <- log2(umi+1)
ct <- sort(unique(samples$cell_type))
cell_label <- list(gaba_fbxl7="GABA_1", gaba_reln="GABA_2", gaba_vip="GABA_3",
                   gaba_sox6_tac1="GABA_4", gaba_sst="GABA_5", pdzd2="Glut_1",
                   col5a2="Glut_2", zfp804b="Glut_3", il1rapl2="Glut_4", dcc="Glut_5", htrc2="Glut_6")
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- gsub(" ", "_", gsub("-", "_", ct))
for(c in ct){
  cell_log_umi[,gsub(" ", "_", gsub("-", "_", c))] <- apply(log_umi[,samples$cell_id[samples$cell_type==c]], 1, mean)
}
colnames(cell_log_umi) <- cell_label[ct]
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Human_Temporal_cortex.txt"), quote=F, row.names=F, sep="\t")


##### GSE103840 #####
### original study: Haring et al. Nat. Neurosci. 2018 (PMID: 29686262)
### expression data was downloaded from GSE103840
gse <- "GSE103840"

#### expression data
umi <- fread(input="gzip -cd Linnarsson/GSE103840_molecule_counts.csv.gz", data.table=F)
samples <- data.frame(cell_id=colnames(umi)[2:ncol(umi)], t(umi[1:6, 2:ncol(umi)]), stringsAsFactors = F)
colnames(samples)[2:ncol(samples)] <- gsub(" ", "_", umi[1:6,1])
umi <- umi[7:nrow(umi),]
genes <- umi[,1]
umi <- as.matrix(apply(umi[, 2:ncol(umi)], 2, as.numeric))
rownames(umi) <- genes

### map symbol to hs ENSG
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
ct <- sort(unique(samples$cluster_name))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$cell_id[samples$cluster_name==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/Linnarsson_", gse, "_Mouse_Dorsal_horn.txt"), quote=F, row.names=F, sep="\t")
