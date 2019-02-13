##### Linnarsson's lab Mouse Brain Atlas data #####
### original study: Zeizel et al. Cell. 2018
### Aggregated data at level 6 and level 5 were downloaded from http://mousebrain.org/
### l5_all.agg.loom for level 5
### l6_r1.agg.loom, l6_r2.agg.loom, l6_r3.agg.loom and l6_r4.agg.loom for level 6
### In total we created 5 data sets
### Note that we did not use individual cell level umi count
### The number of cell types (clusters) are
### 265 clusters for level 5
### 4 cell types for level 6 rank 1
### 6 cell types for level 6 rank 2
### 16 cell types for level 6 rank 3
### 39 cell types for level 6 rank 4
### Prior to this process on R, loom (hdf5) file was converted to plain text file using
### loom2text.py (available in the same github repository)

##### level 5 #####
pre <- "l5_all"
umi <- fread(input=paste0("gzip -cd Linnarsson/MouseBrainAtlas/",pre,"_matrix.txt.gz"), data.table=F, header=F)
samples <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_col_attrs.txt"), data.table=F)
genes <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_row_attrs.txt"), data.table=F)
colnames(umi) <- samples$ClusterName
rownames(umi) <- genes$Accession
log_umi <- log2(umi+1)
load("genome/mm2hs.RData")
genes$hs.ensg <- mm2hs$hs.ensg[match(genes$Accession, mm2hs$mm.ensg)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
genes <- genes[!is.na(genes$hs.ensg),]
log_umi <- log_umi[genes$Accession,]
log_umi$Average <- apply(log_umi, 1, mean)
log_umi <- data.frame(GENE=genes$hs.ensg, log_umi, stringsAsFactors = F, row.names=NULL)
write.table(log_umi, "processed_data/Linnarsson_MouseBrainAtlas_level5.txt", quote=F, row.names=F, sep="\t")

##### level 6 rank 1 #####
pre <- "l6_r1"
umi <- fread(input=paste0("gzip -cd Linnarsson/MouseBrainAtlas/",pre,"_matrix.txt.gz"), data.table=F, header=F)
samples <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_col_attrs.txt"), data.table=F)
genes <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_row_attrs.txt"), data.table=F)
colnames(umi) <- gsub(" ", "_", samples$TaxonomyRank1)
rownames(umi) <- genes$Accession
log_umi <- log2(umi+1)
genes$hs.ensg <- mm2hs$hs.ensg[match(genes$Accession, mm2hs$mm.ensg)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
genes <- genes[!is.na(genes$hs.ensg),]
log_umi <- log_umi[genes$Accession,]
log_umi <- log_umi[,sort(colnames(log_umi))]
log_umi$Average <- apply(log_umi, 1, mean)
log_umi <- data.frame(GENE=genes$hs.ensg, log_umi, stringsAsFactors = F, row.names=NULL)
write.table(log_umi, "processed_data/Linnarsson_MouseBrainAtlas_level6_rank1.txt", quote=F, row.names=F, sep="\t")

##### level 6 rank 2 #####
pre <- "l6_r2"
umi <- fread(input=paste0("gzip -cd Linnarsson/MouseBrainAtlas/",pre,"_matrix.txt.gz"), data.table=F, header=F)
samples <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_col_attrs.txt"), data.table=F)
genes <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_row_attrs.txt"), data.table=F)
colnames(umi) <- gsub(" ", "_", gsub("-", "_", samples$TaxonomyRank2))
rownames(umi) <- genes$Accession
log_umi <- log2(umi+1)
genes$hs.ensg <- mm2hs$hs.ensg[match(genes$Accession, mm2hs$mm.ensg)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
genes <- genes[!is.na(genes$hs.ensg),]
log_umi <- log_umi[genes$Accession,]
log_umi <- log_umi[,sort(colnames(log_umi))]
log_umi$Average <- apply(log_umi, 1, mean)
log_umi <- data.frame(GENE=genes$hs.ensg, log_umi, stringsAsFactors = F, row.names=NULL)
write.table(log_umi, "processed_data/Linnarsson_MouseBrainAtlas_level6_rank2.txt", quote=F, row.names=F, sep="\t")

##### level 6 rank 3 #####
pre <- "l6_r3"
umi <- fread(input=paste0("gzip -cd Linnarsson/MouseBrainAtlas/",pre,"_matrix.txt.gz"), data.table=F, header=F)
samples <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_col_attrs.txt"), data.table=F)
genes <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_row_attrs.txt"), data.table=F)
colnames(umi) <- gsub(" ", "_", gsub(",", "", gsub("-", "", samples$TaxonomyRank3)))
rownames(umi) <- genes$Accession
log_umi <- log2(umi+1)
genes$hs.ensg <- mm2hs$hs.ensg[match(genes$Accession, mm2hs$mm.ensg)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
genes <- genes[!is.na(genes$hs.ensg),]
log_umi <- log_umi[genes$Accession,]
log_umi <- log_umi[,sort(colnames(log_umi))]
log_umi$Average <- apply(log_umi, 1, mean)
log_umi <- data.frame(GENE=genes$hs.ensg, log_umi, stringsAsFactors = F, row.names=NULL)
write.table(log_umi, "processed_data/Linnarsson_MouseBrainAtlas_level6_rank3.txt", quote=F, row.names=F, sep="\t")

##### level 6 rank 4 #####
pre <- "l6_r4"
umi <- fread(input=paste0("gzip -cd Linnarsson/MouseBrainAtlas/",pre,"_matrix.txt.gz"), data.table=F, header=F)
samples <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_col_attrs.txt"), data.table=F)
genes <- fread(paste0("Linnarsson/MouseBrainAtlas/",pre,"_row_attrs.txt"), data.table=F)
colnames(umi) <- samples$TaxonomySymbol
rownames(umi) <- genes$Accession
log_umi <- log2(umi+1)
genes$hs.ensg <- mm2hs$hs.ensg[match(genes$Accession, mm2hs$mm.ensg)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
genes <- genes[!is.na(genes$hs.ensg),]
log_umi <- log_umi[genes$Accession,]
log_umi <- log_umi[,sort(colnames(log_umi))]
log_umi$Average <- apply(log_umi, 1, mean)
log_umi <- data.frame(GENE=genes$hs.ensg, log_umi, stringsAsFactors = F, row.names=NULL)
write.table(log_umi, "processed_data/Linnarsson_MouseBrainAtlas_level6_rank4.txt", quote=F, row.names=F, sep="\t")
