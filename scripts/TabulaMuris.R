##### Tabula Muris data #####
### original study: The Tabula Muris Consortium et al. bioRxiv. 2018
### data was downloaded from 
### https://figshare.com/articles/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells/5715040 for FACS
### https://figshare.com/articles/Single-cell_RNA-seq_data_from_microfluidic_emulsion/5715025 for droplet
### FACS and droplet data sets were processed separately
### FACS data set contains 115 cell types from 20 tissues
### droplet data set contains 75 cells from 12 tissues
### We created data set with all tissue together and each tissue type separately for both FACS and droplet
### In total, 1 (all) + 20 (tissues) data sets for FACS
### and 1 (all) + 12 (tissues) data sets for droplet
### Cells with label "unknown" were excluded

library(data.table)
library(Matrix)

##### FACS data #####
### samples
samples <- fread("TabulaMuris/annotations_FACS.csv", data.table=F)
samples <- samples[!grepl("unknown", samples$cell_ontology_class),]

### expression
count <- fread(input="gzip -cd TabulaMuris/FACS_all.csv.gz", data.table=F)
g <- count[,1]
count <- as.matrix(count[,2:ncol(count)])
rownames(count) <- g
count <- count[,samples$cell] # not all samples in the count matrix exist in annotation file

### map symbol to hs ENSG
load("genome/mm2hs.RData")
genes <- data.frame(symbol = rownames(count), stringsAsFactors = F)
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

### convert to CPM
lib_size <- apply(count, 2, sum)/1000000
log_cpm <- log2(t(apply(count, 1, function(x){x/lib_size}))+1)

### process
ts_cell_log_cpm <- matrix(nrow=nrow(log_cpm), ncol=length(cell_label))
colnames(ts_cell_log_cpm) <- cell_label
ts <- unique(samples$tissue)
## all tissues together
for(t in ts){
  ct <- unique(samples$cell_ontology_class[samples$tissue==t])
  for(c in ct){
    ts_cell_log_cpm[,paste(gsub(" ","_",t),gsub(" ", "_",c), sep=":")] = apply(log_cpm[,samples$cell[samples$tissue==t & samples$cell_ontology_class==c]], 1, mean)
  }
}
ts_cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(count), genes$symbol)], ts_cell_log_cpm)
ts_cell_log_cpm$Average <- apply(ts_cell_log_cpm[,2:ncol(ts_cell_log_cpm)], 1, mean)
write.table(ts_cell_log_cpm, "processed_data/TabulaMuris_FACS_all.txt", quote=F, row.names=F, sep="\t")

## per tissue
for(t in ts){
  tmp_label <- unique(paste(gsub(" ","_",samples$tissue[samples$tissue==t]),gsub(" ", "_",samples$cell_ontology_class[samples$tissue==t]), sep=":"))
  tmp_cpm <- log_cpm[,samples$cell[samples$tissue==t]]
  ct <- unique(samples$cell_ontology_class[samples$tissue==t])
  cell_log_cpm <- matrix(nrow=nrow(tmp_cpm), ncol=length(ct))
  colnames(cell_log_cpm) <- gsub(" ", "_", ct)
  for(c in ct){
    cell_log_cpm[,gsub(" ", "_", c)] <- apply(tmp_cpm[,samples$cell[samples$tissue==t & samples$cell_ontology_class==c]], 1, mean)
  }
  cell_log_cpm <- data.frame(GENE=genes$hs.ensg[match(rownames(tmp_cpm), genes$symbol)], cell_log_cpm)
  cell_log_cpm$Average <- apply(cell_log_cpm[,2:ncol(cell_log_cpm)], 1, mean)
  write.table(cell_log_cpm, paste0("processed_data/TabulaMuris_FACS_", gsub(" ", "_", t), ".txt"), quote=F, row.names=F, sep="\t")
}

##### droplet #####
### samples 
samples <- fread("TabulaMuris/annotations_droplets.csv", data.table=F)
samples <- samples[!grepl("unknown", samples$cell_ontology_class),]

### expression
folders <- list.files("TabulaMuris/droplet")
umi <- c()
n_cells <- 0
for(f in folders){
  cat(f,"\n")
  pre <- sub(".+-(.+)", "\\1", f)
  barcode <- fread(paste0("TabulaMuris/droplet/", f, "/barcodes.tsv"), header=F)[[1]]
  barcode <- sub("(.+)-1", "\\1", barcode)
  genes <- fread(paste0("TabulaMuris/droplet/", f, "/genes.tsv"), header=F, data.table=F)
  tmp <- readMM(paste0("TabulaMuris/droplet/", f, "/matrix.mtx"))
  n_cells <- n_cells+ncol(tmp)
  cat("readMM ok\n")
  cols <- paste(pre, barcode, sep="_")
  tmp <- tmp[,which(cols %in% samples$cell)]
  cols <- cols[cols %in% samples$cell]
  tmp <- as.matrix(tmp)
  rownames(tmp) <- genes[,1]
  colnames(tmp) <- cols
  cat("names matrix ok\n")
  if(length(umi)==0){
    umi <- tmp
  }else{
    umi <- cbind(umi, tmp[rownames(umi),])
  }
  rm(tmp)
}
rm(genes, barcode)
umi <- umi[,samples$cell] # not all samples in umi matrix exist in annotation file

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
ts_cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(cell_label))
colnames(ts_cell_log_umi) <- cell_label
ts <- unique(samples$tissue)
## all tissues together
for(t in ts){
  ct <- unique(samples$cell_ontology_class[samples$tissue==t])
  for(c in ct){
    ts_cell_log_umi[,paste(gsub(" ","_",t),gsub(" ", "_",c), sep=":")] = apply(log_umi[,samples$cell[samples$tissue==t & samples$cell_ontology_class==c]], 1, mean)
  }
}
ts_cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(umi), genes$symbol)], ts_cell_log_umi)
ts_cell_log_umi$Average <- apply(ts_cell_log_umi[,2:ncol(ts_cell_log_umi)], 1, mean)
write.table(ts_cell_log_umi, "processed_data/TabulaMuris_droplet_all.txt", quote=F, row.names=F, sep="\t")

## per tissue
for(t in ts){
  tmp_label <- unique(paste(gsub(" ","_",samples$tissue[samples$tissue==t]),gsub(" ", "_",samples$cell_ontology_class[samples$tissue==t]), sep=":"))
  tmp_umi <- log_umi[,samples$cell[samples$tissue==t]]
  ct <- unique(samples$cell_ontology_class[samples$tissue==t])
  cell_log_umi <- matrix(nrow=nrow(tmp_umi), ncol=length(ct))
  colnames(cell_log_umi) <- gsub(" ", "_", ct)
  for(c in ct){
    cell_log_umi[,gsub(" ", "_", c)] <- apply(tmp_umi[,samples$cell[samples$tissue==t & samples$cell_ontology_class==c]], 1, mean)
  }
  cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(tmp_umi), genes$symbol)], cell_log_umi)
  cell_log_umi$Average <- apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
  write.table(cell_log_umi, paste0("processed_data/TabulaMuris_droplet_", gsub(" ", "_", t), ".txt"), quote=F, row.names=F, sep="\t")
}