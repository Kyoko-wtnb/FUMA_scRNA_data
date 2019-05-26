##### PsychENCODE data #####
### original study: Wang, D. et al. Science. 2018
### data was downloaded from http://resource.psychencode.org/
### Two datasets were created,
### one for developmental samples (DER-20_Single_cell_expresion_processed_TPM.tsv)
### and another for adult samples (DER-22_Single_cell_expression_raw_UMI.tsv)
### From adult dataset, 32 cels with cell type label "NA" were excluded.

##### Developmental PEC #####
### samples and expression
samples <- as.character(fread(cmd="gzip -cd PsychENCODE/DER-20_Single_cell_expression_processed_TPM.tsv.gz | head -1", header=F)[1,])
tpm <- fread(cmd="gzip -cd PsychENCODE/DER-20_Single_cell_expression_processed_TPM.tsv.gz", data.table=F)
g <- tpm[,1]
tpm <- as.matrix(tpm[,-1])
samples <- data.frame(id=1:ncol(tpm), celltype=sub("-","_",samples[-1]), stringsAsFactors = F)
colnames(tpm) <- 1:ncol(tpm)
rownames(tpm) <- g

### map to ENSG
ENSG <- fread("BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(tpm), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
genes$hs.ensg[is.na(genes$hs.ensg)] <- sapply(genes$symbol[is.na(genes$hs.ensg)], function(x){
  n <- grepl(paste0("|",x,"|"), ENSG$alias)
  if(length(n)==1){ENSG$ensembl_gene_id[n]}
  else{NA}
})
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
tpm <- tpm[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

### process
tpm_log <- log2(tpm+1)
ct <- sort(unique(samples$celltype))
cell_log_tpm <- matrix(nrow=nrow(tpm_log), ncol=length(ct))
colnames(cell_log_tpm) <- ct
for(c in ct){
  cell_log_tpm[,c] = apply(tpm_log[,samples$id[samples$celltype==c]], 1, mean)
}
cell_log_tpm <- data.frame(GENE=genes$hs.ensg, cell_log_tpm)
cell_log_tpm$Average <- apply(cell_log_tpm[,-1], 1, mean)
write.table(cell_log_tpm, "processed_data/PsychENCODE_Developmental.txt", quote=F, row.names=F, sep="\t")

##### Adult PEC #####
### samples and expression
samples <- as.character(fread(cmd="gzip -cd Data/PsychENCODE/DER-22_Single_cell_expression_raw_UMI.tsv.gz | head -1", header=F)[1,])
umi <- fread(cmd="gzip -cd PsychENCODE/DER-22_Single_cell_expression_raw_UMI.tsv.gz | sed 1d", header=F, data.table=F)
g <- umi[,1]
umi <- as.matrix(umi[,-1])
samples <- data.frame(id=1:ncol(umi), celltype=samples, stringsAsFactors = F)
colnames(umi) <- 1:ncol(umi)
rownames(umi) <- g

### map to ENSG
ENSG <- fread("BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(umi), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
genes$hs.ensg[is.na(genes$hs.ensg)] <- sapply(genes$symbol[is.na(genes$hs.ensg)], function(x){
  n <- grepl(paste0("|",x,"|"), ENSG$alias)
  if(length(n)==1){ENSG$ensembl_gene_id[n]}
  else{NA}
})
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
umi <- umi[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

### process
umi_log <- log2(umi+1)
ct <- sort(unique(samples$celltype[samples$celltype!="NA"]))
cell_log_umi <- matrix(nrow=nrow(umi_log), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] = apply(umi_log[,samples$id[samples$celltype==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg, cell_log_umi)
cell_log_umi$Average <- apply(cell_log_umi[,-1], 1, mean)
write.table(cell_log_umi, "processed_data/PsychENCODE_Adult.txt", quote=F, row.names=F, sep="\t")
