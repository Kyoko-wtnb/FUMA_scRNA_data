##### 10X PBMC data set #####
### original study: Zheng et al. NAt. Commun. 2017
### Data was downloaded from 10X Genomics website
### https://community.10xgenomics.com/t5/Data-Sharing/10x-Single-Cell-3-Paper-Zheng-et-al-2016-Datasets/td-p/231

library(data.table)
library(Matrix)

##### expression data #####
barcode <- fread("10x/filtered_matrices_mex/hg19/barcodes.tsv", header=F)[[1]]
genes <- fread("10x/filtered_matrices_mex/hg19/genes.tsv", data.table=F, header=F)
umi <- readMM("10x/filtered_matrices_mex/hg19/matrix.mtx")
umi <- as.matrix(umi)
colnames(umi) <- barcode
rownames(umi) <- genes$V1

##### samples #####
samples <- fread("10x/68k_pbmc_barcodes_annotation.tsv", data.table=F)
samples$celltype <- gsub(" ", "_", samples$celltype)

##### map genes to ENSG #####
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(umi), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
umi <- umi[genes$symbol[!is.na(genes$hs.ensg)],]
genes <- genes[!is.na(genes$hs.ensg),]

##### process #####
log_umi <- log2(umi+1)
ct <- sort(unique(samples$celltype))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$barcodes[samples$celltype==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average <- apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, "processed_data/PBMC_10x_68k.txt", quote=F, row.names=F, sep="\t")
