##### DroNc data sets #####
### original study: Habib et al. Nat. Meth. 2017 (PMID: 28846088)
### Data is downloadeed from Broadinstitute Single Cell Portal at https://portals.broadinstitute.org/single_cell
### For human data set, expression data was fownloaded from GTEx at https://www.gtexportal.org/home/

library(data.table)

##### human data set #####
### cluster 1-14 and 16 are assigned in the original paper (my matching the tSNE plot), exclude cells with other clusters
### samples
samples <- fread(input="gzip -cd DroNc/GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.clusters.txt.gz", data.table=F)
colnames(samples) <- c("cell_id", "clst")
notes <- paste0("Number of samples in clst file: ", nrow(samples))
samples <- samples[samples$clst %in% c(1:14,16),]

### expression
umi <- fread(input="gzip -cd DroNc/GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.umi_counts.txt.gz", data.table=F)
g <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- g
umi <- umi[,samples$cell_id]

# manually assigne cell labels for each cluster
cell_label <- data.frame(clst=c(1:14,16), 
                         cell_type=c("exPFC1","exPFC2", "exCA1", "exCA3", "GABA1", "GABA2", "exDG", "ASC1", "ASC2", "ODC1", "ODC2", "OPC", "MG", "NSC", "END"),
                         stringsAsFactors = F)
samples$cell_type <- cell_label$cell_type[match(samples$clst, cell_label$clst)]

### map genes to ENSG
ENSG <- fread("~/Documents/VU/Data/BioMart/v92/ENSG.genes.txt", data.table = F)
genes <- data.frame(symbol=rownames(umi), stringsAsFactors = F)
genes$hs.ensg <- ENSG$ensembl_gene_id[match(genes$symbol, ENSG$external_gene_name)]
dup <- unique(genes$hs.ensg[duplicated(genes$hs.ensg)])
genes$hs.ensg[genes$hs.ensg %in% dup] <- NA
umi <- umi[genes$symbol[!is.na(genes$hs.ensg)],]
cell_zero_lt_80 <- cell_zero_lt_80[genes$symbol[!is.na(genes$hs.ensg)],]

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
write.table(cell_log_umi, paste0("processed_data/DroNc_Human_Hippocampus.txt"), quote=F, row.names=F, sep="\t")

##### mouse data set #####
### samples
samples <- fread(input="gzip -cd DroNc/Mouse_Meta_Data_with_cluster.txt.gz",data.table=F)
samples <- samples[2:nrow(samples),]
# remove un-informative cell types (the ones not assigned in the original study)
samples <- samples[!grepl("Unclassified", samples$Cluster) & samples$Cluster!="Doublets" & samples$Cluster!="ChP",]

### expression
umi <- fread(input="gzip -cd DroNc/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt.gz", data.table=F)
g <- umi[,1]
umi <- as.matrix(umi[,2:ncol(umi)])
rownames(umi) <- g
umi <- umi[,samples$NAME]

### map gsenes to hs ENSG
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
cell_zero_lt_80 <- cell_zero_lt_80[genes$symbol[!is.na(genes$hs.ensg)],]

### process
log_umi <- umi #already transformed
ct <- sort(unique(samples$Cluster))
cell_log_umi <- matrix(nrow=nrow(log_umi), ncol=length(ct))
colnames(cell_log_umi) <- ct
for(c in ct){
  cell_log_umi[,c] <- apply(log_umi[,samples$NAME[samples$Cluster==c]], 1, mean)
}
cell_log_umi <- data.frame(GENE=genes$hs.ensg[match(rownames(log_umi), genes$symbol)], cell_log_umi)
cell_log_umi$Average = apply(cell_log_umi[,2:ncol(cell_log_umi)], 1, mean)
write.table(cell_log_umi, paste0("processed_data/DroNc_Mouse_Hippocampus.txt"), quote=F, row.names=F, sep="\t")

