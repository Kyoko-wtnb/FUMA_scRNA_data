#### prepare genes #####
#
# 6 Apr 2018
########################

library(data.table)
library(GenomicFeatures)

##### mm10 gtf file #####
txdb <- makeTxDbFromGFF("genome/UCSC_mm10.gft.gz", format="gtf")
exons.list.per.genes <- exonsBy(txdb, by="gene")
exonic.gene.sizes <- lapply(exons.list.per.genes, function(x){sum(width(reduce(x)))})
exonic.gene.sizes <- data.frame(gene=names(exonic.gene.sizes), Length=as.numeric(unlist(exonic.gene.sizes)))
rm(txdb, exons.list.per.genes)
Xref <- fread(input="gzip -cd genome/UCSC_mm10_Xref.txt.gz", data.table=F)
colnames(Xref)[1] <- "mm10.knownGene.name"
save(list=c("exonic.gene.sizes", "Xref"), file="processed/UCSC_mm10.RData")

##### mm10 mouse to human ENSG #####
library(biomaRt)
human_mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mouse_mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
ncbi_mouse <- fread(input="gzip -cd ~/Documents/VU/Data/NCBIgenes/Mus_musculus.gene_info.gz", data.table=F)
ncbi_mouse <- ncbi_mouse[,c(2,3,5,6)]
ncbi_mouse$ensg <- sub(".*Ensembl:(ENSMUSG\\d+).*", "\\1", ncbi_mouse$dbXrefs)
ncbi_mouse <- ncbi_mouse[,c(1,2,3,5)]
ncbi_mouse <- ncbi_mouse[grepl("^ENSM", ncbi_mouse$ensg),]

mm2hs <- getLDS(attributes=c("ensembl_gene_id"),
       filters="ensembl_gene_id", values=ncbi_mouse$ensg, mart=mouse_mart,
       attributesL=c("ensembl_gene_id"), martL=human_mart)
colnames(mm2hs) <- c("mm.ensg", "hs.ensg")

save(list=c("ncbi_mouse", "mm2hs"), file="genome/mm2hs.RData")

