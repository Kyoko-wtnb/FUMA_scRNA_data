# FUMA scRNA-seq data sets

This is a github repository for pre-process script of scRNA-seq data sets used on FUMA web application ([http://fuma.ctglab.nl](http://fuma.ctglab.nl)).
Processed data sets can also be downloaded from this repository to run MAGMA by yourself.

## updates
* 19th May 2019:  
Added: PsychENCODE, GSE97478, GSE106707
* 13th Feb 2019:  
Added: Allene Brain Atlas Cell Type second release,10X dataset   
Updated: Tabula Muris FACS
* 19th Jul 2018:  
First release

## Preparation of genes
### Human Ensembl genes (GRCh37 v92)
The script is available under scripts/human_ENSG.R.
The output file is ENSG.genes.txt

### Mapping of mouse gene ID to human gene ID
The script is available under scripts/mm2hs.R.
The output file is mm2hs.RData

## Data sets
The list of data sets are available at [FUMA tutorial page](http://fuma.ctglab.nl/tutorial#celltype).

Preprocess of each data set is available in one of the following R scripts.
1. Allen_Brain_Atlas_Cell_Type.R  
2 human and 3 mousebrain data sets from Allen Brain Atlas.
2. Allen_Brain_Atlas_Cell_Type_prev.R  
3 mouse data sets from Allen Brain Atlas (previously released datasets).
3. DroNc.R  
2 data sets (human and mouse) from Habib et al. Nat. Meth. (2017)
4. DropViz.R  
Data set from [http://dropviz.org](http://dropviz.org).
5. Linnarsson_MouseBrainAtlas.R  
Data sets from [http://mousebrain.org](http://mousebrain.org).
6. TabulaMuris.R  
2 data sets (FACS and droplet) from The Tabula Muris Consortium.
6. MouseCellAtlas.R  
Data set from [http://bis.zju.edu.cn/MCA/](http://bis.zju.edu.cn/MCA/).
8. Linnersson_lab.R  
Other data sets from Linnarsson's group.
9. 10X.R  
PBMC data set downloaded from 10X Genomics.
10. PsychENCODE.R  
2 data sets from PsychENCODE.
10. GEO.R  
Everything else.

## Citation
A publication related to this study is under preparation.

## Additional process prior to the process in R scripts
### 1. GSE67602 (Linnarsson's lab)
Metadata for cells were extracted from family soft file by the following commands.
```
echo -e "cell_id\tcell_type" > GSE67602_celltype.txt
gzip -cd GSE67602_family.soft.gz | grep Sample_title | sed 's/!Sample_title = //' > title.txt
gzip -cd GSE67602_family.soft.gz | grep "cell type level 1" | sed 's/!Sample_characteristics_ch1 = cell type level 1: //' > celltype.txt
paste title.txt celltype.txt >>GSE67602_celltype.txt
```
### 2. GSE75330 (Linnarsson's lab)
Metadata for cells were extracted from family soft file by the following commands.
```
echo -e "cell_id\ttissue\tage\ttreatment\tcell_type" > GSE75330_metadata.txt
gzip -cd GSE75330_family.soft.gz | grep Sample_title | sed 's/^!Sample_title = //' >title.txt
gzip -cd GSE75330_family.soft.gz | grep Sample_source_name_ch1 | sed 's/^!Sample_source_name_ch1 = //' >tissue.txt
gzip -cd GSE75330_family.soft.gz | grep "= age:" | sed 's/^!Sample_characteristics_ch1 = age: //' >age.txt
gzip -cd GSE75330_family.soft.gz | grep "= treatment:" | sed 's/^!Sample_characteristics_ch1 = treatment: //' >treatment.txt
gzip -cd GSE75330_family.soft.gz | grep "inferred cell type:" | sed 's/^!Sample_characteristics_ch1 = inferred cell type: //' >celltype.txt
paste title.txt tissue.txt age.txt treatment.txt celltype.txt >>GSE75330_metadata.txt
```
### 3. GSE78845 (Linnarsson's lab)
Metadata for cells were extracted from family soft file by the following commands.
```
echo -e "cell_id\tcell_type\tage" >GSE78845_metadata.txt
gzip -cd GSE78845_family.soft.gz | grep Sample_title | sed 's/^!Sample_title = //' >title.txt
gzip -cd GSE78845_family.soft.gz | grep celltype | sed 's/^!Sample_characteristics_ch1 = celltype identifier: //' >celltype.txt
gzip -cd GSE78845_family.soft.gz | grep "= age:" | sed 's/^!Sample_characteristics_ch1 = age: //' >age.txt
paste title.txt celltype.txt age.txt >>GSE78845_metadata.txt
```
### 4. GSE95315 (Linnarsson's lab)
Metadata for cells were extracted from family soft file by the following commands.
```
echo -e "cell_id\tcell_type\tpostnatal_day" >GSE95315_metadata.txt
gzip -cd GSE95315_family.soft.gz | grep Sample_title | sed 's/^!Sample_title = //' >title.txt
gzip -cd GSE95315_family.soft.gz | grep "cell cluster:" | sed 's/^!Sample_characteristics_ch1 = cell cluster: //' >celltype.txt
gzip -cd GSE95315_family.soft.gz | grep "postnatal day:" | sed 's/^!Sample_characteristics_ch1 = postnatal day: //' >age.txt
paste title.txt celltype.txt age.txt >>GSE95315_metadata.txt
```
### 5. GSE95752 (Linnarsson's lab)
Metadata for cells were extracted from family soft file by the following commands.
```
echo -e "cell_id\tcell_type\tpostnatal_day" >GSE95752_metadata.txt
gzip -cd GSE95752_family.soft.gz | grep Sample_title | sed 's/^!Sample_title = //' >title.txt
gzip -cd GSE95752_family.soft.gz | grep "cell cluster:" | sed 's/^!Sample_characteristics_ch1 = cell cluster: //' >celltype.txt
gzip -cd GSE95752_family.soft.gz | grep "postnatal day:" | sed 's/^!Sample_characteristics_ch1 = postnatal day: //' >age.txt
paste title.txt celltype.txt age.txt >>GSE95752_metadata.txt
```
### 6. GSE81547
Metadata for cells were extracted from family soft file and expression data was concatanaged into a single matrix by the following commands.
```
# metadata
gzip -cd GSE81547_family.soft.gz | grep "SAMPLE = " | sed 's/^\^SAMPLE = //' >sample.txt
gzip -cd GSE81547_family.soft.gz | grep "donor_age:" | sed 's/^!Sample_characteristics_ch1 = donor_age: //' >age.txt
gzip -cd GSE81547_family.soft.gz | grep "inferred_cell_type:" | sed 's/^!Sample_characteristics_ch1 = inferred_cell_type: //' >celltype.txt
echo -e "cell_id\tage\tcell_type" >GSE81547_metadata.txt
paste sample.txt age.txt celltype.txt >> GSE81547_metadata.txt
rm sample.txt age.txt celltype.tt

# expression matrix
mkdir GSE81547_raw
mv GSE81547_RAW.tar GSE81547_raw
cd GSE81547_raw
tar -xvf GSE81547_RAW.tar
files=($(ls GSM*))
f=${files[0]}
pre=${f%%_*}
echo -e "gene\t$pre" >../GSE81547_count.txt
gzip -cd $f >>../GSE81547_count.txt
echo ${#files[@]}
for i in {1..2543}; do
        f=${files[$i]}
        pre=${f%%_*}
        echo $pre >tmp.txt
        gzip -cd $f | cut -f 2 >>tmp.txt
        paste ../GSE81547_count.txt tmp.txt >exp.txt
        mv exp.txt ../GSE81547_count.txt
done
cd ../
rm GSE81547_raw/*.gz
gzip GSE81547_count.txt
```
### 7. GSE67835
Metadata for cells were extracted from family soft file and expression data was concatanaged into a single matrix by the following commands.
```
# metadata
echo -e "cell_id\tcell_type\tage" > GSE67835_metadata.txt
gzip -cd GSE67835_family.soft.gz | grep "SAMPLE = " | sed 's/^\^SAMPLE = //' >sample.txt
gzip -cd GSE67835_family.soft.gz | grep "= cell type:" | sed 's/!Sample_characteristics_ch1 = cell type: //' > celltype.txt
gzip -cd GSE67835_family.soft.gz | grep "= age: " | sed 's/!Sample_characteristics_ch1 = age: //'>age.txt
paste sample.txt celltype.txt age.txt >>GSE67835_metadata.txt

# expression matrix
mkdir GSE67835_raw
mv GSE67835_RAW.tar GSE67835_raw
cd GSE67835_raw
tar -xvf GSE67835_RAW.tar
files=($(ls GSM*))
f=${files[0]}
pre=${f%%_*}
echo -e "gene\t$pre" >../GSE67835_count.txt
gzip -cd $f >>../GSE67835_count.txt
echo ${#files[@]}
for i in {1..465}; do
	f=${files[$i]}
	pre=${f%%_*}
	echo $pre >tmp.txt
	gzip -cd $f | cut -f 2 >>tmp.txt
	paste ../GSE67835_count.txt tmp.txt >exp.txt
	mv exp.txt ../GSE67835_count.txt
done
cd ../
rm GSE67835_raw/*.gz
gzip GSE67835_count.txt
```
### 8. GSE89232
Metadata for cells were extracted from family soft file by the following commands.
```
gzip -cd GSE89232_family.soft.gz | grep Sample_title | sed 's/^!Sample_title = //' >title.txt
gzip -cd GSE89232_family.soft.gz | grep "cell type:" | grep -v "(" | sed 's/!Sample_characteristics_ch1 = cell type: //' >celltype.txt
echo -e "cell_id\tcell_type" >GSE89232_celltype.txt
paste title.txt celltype.txt >>GSE89232_celltype.txt
```

### 9. GSE97478
Metadata for cells were extracted from family soft file by the following commands.
```
gzip -cd GSE97478_family.soft.gz | grep Sample_title | sed 's/^!Sample_title = //' >title.txt
gzip -cd GSE97478_family.soft.gz | grep "cell type:" | grep -v "(" | sed 's/!Sample_characteristics_ch1 = cell type: //' >celltype.txt
echo -e "cell_id\tcell_type" >GSE97478_celltype.txt
paste title.txt celltype.txt >>GSE97478_celltype.txt
```

### 9. GSE106707
Metadata for cells were extracted from family soft file by the following commands.
```
gzip -cd GSE106707_family.soft.gz | grep Sample_title | sed 's/^!Sample_title = //' >title.txt
gzip -cd GSE106707_family.soft.gz | grep "cell type:" | grep -v "(" | sed 's/!Sample_characteristics_ch1 = cell type: //' >celltype.txt
gzip -cd GSE106707_family.soft.gz | grep "postnatal days:" | grep -v "(" | sed 's/!Sample_characteristics_ch1 = postnatal days: //' >pd.txt
echo -e "cell_id\tcell_type\tpd" >GSE106707_celltype.txt
paste title.txt celltype.txt pd.txt >>GSE106707_celltype.txt
```
