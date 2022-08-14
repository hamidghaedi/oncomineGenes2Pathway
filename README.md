# oncomineGenes2Pathway
Mapping genes in the Oncomine Comprehensive Assay v3 (ocav3) to KEGG pathways

## Genes on the ocav3 panel:

### Hotspot genes:

AKT1, AKT2, AKT3, ALK, AR, ARAF, AXL, BRAF, BTK, CBL, CCND1, CDK4, CDK6, CHEK2, CSF1R, CTNNB1, DDR2, EGFR, ERBB2, ERBB3, ERBB4, ERCC2, ESR1, EZH2, FGFR1, FGFR2, FGFR3, FGFR4, FLT3, FOXL2, GATA2, GNA11, GNAQ, GNAS, H3F3A, HIST1H3B, HNF1A, HRAS, IDH1, IDH2, JAK1, JAK2, JAK3, KDR, KIT, KNSTRN, KRAS, MAGOH, MAP2K1, MAP2K2, MAP2K4, MAPK1, MAX, MDM4, MED12, MET, MTOR, MYC, MYCN, MYD88, NFE2L2, NRAS, NTRK1, NTRK2, NTRK3, PDGFRA, PDGFRB, PIK3CB, PIK3CA, PPP2R1A, PTPN11, RAC1, RAF1, RET, RHEB, RHOA, ROS1, SF3B1, SMAD4, SMO, SPOP, SRC, STAT3, TERT, TOP1, U2AF1, XPO1


### full-length gene 

ARID1A, ATM, ATR, ATRX, BAP1, BRCA1, BRCA2, CDK12, CDKN1B, CDKN2A, CDKN2B, CHEK1, CREBBP, FANCA, FANCD2, FANCI, FBXW7, MLH1, MRE11, MSH6, MSH2, NBN, NF1, NF2, NOTCH1, NOTCH2, NOTCH3, PALB2, PIK3R1, PMS2, POLE, PTCH1, PTEN, RAD50, RAD51, RAD51B, RAD51C, RAD51D, RNF43, RB1, SETD2, SLX4, SMARCA4, SMARCB1, STK11, TP53, TSC1, TSC2

I use the Broad institute  gene list called as MSigDB Collections to mapp genes to :

1. Hallmark pathways
2. oncogenic signature gene sets


```r
# loading packages
library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)


db <- org.Hs.egPATH
# Get the entrez gene identifiers that are mapped to a KEGG pathway ID
mapped_genes <- mappedkeys(db)
# Convert to a list
mapped_genesList <- as.list(db[mapped_genes])
# converting list to dataframe
df <- plyr::ldply (mapped_genesList, data.frame)
mappedDF = data.frame(ENTREZID = as.numeric(df[,1]), KEGG_ID = paste0("map",df[,2]))
# gen2pathway dataset
gen2path = aggregate(. ~ ENTREZID, mappedDF, FUN = function(x) 
  toString(x), na.action = NULL)
  
# merging hotspot and full-length genes:
#NB two genes "H3F3A" and "HIST1H3B" were changed to theri official gene symbols; "H3-3A" and "H3C2" , respectively.

geneList = c(
c("ARID1A", "ATM", "ATR", "ATRX", "BAP1", "BRCA1", "BRCA2", "CDK12", "CDKN1B", "CDKN2A", "CDKN2B", "CHEK1", "CREBBP", "FANCA", "FANCD2", "FANCI", "FBXW7", "MLH1", "MRE11", "MSH6", "MSH2", "NBN", "NF1", "NF2", "NOTCH1", "NOTCH2", "NOTCH3", "PALB2", "PIK3R1", "PMS2", "POLE", "PTCH1", "PTEN", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RNF43", "RB1", "SETD2", "SLX4", "SMARCA4", "SMARCB1", "STK11", "TP53", "TSC1", "TSC2")
, c("GNAS", "H3-3A", "H3C2", "HNF1A", "HRAS", "IDH1", "IDH2", "JAK1", "JAK2", "JAK3", "KDR", "KIT", "KNSTRN", "KRAS", "MAGOH", "MAP2K1", "MAP2K2", "MAP2K4", "MAPK1", "MAX", "MDM4", "MED12", "MET", "MTOR", "MYC", "MYCN", "MYD88", "NFE2L2", "NRAS", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PDGFRB", "PIK3CB", "PIK3CA", "PPP2R1A", "PTPN11", "RAC1", "RAF1", "RET", "RHEB", "RHOA", "ROS1", "SF3B1", "SMAD4", "SMO", "SPOP", "SRC", "STAT3", "TERT", "TOP1", "U2AF1", "XPO1")

)

# finding ENTREZ gene ID for all genes

gene.df <- bitr(geneList , fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db)
# changing column type
gene.df$ENTREZID = as.character(gene.df$ENTREZID)


pathways.hallmark <- fgsea::gmtPathways("~/mysigdb/h.all.v7.2.symbols.gmt")
# converting to a dataframe
pathways.hallmark <- plyr::ldply (pathways.hallmark, data.frame)
# naming the columns
colnames(pathways.hallmark) <- c("pathName", "ENTREZID")

oncoSig.pathways <- fgsea::gmtPathways("~/mysigdb/c6.all.v7.2.")
# converting to a dataframe
oncoSig.pathways <- plyr::ldply (oncoSig.pathways, data.frame)
# naming the columns
colnames(oncoSig.pathways) <- c("oncSigPathName", "ENTREZID")

# joing incomine genes with hallmark pathways
ocavHallmarkPath = dplyr::left_join(gene.df, pathways.hallmark)
ocavOncoSingPath = dplyr::left_join(gene.df, oncoSig.pathways)
```


