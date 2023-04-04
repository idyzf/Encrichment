library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(openxlsx)
library(stringr)
rm(list=ls())

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE, force=TRUE)
library(organism, character.only = TRUE)

# reading in data
df = read.xlsx("statistic_fc.xlsx", sheet=2)

# we want the log2 fold change 
original_gene_list <- df$FoldChange

# name the vector
names(original_gene_list) <- df$.y.

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.7, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
 #dev.off()


require(DOSE)
dotplot(gse, showCategory=10, split=".sign", font.size=8) + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)


x2 <- pairwise_termsim(gse) 
emapplot(x2)
emapplot_cluster(x2)


######################KEGG####################
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$.y. %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y
names(kegg_gene_list) <- str_replace( names(kegg_gene_list),
                                      pattern = ".[0-9]+$",
                                      replacement = "")

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 25,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

