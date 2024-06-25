rbp <-read.table(gzfile("/Users/aleynaeray/Downloads/quantification_gencode.counts.txt.gz"),header = T)   
rownames(rbp) <- rbp$transcript
transcript_ids <- rbp$transcript
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
head(listAttributes(ensembl, page="feature_page"))
anno <- getBM(attributes = c('ensembl_transcript_id','ensembl_transcript_id_version','ensembl_gene_id','ensembl_gene_id_version'),
              values = transcript_ids,            
              filters = 'ensembl_transcript_id_version', 
      mart = mart)
##
idx <- rbp$transcript %in% anno$ensembl_transcript_id_version
rbp_anno <- rbp[idx,]
cts <- rbp_anno[,81:93] 

counts <- data.frame(gene_id=anno$ensembl_gene_id_version,
                     feature_id=anno$ensembl_transcript_id_version,
                     cts)

#filter out lowly expressed transcripts
library(DRIMSeq)
samps <- data.frame(sample_id = colnames(counts[3:15]), condition=c(1,2,1,1,2,1,2,2,1,2,2,1,2))
d <- dmDSdata(counts=counts, samples=samps)
methods(class=class(d))
counts(d[1,])[,1:4]


n <- 12
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
d
## An object of class dmDSdata 
## with 7764 genes and 12 samples
## * data accessors: counts(), samples()
