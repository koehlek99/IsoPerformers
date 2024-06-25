library(biomaRt)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DRIMSeq)

rbp <-read.table(gzfile("data/quantification_gencode.counts.txt"),header = T)   
#rbp$transcript <- unlist(lapply(strsplit(rbp$transcript,split='\\.'),'[[',1))

samples_exp = names(rbp)[grepl('exp', names(rbp))]
samples_ctrl = names(rbp)[grepl('ctr', names(rbp))]
rbp = rbp %>% select(transcript, samples_exp, samples_ctrl)
rbp = rbp[rowSums(rbp[c(samples_exp, samples_ctrl)])!=0,]

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
anno <- getBM(attributes = c('ensembl_transcript_id_version','ensembl_gene_id', 'hgnc_symbol'),
              values = rbp$transcript,            
              filters = 'ensembl_transcript_id_version', 
              mart = mart)

rbp = merge(rbp, anno, by.x = "transcript", by.y = "ensembl_transcript_id_version", all.x = TRUE)
rbp = rbp[c("ensembl_gene_id", "hgnc_symbol","transcript", samples_exp, samples_ctrl)]

##number of transcripts without gene id
rbp %>% subset(is.na(ensembl_gene_id)) %>% nrow()
rbp %>% subset(is.na(hgnc_symbol)) %>% nrow()

rbp_filt = na.omit(rbp)
rbp_filt = rbp_filt %>% subset(hgnc_symbol!="")

##Number of unique genes
length(unique(rbp_filt$ensembl_gene_id))

##Number of unique transcript
length(unique(rbp_filt$transcript))

NrTranscripts = rbp_filt %>% group_by(hgnc_symbol) %>% summarise(NrTranscript = n_distinct(transcript))
ggplot(NrTranscripts) + geom_boxplot(aes(y=NrTranscript))
NrTranscripts %>% arrange(-NrTranscript)

rbp_filt %>% subset(hgnc_symbol=="MUC20-OT1") %>% select(transcript)

PTBP1_counts = rbp_filt %>% subset(hgnc_symbol=="PTBP1")
PTBP1_counts_long = pivot_longer(PTBP1_counts, cols = c(samples_exp, samples_ctrl))
PTBP1_counts_long$condition = ifelse(grepl('exp', PTBP1_counts_long$name), 'knockout', 'control')
png('plots/PTBPC_transcriptExp.png')
ggplot(PTBP1_counts_long) + geom_boxplot(aes(x=condition, y = value, fill = condition)) + facet_wrap(~transcript)
dev.off()

##PCA of samples and look at gene loadings, DESeq object > plotPCA()
dds <- DESeqDataSetFromMatrix(rbp, colData, design)
vsd <- vst(dds)
assay(vsd) # scaled, vst data

pca <- prcomp(t(assay(vsd)), center=TRUE)
pca$rotation




#filter out lowly expressed transcripts
sample_dat <- data.frame(sample_id = c(samples_ctrl, samples_exp))
sample_dat$condition = ifelse(grepl('exp', sample_dat$sample_id), 'knockout', 'control')

counts_drim = rbp_filt[-2] %>% rename(gene_id = ensembl_gene_id, feature_id = transcript)

d <- dmDSdata(counts=counts_drim, samples=sample_dat)
methods(class=class(d))
counts(d[1,])[,1:4]


n <- 13
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=13)

## An object of class dmDSdata 
## with 7764 genes and 12 samples
## * data accessors: counts(), samples()
