library(biomaRt)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DRIMSeq)

rbp <-read.table(gzfile("data/quantification_gencode.counts.txt.gz"),header = T)   
rbp$transcript <- unlist(lapply(strsplit(rbp$transcript,split='\\.'),'[[',1))

samples_exp = names(rbp)[grepl('exp', names(rbp))]
samples_ctrl = names(rbp)[grepl('ctr', names(rbp))]
rbp = rbp %>% dplyr::select(transcript, samples_exp, samples_ctrl)
rbp = rbp[rowSums(rbp[c(samples_exp, samples_ctrl)])!=0,]

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
txp <- transcripts(edb, return.type="data.frame")

library(tibble)
txp <- txp |> as_tibble()

library(dplyr)
tx2gene <- txp |>
  dplyr::select(tx_id, gene_id)


rbp = merge(rbp, tx2gene, by.x = "transcript", by.y = "tx_id", all.x = TRUE)
rbp = rbp[c("gene_id","transcript", samples_exp, samples_ctrl)]

##number of transcripts without gene id
rbp %>% subset(is.na(gene_id)) %>% nrow()

rbp_filt = na.omit(rbp)
#rbp_filt = rbp_filt %>% subset(hgnc_symbol!="")

##Number of unique genes
length(unique(rbp_filt$gene_id))

##Number of unique transcript
length(unique(rbp_filt$transcript))

NrTranscripts = rbp_filt %>% group_by(gene_id) %>% summarise(NrTranscript = n_distinct(transcript))
ggplot(NrTranscripts) + geom_boxplot(aes(y=NrTranscript))
NrTranscripts %>% arrange(-NrTranscript)

#rbp_filt %>% subset(hgnc_symbol=="MUC20-OT1") %>% select(transcript)

PTBP1_counts = rbp_filt %>% subset(gene_id=="ENSG00000011304")
PTBP1_counts_long = pivot_longer(PTBP1_counts, cols = c(samples_exp, samples_ctrl))
PTBP1_counts_long$condition = ifelse(grepl('exp', PTBP1_counts_long$name), 'knockout', 'control')
png('plots/PTBPC_transcriptExp.png')
ggplot(PTBP1_counts_long) + geom_boxplot(aes(x=condition, y = value, fill = condition)) + facet_wrap(~transcript)
dev.off()


#filter out lowly expressed transcripts
sample_dat <- data.frame(sample_id = c(samples_ctrl, samples_exp))
sample_dat$condition = ifelse(grepl('exp', sample_dat$sample_id), 'knockout', 'control')



library(DESeq2)
##PCA of samples and look at gene loadings, DESeq object > plotPCA()
dds <- DESeqDataSetFromMatrix(countData = rbp_filt[-c(1,2)], colData = sample_dat, design =  ~ condition)
vsd <- vst(dds)
head(assay(vsd)) # scaled, vst data

ntop <- 1000 # pick a number

# calculate the variance for each feature
rv <- rowVars(assay(vsd))

# select the ntop feature by variance
select <- head(order(rv, decreasing=TRUE), ntop)

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

#calculate explained variance
variance <- pca$sdev^2
explained_variance <- variance / sum(variance)
plot(explained_variance)

rotation <- as.data.frame(pca$rotation)
pcs <- as.data.frame(pca$x)
pcs$condition <- ifelse(grepl('exp', rownames(pcs)), 'knockout', 'control')
pcs$sample <- rownames(pcs) 

png('plots/PCA_highVar.png')
ggplot(pcs, aes(x=PC1, y = PC2, color = condition)) + geom_point() + 
  geom_text(aes(label = sample), vjust = -1, hjust = 0.5)
dev.off()

##get transcripts with highest eigenvectors for PC1
transcript_ids <- rbp_filt[select, ]$transcript
ordered_PC1 <- pca$x[order(pca$x[,"PC1"],decreasing = T),]
rownames(ordered_PC1)[1:10]
transcript_IDs_PC1 = list()
for (i in rownames(ordered_PC1)[1:10]){
  i = as.integer(i)
  transcript_IDs_i = rbp_filt$transcript[i]
  transcript_IDs_PC1 <- c(transcript_IDs_PC1, transcript_IDs_i)
}
unlist(transcript_IDs_PC1)

##get transcripts with highest eigenvectors for PC2
ordered_PC2 <- pca$x[order(pca$x[,"PC2"],decreasing = T),]
rownames(ordered_PC2)[1:10]
transcript_IDs_PC2 = list()
for (i in rownames(ordered_PC2)[1:10]){
  i = as.integer(i)
  transcript_IDs_i = rbp_filt$transcript[i]
  transcript_IDs_PC2 <- c(transcript_IDs_PC2, transcript_IDs_i)
}
unlist(transcript_IDs_PC2)






##To-Do
counts_drim = rbp_filt[-2] %>% rename(gene_id = ensembl_gene_id, feature_id = transcript)

d <- dmDSdata(counts=counts_drim, samples=sample_dat)
methods(class=class(d))
counts(d[1,])[,1:4]
rbp_filt$transcript[84066]
head(transcript_ids)
n <- 13
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=13)

## An object of class dmDSdata 
## with 7764 genes and 12 samples
## * data accessors: counts(), samples()
