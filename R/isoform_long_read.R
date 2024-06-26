library(biomaRt)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DRIMSeq)

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
txp <- transcripts(edb, return.type="data.frame")
library(tibble)
txp <- txp |> as_tibble()
library(dplyr)
tx2gene <- txp |>
  dplyr::select(transcript=tx_id, gene=gene_id)

rbp <-read.table(gzfile("data/quantification_gencode.counts.txt.gz"),header = T)   
rbp$transcript <- unlist(lapply(strsplit(rbp$transcript,split='\\.'),'[[',1))
samples_exp = names(rbp)[grepl('exp', names(rbp))]
samples_ctrl = names(rbp)[grepl('ctr', names(rbp))]
rbp = rbp %>% dplyr::select(transcript, samples_exp, samples_ctrl)
rbp = rbp[rowSums(rbp[c(samples_exp, samples_ctrl)])!=0,]
rbp = merge(rbp, tx2gene, by.x = "transcript", by.y = "tx_id", all.x = TRUE)
rbp = rbp[c("gene_id","transcript", samples_exp, samples_ctrl)]

##########################################
## an alternative approach to the above ##
##########################################
library(readr)
counts <- read_delim("data/quantification_gencode.counts.txt.gz")
counts <- counts |>
  mutate(transcript = sub("\\..*","",transcript)) |>
  left_join(tx2gene) |>
  dplyr::select(transcript, gene, contains("exp"), contains("ctr"))
counts <- counts %>%
  mutate(rowsum = rowSums(dplyr::select(., contains(c("exp","ctr"))))) |>
  dplyr::filter(rowsum > 0) |>
  dplyr::select(-rowsum)
###########################################

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


#create colData df
sample_dat <- data.frame(sample_id = c(samples_ctrl, samples_exp))
sample_dat$condition = ifelse(grepl('exp', sample_dat$sample_id), 'knockout', 'control')


library(DESeq2)
##PCA of samples and look at gene loadings, DESeq object > plotPCA()
dds <- DESeqDataSetFromMatrix(countData = rbp_filt[-c(1,2)], colData = sample_dat, design =  ~ condition)
vsd <- vst(dds)
assay(vsd) # scaled, vst data

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

png('plots/explainedVar.png')
plot(explained_variance)
dev.off()

rotation <- as.data.frame(pca$rotation)
pcs <- as.data.frame(pca$x)
pcs$condition <- ifelse(grepl('exp', rownames(pcs)), 'knockout', 'control')
pcs$sample <- rownames(pcs) 

png('plots/PCA_highVar.png')
ggplot(pcs, aes(x=PC1, y = PC2, color = condition)) + geom_point() + 
  geom_text(aes(label = sample),show.legend=FALSE, vjust = -1, hjust = 0.5)
dev.off()


##get features (transcripts) with highest loadings for PC1 and 2
pc1_idx = rownames(rotation %>% arrange(-PC1) %>% head(10))
pc2_idx = rownames(rotation %>% arrange(-PC2) %>% head(10))

##extract corresponding gene ids 
top10_pc1_trans = rbp_filt[pc1_idx,'transcript']
top10_pc1_gene = rbp_filt[pc1_idx,'gene_id']

top10_pc2_trans = rbp_filt[pc2_idx,'transcript']
top10_pc2_gene = rbp_filt[pc2_idx,'gene_id']

##plots all transcripts per gene 
plotTranscriptExp = function(gene, counts, plots){
  tmp_count = counts %>% subset(gene_id==gene)
  tmp_count = pivot_longer(tmp_count, cols = c(samples_exp, samples_ctrl))
  tmp_count$condition = ifelse(grepl('exp', tmp_count$name), 'knockout', 'control')
  
  plots[[gene]] =  ggplot(data.frame(tmp_count)) + 
                  geom_boxplot(aes(x=condition, y = value, fill = condition)) + 
                  facet_wrap(~transcript) + 
                  ggtitle('Transcript isoforms of ', gene)
  return(plots)
}

p1 = lapply(top10_pc1_gene, plotTranscriptExp, counts = rbp_filt, plots = list())
p2 = lapply(top10_pc2_gene, plotTranscriptExp, counts = rbp_filt, plots = list())

pdf(paste0('plots/topTranscripts_PC1.pdf'))
p1
dev.off()


pdf(paste0('plots/topTranscripts_PC2.pdf'))
p2
dev.off()



##DrimSeq filtering
counts_drim = rbp_filt %>% dplyr::rename(feature_id=transcript)
d <- dmDSdata(counts=counts_drim, samples=sample_dat)

n <- 13
n.small <- 6
d_filt <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=13)

##we're filtering out the huge majority of transcripts and genes, super stringent!

length(unique(rbp_filt$gene_id)) 
length(unique(rbp_filt$transcript)) 

length(unique(counts(d_filt)$gene_id)) 
length(unique(counts(d_filt)$feature_id)) 

filterStats = data.frame(DrimFiltered = c('No', 'Yes'), 
           NrGenes = c(length(unique(rbp_filt$gene_id)), length(unique(counts(d_filt)$gene_id)) ),
           NrTranscripts = c(length(unique(rbp_filt$transcript)), length(unique(counts(d_filt)$feature_id))))
filterPlot = ggplot(filterStats %>% pivot_longer(cols = c('NrGenes', 'NrTranscripts'))) + 
              geom_bar(aes(x=name, y = value, fill = DrimFiltered), stat = 'identity', position = 'dodge')
filterPlot

png('plots/statsFiltering.png')
filterPlot
dev.off()



##DESeq2- To-Do?


