# load necessary libraries
library(biomaRt)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DRIMSeq)
library(EnsDb.Hsapiens.v86)
library(tibble)

# get gene symbols from ensembldb
edb <- EnsDb.Hsapiens.v86
txp <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"), return.type="data.frame")
txp <- txp |> as_tibble()
tx2gene <- txp |>
  dplyr::select(transcript=tx_id, gene=gene_id, symbol = gene_name)

# load count matrix
rbp <-read.table(gzfile("data/quantification_gencode.counts.txt.gz"),header = T)    # nolint
rbp$transcript <- unlist(lapply(strsplit(rbp$transcript,split='\\.'),'[[',1))

# define what samples are control and what are knockdowns (exp)
samples_exp <- names(rbp)[grepl('exp', names(rbp))]
samples_ctrl <- names(rbp)[grepl('ctr', names(rbp))]

# add experiment information (ctr, exp) to rbp
rbp <- rbp |> dplyr::select(transcript, samples_exp, samples_ctrl)
# filter transcripts that have no expression in any sample
rbp <- rbp[rowSums(rbp[c(samples_exp, samples_ctrl)])!=0,]

# add gene symbol based on transcript
rbp <- merge(rbp, tx2gene, by.x = "transcript", by.y = "transcript", all.x = TRUE) # nolint
# reorder columns
rbp <- rbp[c("gene", "symbol", "transcript", samples_exp, samples_ctrl)]

##########################################
## an alternative approach to the above ##
##########################################
#library(readr)
#rbp <- read_delim("data/quantification_gencode.counts.txt.gz")
#rbp <- rbp |>
#  mutate(transcript = sub("\\..*","",transcript)) |>
#  left_join(tx2gene) |>
#  dplyr::select(transcript, gene, contains("exp"), contains("ctr"))
#rbp <- rbp |>
#  mutate(rowsum = rowSums(dplyr::select(., contains(c("exp","ctr"))))) |>
#  dplyr::filter(rowsum > 0) |>
#  dplyr::select(-rowsum)
###########################################

## number of transcripts without gene id
rbp |> subset(is.na(gene)) |> nrow()
# filter these out
rbp_filt <- na.omit(rbp)
#rbp_filt = rbp_filt %>% subset(hgnc_symbol!="")
length(unique(rbp_filt$gene))

## Number of unique transcript
length(unique(rbp_filt$transcript))

NrTranscripts = rbp_filt |> group_by(gene)|> summarise(NrTranscript = n_distinct(transcript)) # nolint
ggplot(NrTranscripts) + geom_boxplot(aes(y=NrTranscript))
NrTranscripts |> arrange(-NrTranscript)

#rbp_filt %>% subset(hgnc_symbol=="MUC20-OT1") %>% select(transcript)

PTBP1_counts = rbp_filt |> subset(gene=="ENSG00000011304") 
PTBP1_counts_long = pivot_longer(PTBP1_counts, cols = c(samples_exp, samples_ctrl))
PTBP1_counts_long$condition = ifelse(grepl('exp', PTBP1_counts_long$name), 'knockout', 'control')
PTBP1_expression_plot <- ggplot(PTBP1_counts_long) + geom_boxplot(aes(x=condition, y = value, fill = condition)) + facet_wrap(~transcript)
PTBP1_expression_plot
# Save the plot to a file
ggsave(filename = "./plots/PTBP1_boxplot.pdf", plot = PTBP1_expression_plot, width = 10, height = 8, dpi = 600)

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
plot(explained_variance)
rotation <- as.data.frame(pca$rotation)
pcs <- as.data.frame(pca$x)
pcs$condition <- ifelse(grepl('exp', rownames(pcs)), 'knockout', 'control')
pcs$sample <- rownames(pcs) 

png('plots/PCA_highVar.png')
ggplot(pcs, aes(x=PC1, y = PC2, color = condition)) + geom_point() 
dev.off()

ggplot(pcs, aes(x=PC1, y = PC2, color = condition)) + geom_point()

##get features (transcripts) with highest loadings for PC1 and 2
pc1_idx <- rownames(rotation |> arrange(-PC1) |> head(10))
pc2_idx <- rownames(rotation |> arrange(-PC2) |> head(10))

##extract corresponding gene ids 
top10_pc1_trans <-  rbp_filt[pc1_idx,'transcript']
top10_pc1_gene <- rbp_filt[pc1_idx,c('symbol', 'gene')]

top10_pc2_trans <-  rbp_filt[pc2_idx,'transcript']
top10_pc2_gene <- rbp_filt[pc2_idx,c('symbol', 'gene')]

##plots all transcripts per gene 
plotTranscriptExp = function(gene, counts, plots){
  # Subset the counts data for the specified gene
  tmp_count <- counts |> dplyr::filter(gene == gene)
  # Pivot the data to a longer format
  tmp_count <- pivot_longer(tmp_count, cols = c(samples_exp, samples_ctrl))
  tmp_count$condition = ifelse(grepl('exp', tmp_count$name), 'knockout', 'control')
  
  plots[[gene]] =  ggplot(data.frame(tmp_count)) + 
                  geom_boxplot(aes(x=condition, y = value, fill = condition)) + 
                  facet_wrap(~transcript) + 
                  ggtitle(paste('Transcript isoforms of ', symbol))
  return(plots)
}
gene_input <- top10_pc1_gene[,2]
test <- rbp_filt |> dplyr::filter(gene == gene_input)
test <- pivot_longer(test, cols = c(samples_exp, samples_ctrl))
test$condition = ifelse(grepl('exp', test$name), 'knockout', 'control')
plots[[gene_input]] =  ggplot(data.frame(test)) + 
                  geom_boxplot(aes(x=condition, y = value, fill = condition)) + 
                  facet_wrap(~transcript) 

















# Initialize an empty list to store plots
plots <- list()

# Generate plots for top 10 genes
for (gene in unique(top10_pc1_gene$gene)) {
  plots <- plotTranscriptExp(gene, rbp_filt, plots)
}



p1 = lapply(top10_pc1_gene, plotTranscriptExp, counts = rbp_filt, plots = list()) # nolint
p1
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


?lapply
