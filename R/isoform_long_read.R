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
dds <- DESeqDataSetFromMatrix(countData = rbp_filt[-c(1,2,3)], colData = sample_dat, design =  ~ condition)
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
#top10_pc1_gene <- rbp_filt[pc1_idx,c('symbol')]

top10_pc2_trans <-  rbp_filt[pc2_idx,'transcript']
top10_pc2_gene <- rbp_filt[pc2_idx,c('symbol', 'gene')]

##plots all transcripts per gene 
plotTranscriptExp = function(geneID, geneName, counts){
  # Subset the counts data for the specified gene
  tmp_count <- counts |> dplyr::filter(gene == geneID)
  # Pivot the data to a longer format
  tmp_count <- pivot_longer(tmp_count, cols = c(samples_exp, samples_ctrl))
  tmp_count$condition = ifelse(grepl('exp', tmp_count$name), 'knockout', 'control')
  
  p = ggplot(data.frame(tmp_count)) + 
                  geom_boxplot(aes(x=condition, y = value, fill = condition)) + 
                  facet_wrap(~transcript) + 
                  ggtitle(paste('Transcript isoforms of ', geneName))
  return(p)
}

# Initialize an empty list to store plots
plots_pc1 <- list()

for(gene in top10_pc1_gene$gene){
  message(gene)
  symbol = top10_pc1_gene[top10_pc1_gene$gene==gene, 'symbol']
  message(symbol)
  plotTranscriptExp(gene, symbol, counts = rbp_filt)
  plots_pc1[[gene]] = plotTranscriptExp(gene, symbol, counts = rbp_filt)
}

plots_pc2 <- list()

for(gene in top10_pc2_gene$gene){
  message(gene)
  symbol = top10_pc2_gene[top10_pc2_gene$gene==gene, 'symbol']
  message(symbol)
  plotTranscriptExp(gene, symbol, counts = rbp_filt)
  plots_pc2[[gene]] = plotTranscriptExp(gene, symbol, counts = rbp_filt)
}


pdf(paste0('plots/topTranscripts_PC1.pdf'))
plots_pc1
dev.off()


pdf(paste0('plots/topTranscripts_PC2.pdf'))
plots_pc2
dev.off()

##DrimSeq filtering
counts_drim = rbp_filt[-2] %>% dplyr::rename(gene_id = gene, feature_id=transcript)
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
           NrGenes = c(length(unique(rbp_filt$gene)), length(unique(counts(d_filt)$gene_id)) ),
           NrTranscripts = c(length(unique(rbp_filt$transcript)), length(unique(counts(d_filt)$feature_id))))
filterPlot = ggplot(filterStats %>% pivot_longer(cols = c('NrGenes', 'NrTranscripts'))) + 
              geom_bar(aes(x=name, y = value, fill = DrimFiltered), stat = 'identity', position = 'dodge')
filterPlot

png('plots/statsFiltering.png')
filterPlot
dev.off()

library(DEXSeq)
sample.data <- DRIMSeq::samples(d_filt)

count.data <- counts(d_filt)[,-c(1:2)]
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d_filt)$feature_id,
                     groupID=counts(d_filt)$gene_id)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet=TRUE)
dxd <- testForDEU(dxd)

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
dxr %>% data.frame() %>% mutate(sign = padj < 0.05) %>% dplyr::count(sign)


png('plots/featureLoadings_PC1.png')
plot((rotation |> arrange(-PC1))$PC1, ylab = "Feature loadings PC1") + abline(v=200, col = 'blue') + abline(v=800, col = 'red') #+ ylab('')
dev.off()


##get features (transcripts) with highest loadings for PC1 and 2
pc1_top200<- rownames(rotation |> arrange(-PC1) |> head(200))
pc1_top800<- rownames(rotation |> arrange(-PC1) |> head(800))

##extract corresponding gene ids 
top200_pc1_trans <-  rbp_filt[pc1_top200,'transcript']
top800_pc1_trans <-  rbp_filt[pc1_top800,'transcript']

library(VennDiagram)
dexseq_transcripts = (dxr %>% subset(padj < 0.05) %>% as.data.frame()%>% arrange(padj))$featureID

dexseq_transcripts_01 = (dxr %>% subset(padj < 0.01) %>% as.data.frame()%>% arrange(padj))$featureID
dexseq_genes_01 = (dxr %>% subset(padj < 0.01) %>% as.data.frame()%>% arrange(padj))$groupID

write.table(dexseq_transcripts_01, 'DEXSeq_sig0.01.txt', quote = F, row.names = F, col.names = F)
write.table(unique(dexseq_genes_01), 'DEXSeq_sig0.01_genes.txt', quote = F, row.names = F, col.names = F)

venn.diagram(x = list(top200_pc1_trans, dexseq_transcripts),
             category.names = c("PC1" , "DEXSeq"),
             filename = 'plots/Venn_PC_DEXseq.png',
             output=TRUE, # Output features,
             imagetype="png" ,
             height = 2000 , 
             width = 2000 , 
             resolution = 300,
             compression = "lzw")

venn.diagram(x = list(top800_pc1_trans, dexseq_transcripts),
             category.names = c("PC1- top transcripts" , "DEXSeq"),
             filename = 'plots/Venn_PCtop800_DEXseq.png',
             output=TRUE)

intersect(dexseq_transcripts, top200_pc1_trans)

dexseq_plt = list() 
top_genes = (dxr %>% data.frame() %>% arrange(padj) %>% head(10))$groupID
for(gene in top_genes){
  message(gene)
  symbol = unique(rbp_filt[rbp_filt$gene==gene, 'symbol'])
  message(symbol)
  dexseq_plt[[gene]] = plotTranscriptExp(gene, symbol, rbp_filt)
}

pdf('plots/DEXSeq_TranscriptExp.pdf')
dexseq_plt
dev.off()

(dxr %>% subset(padj < 0.01) %>% subset(groupID=='ENSG00000140416'))

tpm1_sign = rbp_filt |> 
  subset(transcript %in% (dxr %>% subset(padj < 0.01) %>% subset(groupID=='ENSG00000140416'))$featureID)
rownames(tpm1_sign) = tpm1_sign$transcript

png('plots/TPM1_heatmap.png')
tpm1_sign[-c(1:3)] %>% 
  as.matrix() %>% 
  heatmap(margins = c(15, 10), cexRow = 1, cexCol = 1)
dev.off()

tpm1_sign_long <- pivot_longer(tpm1_sign, cols = c(samples_exp, samples_ctrl))
tpm1_sign_long$condition = ifelse(grepl('exp', tpm1_sign_long$name), 'knockout', 'control')
p = ggplot(data.frame(tpm1_sign_long)) + 
  geom_boxplot(aes(x=condition, y = value, fill = condition)) + 
  facet_wrap(~transcript) + 
  ggtitle(paste('Transcript isoforms of TPM1'))

png('plots/TPM1_exp.png')
p
dev.off()

library(data.table)
gtf = fread('data/hg38.ensGene.gtf', sep = '\t') 
transcript_id <- str_extract(gtf$V9, 'transcript_id "\\S+"')
transcript_id <- str_replace(transcript_id, 'transcript_id ', '')
transcript_id <- str_replace_all(transcript_id, '"', '')
gtf$transcript_name = transcript_id

tpm_exons <- gtf %>% 
  subset(transcript_name %in% rownames(tpm1_sign)) %>% 
  dplyr::filter(V3 == "exon") %>% 
  rename(V1 = "seqnames", V3 = "type", V4 = "start", V5 = "end", V7 = "strand") %>% 
  dplyr::select(seqnames,  start, end, strand, type, transcript_name)

tpm1_sign$transcript
tpm1_sign$transcript 

tpm_exons$transcript_name = factor(tpm_exons$transcript_name, 
                                   levels = rev(c("ENST00000403994", "ENST00000560970", "ENST00000357980", "ENST00000560445", "ENST00000358278")))
library(ggtranscript)
transcriptStructures = tpm_exons %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes()
  ) +
  geom_intron(
    data = to_intron(tpm_exons, "transcript_name"),
    aes(strand = strand)
  )


png("plots/transStructures.png")
transcriptStructures
dev.off()
