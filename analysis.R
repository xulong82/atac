library(dplyr)
library(genomation)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(ggplot2)
library(biomaRt)

rm(list = ls())
setwd("~/Dropbox/GitHub/ATAC/")

# === ATAC broadPeak
(files = list.files("broadPeak/"))
peaks = lapply(files, function(x) readBroadPeak(paste0("broadPeak/", x)))
names(peaks) = gsub("_S[1-9]_.*", "", files)

x = peaks[[1]]
summary(mcols(x)$score)  # is filter required?
sapply(peaks, length)

# === ATAC by gene structure
gene.parts = readTranscriptFeatures("bed/hg19_refseq_ucsc.bed")
annotateWithGeneParts(x, gene.parts)
names(gene.parts)
exons <- reduce(gene.parts$exons)
annotateWithFeature(target = x, feature = gene.parts$exons)
promoters <- reduce(gene.parts$promoters)
score_matrix <- ScoreMatrix(target = x, windows = gene.parts$promoters)
score_matrix <- ScoreMatrixBin(target = x, windows = promoters, bin.num = 50)
pdf("pdf/score_matrix.pdf", width = 10, height = 6)
heatMatrix(score_matrix)
dev.off()
plotMeta(score_matrix, xcoords = c(-1000, 1000))

annot.list = annotateWithGeneParts(GRangesList(peaks), gene.parts)
annot.list[[1]]@ perc.of.OlapFeat
dt <- sapply(annot.list, function(x) x@ perc.of.OlapFeat)
barplot(dt)
ggplot.dt <- data.frame(feature = rep(rownames(dt), 9), sample = rep(colnames(dt), each = 3), value = c(dt))
ggplot.dt$sample2 = gsub("-.*", "", ggplot.dt$sample)
pdf("pdf/gene_parts.pdf", width = 15, height = 10)
ggplot(ggplot.dt, aes(x = sample, y = value, fill = feature)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + xlab("") + ylab("Percentage") + 
  scale_fill_manual(values = c("grey70", "firebrick1", "dodgerblue3")) +
  theme(panel.border = element_blank(),
      axis.line = element_line(color = 'grey30'),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

atac_exons_gr <- lapply(peaks, function(x) subsetByOverlaps(gene.parts$exons, x))
atac_exons_genes <- lapply(atac_exons_gr, function(x) unique(x$name))

promoters = gene.parts$TSSes
start(promoters) = start(promoters) - 1e3
end(promoters) = end(promoters) + 1e3
atac_promoters_gr <- lapply(peaks, function(x) subsetByOverlaps(promoters, x))
atac_promoters_genes <- lapply(atac_promoters_gr, function(x) unique(x$name))

atac_promoters_genes <- lapply(atac_promoters_genes, function(x) {
  symbol = getBM("external_gene_name", "refseq_mrna", x, ensembl)
  symbol$external_gene_name %>% unique
})

refseq = atac_exons_genes[[1]]
refseq = atac_promoters_genes[[1]]

ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
symbol = getBM("external_gene_name", "refseq_mrna", refseq, ensembl)
symbol = symbol$external_gene_name %>% unique

# identify shared genes for astrocyte and neuron
# run GO and KEGG enrichment
# IPA pathway analysis

# === ATAC by ADSP variants
load("../Adsp/data/glmList.rdt"); list <- glmList
for(obj in names(list)) assign(obj, list[[obj]])
vep$Symbol <- gsub(".*SYMBOL=(.*)", "\\1", vep$Symbol)

gwas_lod <- filter(gwas, LOD > 15) # permutation cut
vep_lod <- filter(vep, UID %in% gwas_lod$UID)
gwas_vep <- cbind(vep_lod, gwas_lod[match(vep_lod$UID, gwas_lod$UID), -1])
gwas_vep <- dplyr::mutate(gwas_vep, start = POS, end = POS)
gwas_gr <- makeGRangesFromDataFrame(gwas_vep, keep.extra.columns = T)
gwas_gr <- renameSeqlevels(gwas_gr, paste0("chr", seqlevels(gwas_gr)))

intersect(gwas_gr, x)
y = subsetByOverlaps(gwas_gr, x)
subsetByOverlaps(gwas_gr, x) %>% reduce

atac_gwas_gr <- lapply(peaks, function(idx) subsetByOverlaps(gwas_gr, idx))

gene_astrocytes = lapply(atac_gwas_gr[grep("astrocytes", names(atac_gwas_gr))], function(x) x$Symbol) %>% unlist %>% unique
gene_neurons = lapply(atac_gwas_gr[grep("neurons", names(atac_gwas_gr))], function(x) x$Symbol) %>% unlist %>% unique

intersect(gene_astrocytes, gene_neurons)

# === Neuron/Astroctytes
neuron_file <- "Astrocytes_vs_neurons.HOMER_sorted_final_header_negative.txt"
astrocyte_file <- "Astrocytes_vs_neurons.HOMER_sorted_final_header_positive.txt"
neuron <- read.delim(paste0("diff/", neuron_file), stringsAsFactors = F)
astrocyte <- read.delim(paste0("diff/", astrocyte_file), stringsAsFactors = F)

neuron <- dplyr::select(neuron, -one_of("chr", "start", "end", "strand", "width"))
neuron_gr <- makeGRangesFromDataFrame(neuron, keep.extra.columns = T)
length(neuron_gr)

astrocyte <- dplyr::select(astrocyte, -one_of("chr", "start", "end", "strand", "width"))
astrocyte_gr <- makeGRangesFromDataFrame(astrocyte, keep.extra.columns = T)
length(astrocyte_gr)

neuron_gwas_gr = subsetByOverlaps(gwas_gr, neuron_gr)
neuron_gwas_gr %>% reduce
astrocyte_gwas_gr = subsetByOverlaps(gwas_gr, astrocyte_gr)
astrocyte_gwas_gr %>% reduce

atac_gwas = list(atac_gwas_gr = atac_gwas_gr, 
                 neuron_gwas_gr = neuron_gwas_gr, 
                 astrocyte_gwas_gr = astrocyte_gwas_gr)
