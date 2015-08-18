library(dplyr)
library(tidyr)
library(genomation)
library(GenomicFeatures)
library(VariantAnnotation)
library(ggplot2)
library(biomaRt)
library(VennDiagram)
library(DiffBind)
library(pheatmap)

rm(list = ls())
setwd("~/Dropbox/GitHub/ATAC/")
ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
source("../../X/function.R")

# atac peaks and similarity
(files = list.files("broadPeak/")) # atac broadPeak (all samples) 
peaks = lapply(files, function(x) readBroadPeak(paste0("broadPeak/", x))) %>% GRangesList
(names(peaks) <- c(paste0("Ast", 1:3), paste0("Neu", 1:6)))

(peaks.union <- reduce(unlist(peaks), ignore.strand = T))
peaks.bi <- sapply(peaks, function(x) countOverlaps(peaks.union, x)) & 1 
(unique <- table(rowSums(peaks.bi))) # unique peak number: 64%
(unqiue1 <- colSums(peaks.bi & rowSums(peaks.bi) == 1)) # unique peak number each sample
sapply(peaks, length) # peak number each sample, filter required?
sapply(peaks, length) - unqiue1 
peaks.bi.select <- peaks.bi[rowSums(peaks.bi) > 1, ] + 0
pheatmap(cor(peaks.bi.select))

# atac peaks by gene structures (all samples)
gene.parts <- readTranscriptFeatures("bed/hg19_refseq_ucsc.bed") 
annot.list <- annotateWithGeneParts(peaks, gene.parts) 
scoreMatrix <- ScoreMatrix(target = peaks[[1]], windows = gene.parts$promoters)

(peaks.gene.parts <- sapply(annot.list, function(x) x@annotation) %>% as.data.frame)
(gene.parts.peaks <- sapply(annot.list, function(x) x@perc.of.OlapFeat) %>% as.data.frame)

temp = peaks.gene.parts 
temp = gene.parts.peaks
graph.dt <- gather(temp, sample, value)
graph.dt$feature <- rep(rownames(temp), 9)
graph.dt$sample <- gsub("-.*-", "", graph.dt$sample) 

pdf("pdf/gene_part1.pdf", width = 9, height = 5)
pdf("pdf/gene_part2.pdf", width = 9, height = 5)
ggplot(graph.dt, aes(x = sample, y = value, fill = feature)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + xlab("") + ylab("Percentage") + 
# scale_fill_manual(values = c("grey70", "grey10", "dodgerblue3", "firebrick1")) +
  scale_fill_manual(values = c("grey70", "dodgerblue3", "firebrick1")) +
  theme(legend.key = element_blank()) 
dev.off()

pdf("pdf/score_matrix.pdf", width = 10, height = 6)
my_palette = colorRampPalette(c("blue", "yellow", "red"))(n = 3e2)
heatMatrix(scoreMatrix, col = my_palette)
dev.off()

promoters <- gene.parts$TSSes
start(promoters) <- start(promoters) - 1e3
end(promoters) <- end(promoters) + 1e3
atac_promoters_gr <- lapply(peaks, function(x) subsetByOverlaps(promoters, x))
atac_promoters_genes <- lapply(atac_promoters_gr, function(x) {
  symbol = getBM("external_gene_name", "refseq_mrna", unique(x$name), ensembl)
  symbol$external_gene_name %>% unique
})

atac_promoters_genes_astrocyte <- unlist(atac_promoters_genes[1:3])
atac_promoters_genes_astrocyte <- names(which(table(atac_promoters_genes_astrocyte) == 3)) 
atac_promoters_genes_neuron <- unlist(atac_promoters_genes[4:9])
atac_promoters_genes_neuron <- names(which(table(atac_promoters_genes_neuron) == 6)) 

vennList <- list(Neuron = atac_promoters_genes_neuron, Astrocyte = atac_promoters_genes_astrocyte)
venn.diagram(vennList, imagetype = "png", file = "pdf/venn.png", width = 1500, height = 1500)

atac_promoters_genes_neuron_gk <- hsGK(atac_promoters_genes_neuron)
atac_promoters_genes_astrocyte_gk <- hsGK(atac_promoters_genes_astrocyte)
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

# atac peaks: neuron/Astroctytes-specific
neuron_file <- "Astrocytes_vs_neurons.HOMER_sorted_final_header_negative.txt"
astrocyte_file <- "Astrocytes_vs_neurons.HOMER_sorted_final_header_positive.txt"
neuron <- read.delim(paste0("diff/", neuron_file), stringsAsFactors = F)
astrocyte <- read.delim(paste0("diff/", astrocyte_file), stringsAsFactors = F)
neuron <- dplyr::select(neuron, -one_of("chr", "start", "end", "strand", "width"))
neuron <- makeGRangesFromDataFrame(neuron, keep.extra.columns = T)
astrocyte <- dplyr::select(astrocyte, -one_of("chr", "start", "end", "strand", "width"))
astrocyte <- makeGRangesFromDataFrame(astrocyte, keep.extra.columns = T)
peaks_diff <- GRangesList(neuron = neuron_gr, astrocyte = astrocyte_gr)

# with genes promoters
atac_diff_promoters_gr <- lapply(peaks_diff, function(x) subsetByOverlaps(promoters, x))
atac_diff_promoters_genes <- lapply(atac_diff_promoters_gr, function(x) {
  symbol = getBM("external_gene_name", "refseq_mrna", unique(x$name), ensembl)
  symbol$external_gene_name %>% unique
})

venn.diagram(atac_diff_promoters_genes, imagetype = "png", file = "pdf/venn2.png", width = 1500, height = 1500)
atac_diff_promoters_genes_gk <- lapply(atac_diff_promoters_genes, hsGK)
lapply(atac_diff_promoters_genes_gk, function(gk) {
  data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])
})

load("../Adsp/data/glmList.rdt") # GWAS
for(obj in names(glmList)) assign(obj, glmList[[obj]])
select <- filter(gwas, LOD > 15) # permutation cut
gwas_gr <- merge(vep, select, by.x = "UID", by.y = "UID")
gwas_gr <- makeGRangesFromDataFrame(gwas_gr, start.field = "POS", end.field = "POS", keep.extra.columns = T)
gwas_gr <- renameSeqlevels(gwas_gr, paste0("chr", seqlevels(gwas_gr)))

# atac by ADSP variants
(atac_gwas_gr <- lapply(peaks, function(idx) subsetByOverlaps(gwas_gr, idx)))
(gene_neurons <- lapply(atac_gwas_gr[grep("Neu", names(atac_gwas_gr))], function(x) unique(x$Symbol)))
(gene_neurons <- names(which(table(unlist(gene_neurons)) == 6)))
(gene_astrocytes <- lapply(atac_gwas_gr[grep("Ast", names(atac_gwas_gr))], function(x) unique(x$Symbol)))
(gene_astrocytes <- names(which(table(unlist(gene_astrocytes)) == 3)))
intersect(gene_astrocytes, gene_neurons)

# atac_diff by ADSP variants
atac_diff_gwas_gr <- lapply(peaks_diff, function(idx) subsetByOverlaps(gwas_gr, idx))
lapply(atac_diff_gwas_gr, function(x) mcols(x)[c("Feature", "Consequence", "Symbol")])
sapply(atac_gwas_gr, length)
