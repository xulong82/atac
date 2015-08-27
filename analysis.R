library(ape)
library(amap)
library(dplyr)
library(tidyr)
library(reshape)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(VennDiagram)

rm(list = ls())
setwd("~/Dropbox/GitHub/ATAC/")
source("../../X/function.R")

# Summary of ATAC peaks
(files = list.files("broadPeak/")) # atac broadPeak (all samples) 
peaks = lapply(files, function(x) readBroadPeak(paste0("broadPeak/", x))) %>% GRangesList
(names(peaks) <- c(paste0("Ast", 1:3), paste0("Neu", 1:6)))
(peaks.All <- reduce(unlist(peaks), ignore.strand = T))

peaks.bi <- sapply(peaks, function(x) countOverlaps(peaks.All, x)) & 1 
table(rowSums(peaks.bi)) # unique peak number: 64%
(number.1 <- colSums(peaks.bi & rowSums(peaks.bi) == 1)) # unique peak number each sample
(number.All <- sapply(peaks, length)) # peak number each sample

gdt = data.frame(sample = names(peaks), unique = number.1, shared = number.All - number.1)
gdt = gather(gdt, type, value, one_of(c("unique", "shared")))
pdf("pdf/peaks.pdf", width = 6, height = 2.5)
ggplot(gdt, aes(x = sample, y = value, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + xlab("") + ylab("Number") + 
  scale_fill_manual(values = c("grey70", "firebrick1")) +
  theme(legend.key = element_blank()) 
dev.off()

peaks.All.select <- peaks.All[rowSums(peaks.bi) > 1, ]
hist(width(peaks.All.select), n = 1e3, col = "red", border = "red")

peaks.bi.select <- peaks.bi[rowSums(peaks.bi) > 1, ] + 0
hc1 <- hcluster(t(peaks.bi.select), method = "pearson", link = "average") %>% as.phylo
mycol = c(rep("blue", 3), rep("red", 6))
pdf("pdf/Hc.pdf", width = 5, height = 3)
plot(hc1, label.offset=1e-2, tip.color = mycol, direction="downward")
dev.off()

peaks.bi.select.Ast = rowSums(peaks.bi.select[, 1:3]) & 1
peaks.bi.select.Neu = rowSums(peaks.bi.select[, 4:9]) & 1
vennList = list(Ast = which(peaks.bi.select.Ast), Neu = which(peaks.bi.select.Neu))
venn.diagram(vennList, imagetype = "png", file = "pdf/venn1.png", width = 500, height = 500, resolution = 200)

# Genetic 
txdb = keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb.gr = GenomicRangesList(cds = cds(txdb), exons = exons(txdb), genes = genes(txdb), promoters = promoters(txdb))
txdb.gr$introns = unlist(intronsByTranscript(txdb))
txdb.gr$five = unlist(fiveUTRsByTranscript(txdb))
txdb.gr$three = unlist(threeUTRsByTranscript(txdb))
txdb.gr = lapply(txdb.gr, function(x) {strand(x) = "*"; reduce(x)})
txdb.gr$intergenic = gaps((txdb.gr$genes + 3e3))
txdb.gr$intergenic = txdb.gr$intergenic[strand(txdb.gr$intergenic) == "*", ]

overlap.bs = sapply(peaks, function(x) sapply(txdb.gr, function(y) sum(width(intersect(x, y)))))
overlap.ct.pk = sapply(peaks, function(x) sapply(txdb.gr, function(y) sum(countOverlaps(x, y) & 1))) # of peak
overlap.ct.gs = sapply(peaks, function(x) sapply(txdb.gr, function(y) sum(countOverlaps(y, x) & 1))) # of feature

(bs2pk = sweep(overlap.bs, 2, sapply(peaks, function(x) sum(width(x))), "/")) # to peak in bp
(bs2gs = sweep(overlap.bs, 1, sapply(txdb.gr, function(x) sum(width(x))), "/")) # to gene structure in bp
(ct2pk = sweep(overlap.ct.pk, 2, sapply(peaks, length), "/")) # to peak in count
(ct2gs = sweep(overlap.ct.gs, 1, sapply(txdb.gr, length), "/")) # to gene structure in count

col = rep("grey30", 8); col[c(3, 7)] = "firebrick1"
overlap_pctg = list(bs2pk = bs2pk, bs2gs = bs2gs, ct2pk = ct2pk, ct2gs = ct2gs)
pdf("pdf/overlap_pctg.pdf", width = 7, height = 4)
lapply(names(overlap_pctg), function(x) { 
  ggplot(melt(overlap_pctg[[x]]), aes(x = X1, y = value, fill = X1)) + geom_boxplot() + 
    scale_fill_manual(values = col) + xlab(x) + theme_bw() + theme(legend.position = "none")
}); dev.off()

openAll = sapply(peaks, function(x) sum(width(x))) / sum(as.numeric(seqlengths(txdb)))
odds_feature = sweep(bs2gs, 2, openAll, "/")
pdf("pdf/odds.pdf", width = 7, height = 4)
ggplot(melt(odds_feature), aes(x = X1, y = value, fill = X1)) + geom_boxplot() +
  scale_fill_manual(values = col) + theme_bw() + theme(legend.position = "none")
dev.off()

# Genes
genes = genes(txdb)
genes.tts = resize(genes, 1)
promoters = promoters(genes.tts, 2000, 200) # upstream:2000; downstream:200
promoters$symbol = select(org.Hs.eg.db, promoters$gene_id, columns=c("SYMBOL"), keytype="ENTREZID")$SYMBOL
open_promoters <- GRangesList(lapply(peaks, function(x) subsetByOverlaps(promoters, x)))
genesAst = names(which(table(unlist(open_promoters[1:3])$symbol) == 3))
genesNeu = names(which(table(unlist(open_promoters[4:9])$symbol) == 6))

vennList <- list(Neuron = genesNeu, Astrocyte = genesAst)
venn.diagram(vennList, imagetype = "png", file = "pdf/venn2.png", width = 500, height = 500, resolution = 200)

genesAst_gk <- hsGK(setdiff(genesAst, genesNeu)); gk = genesAst_gk
genesNeu_gk <- hsGK(setdiff(genesNeu, genesAst)); gk = genesNeu_gk
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

# Cell-specific open chromatin 
Neu_file <- "diff/Astrocytes_vs_neurons.HOMER_sorted_final_header_negative.txt"
Ast_file <- "diff/Astrocytes_vs_neurons.HOMER_sorted_final_header_positive.txt"
peaks_diff = lapply(c(Neu_file, Ast_file), function(file) {
  x1 = read.delim(file, stringsAsFactors = F)
  dplyr::select(x1, -one_of("chr", "start", "end", "strand", "width")) %>% makeGRangesFromDataFrame
}) %>% GRangesList
names(peaks_diff) = c("Neu", "Ast")

open_promoters_diff <- GRangesList(lapply(peaks_diff, function(x) subsetByOverlaps(promoters, x)))
open_promoters_diff_genes = lapply(open_promoters_diff, function(x) {y = x$symbol; y[! is.na(y)]})
venn.diagram(open_promoters_diff_genes, imagetype="png", file="pdf/venn3.png", width=500, height=500, resolution=200)
myGK <- lapply(open_promoters_diff_genes, hsGK)
lapply(myGK, function(x) {data.frame(KEGG=x$KEGG$Term[1:20], BP=x$GO$BP$Term[1:20], MF=x$GO$BP$Term[1:20])})

# GWAS
load("../Adsp/data/glmList.rdt")
select <- filter(glmList$gwas, LOD > 15) # permutation cut
select_gr <- makeGRangesFromDataFrame(select, start.field = "POS", end.field = "POS", keep.extra.columns = T)
select_gr <- renameSeqlevels(select_gr, paste0("chr", seqlevels(select_gr)))
open_select <- GRangesList(lapply(peaks, function(idx) subsetByOverlaps(select_gr, idx)))

(openVar_Ast = names(which(table(unlist(open_select[1:3])$UID) == 3)))
(openVar_Neu = names(which(table(unlist(open_select[4:9])$UID) == 6)))
openVar_Ast_gwas = merge(glmList$vep, select[select$UID %in% openVar_Ast, ], by.x = "UID", by.y = "UID")
openVar_Neu_gwas = merge(glmList$vep, select[select$UID %in% openVar_Neu, ], by.x = "UID", by.y = "UID")

openVar_diff = lapply(peaks_diff, function(idx) subsetByOverlaps(select_gr, idx))
openVar_diff_gwas = lapply(openVar_diff, function(x) merge(glmList$vep, select[select$UID %in% x$UID, ], by.x = "UID", by.y = "UID"))
lapply(openVar_diff_gwas, function(x) unique(x$Symbol))
