# rm(list=ls())

# This is very much a hacked version of the plots in RACER package (https://github.com/oliviasabik/RACER). In particular, the genetic annotation in the bottom panel is based on the "hg19" data object included in the RACER package.

library(coloc)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
# library(patchwork)
library(RACER)
library(reshape2)

genechr <- 19
hg19start <- 45409011
hg19end <- 45412650
genewindow <- 100

load("data/ldl_ad_APOE_coloc.RData")
# R objects required:
# 'coloc_results' -- colocalization results object
# 'd' -- a combined data frame of the GWAS summary statistics for trait 1 and trait 2
# 'D1' -- coloc dataset for trait 1
# 'D2' -- coloc dataset for trait 2
# 'LDmat' -- LD correlation matrix for the variants used

leadSNP_trait1 <- coloc_results[["results"]][
  which.max(coloc_results[["results"]]$lABF.df1),"snp"]

leadSNP_trait2 <- coloc_results[["results"]][
  which.max(coloc_results[["results"]]$lABF.df2),"snp"]

d_ld <- data.frame(rsid = colnames(LDmat), LDRsq = LDmat[leadSNP_trait1,]^2)
setDT(d_ld)

d_ld[,"LD" := cut(LDRsq, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
 labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))]

d <- d[d_ld, on = c("rsid" = "rsid"), nomatch = NULL]
d$lab <- ifelse(d$rsid %in% leadSNP_trait1, leadSNP_trait1, "")



d_ld2 <- data.frame(rsid = colnames(LDmat), LDRsq2 = LDmat[leadSNP_trait2,]^2)
setDT(d_ld2)

d_ld2[,"LD2" := cut(LDRsq2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                  labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))]

d <- d[d_ld2, on = c("rsid" = "rsid"), nomatch = NULL]
d$lab2 <- ifelse(d$rsid %in% leadSNP_trait2, leadSNP_trait2, "")



p1 <- ggplot(d, aes(x = POS, y = -log10(P_ldl), col = LD)) +
  geom_point() +
  scale_colour_manual(name = bquote(italic(R)^2~" with rs7412"),
   values = c(`0.8-1.0` = "red", `0.6-0.8` = "darkorange1",
              `0.4-0.6` = "green1", `0.2-0.4` = "skyblue1",
              `0.0-0.2` = "navyblue", `NA` = "grey"),
   labels = c(rev(levels(d$LD)), "Not available"),
   drop = FALSE) +
  labs(x = "", y = bquote(-log[10](italic(p))),
       title = "LDL-cholesterol") +
  geom_text_repel(aes(label = lab), min.segment.length = 0,
                  box.padding = 0.5, max.overlaps = Inf, col = "black") +
  geom_text_repel(aes(label = lab2), min.segment.length = 0,
                  box.padding = 0.5, max.overlaps = Inf, col = "black") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(45350000, 45400000, 45450000, 45500000),
                     labels = c(45350, 45400, 45450, 45500))

p2 <- ggplot(d, aes(x = POS, y = -log10(P_ad), col = LD)) +
  geom_point() +
  scale_colour_manual(name = bquote(italic(R)^2~" with rs7412"),
   values = c(`0.8-1.0` = "red", `0.6-0.8` = "darkorange1",
              `0.4-0.6` = "green1", `0.2-0.4` = "skyblue1",
              `0.0-0.2` = "navyblue", `NA` = "grey"),
   labels = c(rev(levels(d$LD)), "Not available"),
   drop = FALSE) +
  labs(x = "", y = bquote(-log[10](italic(p))),
         title = "AD risk") +
  geom_text_repel(aes(label = lab), min.segment.length = 0,
                  box.padding = 0.5, max.overlaps = Inf, col = "black") +
  geom_text_repel(aes(label = lab2), min.segment.length = 0,
                  box.padding = 0.5, max.overlaps = Inf, col = "black") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(45350000, 45400000, 45450000, 45500000),
                     labels = c(45350, 45400, 45450, 45500))

# p1 / p2


# add gene info
utils::data(hg19)
colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", 
                   "TRX_END", "LENGTH", "GENE_NAME", 
                   "TYPE")
gene_sub = hg19[hg19$CHR == genechr, ]
gene_sub = gene_sub[gene_sub$TYPE == "protein_coding",]

start <- hg19start - 1 - genewindow*1000
end <- hg19end + 1 + genewindow*1000

gene_sub = subset(gene_sub, gene_sub$TRX_START > (start - 
                                                    5000))
gene_sub = subset(gene_sub, gene_sub$TRX_END < (end + 5000))
gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID), ]
gene_sub = gene_sub[, c(3, 4, 6)]
gene_sub = reshape2::melt(gene_sub, id.vars = "GENE_NAME")
gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")

p_genes <- ggplot2::ggplot(gene_sub,
                    ggplot2::aes_string(x = "value", y = "y_value")) +
  ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) +
  ggplot2::theme_bw() +
  ggplot2::geom_text(data = plot_lab,
    ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
    hjust = -0.1, vjust = 0.3, size = 1.5) +
  ggplot2::theme(
    axis.title.y = ggplot2::element_text(color = "white", size = 28),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()) + 
  ggplot2::xlab(paste0("Chromosome ", genechr, " Position (kB)")) + 
  ggplot2::ylim(0, (max(gene_sub$y_value) + 1)) +
  scale_x_continuous(limits = c(start, end),
                     breaks = c(45350000, 45400000, 45450000, 45500000),
                     labels = c(45350, 45400, 45450, 45500))

p_out <- ggpubr::ggarrange(p1, p2, p_genes,
 heights = c(1.5, 1.5, 1), nrow = 3, 
 ncol = 1, common.legend = TRUE, legend = "right", align = "hv")

ggsave("results/coloc_plot_ldl_ad_APOE.png", p_out,
       height = 170, width = 170*(3/4), units = "mm", dpi = 300)
