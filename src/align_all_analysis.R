library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
theme_set(theme_classic())
setwd(dir = "/home/julien/Documents/Cours/M2/projet_long/compseq/results")

circRNA_GDP5L_1 <- c(625,1323) 
circRNA_GDP5L_2 <- c(28003,28270)
circRNA_GDP5L_3 <- c(29424,29890)
circRNA_GDP7L_1 <- c(13629,13930)
circRNA_GDP7L_2 <- c(18694,19022)
circRNA_GDP8L_1 <- c(12120,12729)
circRNA_GDP8L_2 <- c(16230,16741)
circRNA_GDP8L_3 <- c(18694,19374)
circRNA_GDP8L_4 <- c(20677,21434)
circRNA_GDP8L_5 <- c(28581,29135)

#--------------------- Distribution conserved motifs all k all sequences ------------#

data <- read.table(file = 'all_k.tsv', sep = '\t', header = TRUE)
filtered_data <- data[!duplicated(data[,c('start','comp_start')]),]
totally_conserved <- subset(filtered_data, count == 53)
totally_conserved$k <- as.factor(totally_conserved$k)

g <- ggplot(totally_conserved, aes(start)) + scale_fill_brewer(palette = "Spectral")
g + geom_histogram(aes(fill=k),
                   binwidth = 1000, 
                   col="black", 
                   size=0.1) +
  labs(title="Distribution on the genome of conserved motifs across all the sequence")

write.table(totally_conserved, file = "totally_conserved.tsv", sep = "\t",
            row.names = F, col.names = T)

#------------------- Mapping on circRNA all k all sequences-----------------------#


data <- read.table(file = 'all_k.tsv', sep = '\t', header = TRUE)
filtered_data <- data[!duplicated(data[,c('start','comp_start')]),]
totally_conserved <- subset(filtered_data, count == 53)
totally_conserved$k <- as.factor(totally_conserved$k)

map_circRNA_GDP5L_1 <- subset(totally_conserved, start > circRNA_GDP5L_1[1] & comp_start < circRNA_GDP5L_1[2])
map_circRNA_GDP5L_2 <- subset(totally_conserved, start > circRNA_GDP5L_2[1] & comp_start < circRNA_GDP5L_2[2])
map_circRNA_GDP5L_3 <- subset(totally_conserved, start > circRNA_GDP5L_3[1] & comp_start < circRNA_GDP5L_3[2])
map_circRNA_GDP7L_1 <- subset(totally_conserved, start > circRNA_GDP7L_1[1] & comp_start < circRNA_GDP7L_1[2])
map_circRNA_GDP7L_2 <- subset(totally_conserved, start > circRNA_GDP7L_2[1] & comp_start < circRNA_GDP7L_2[2])
map_circRNA_GDP8L_1 <- subset(totally_conserved, start > circRNA_GDP8L_1[1] & comp_start < circRNA_GDP8L_1[2])
map_circRNA_GDP8L_2 <- subset(totally_conserved, start > circRNA_GDP8L_2[1] & comp_start < circRNA_GDP8L_2[2])
map_circRNA_GDP8L_3 <- subset(totally_conserved, start > circRNA_GDP8L_3[1] & comp_start < circRNA_GDP8L_3[2])
map_circRNA_GDP8L_4 <- subset(totally_conserved, start > circRNA_GDP8L_4[1] & comp_start < circRNA_GDP8L_4[2])
map_circRNA_GDP8L_5 <- subset(totally_conserved, start > circRNA_GDP8L_5[1] & comp_start < circRNA_GDP8L_5[2])

g <- ggplot(map_circRNA_GDP5L_3, aes(start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP5L_3)),
                   binwidth = 10, 
                   col="black", 
                   size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP5L_3[1], circRNA_GDP5L_3[2]) +
  labs(title="Mapping of motifs on circRNA_GDP5L_3")
gg <- ggplot(map_circRNA_GDP5L_3, aes(comp_start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP5L_3)),
                   binwidth = 10, 
                   col="black", 
                   size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP5L_3[1], circRNA_GDP5L_3[2])
plot_grid(g, gg, labels=c("A", "B"), ncol = 1, nrow = 2)



g <- ggplot(map_circRNA_GDP7L_1, aes(start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP7L_1)),
                 binwidth = 1, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP7L_1[1], circRNA_GDP7L_1[2]) +
  labs(title="Mapping of motifs on circRNA_GDP7L_1")
gg <- ggplot(map_circRNA_GDP7L_1, aes(comp_start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP7L_1)),
                 binwidth = 1, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP7L_1[1], circRNA_GDP7L_1[2])
plot_grid(g, gg, labels=c("A", "B"), ncol = 1, nrow = 2)



g <- ggplot(map_circRNA_GDP8L_1, aes(start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_1)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_1[1], circRNA_GDP8L_1[2]) +
  labs(title="Mapping of motifs on circRNA_GDP8L_1")
gg <- ggplot(map_circRNA_GDP8L_1, aes(comp_start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_1)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_1[1], circRNA_GDP8L_1[2])
plot_grid(g, gg, labels=c("A", "B"), ncol = 1, nrow = 2)


g <- ggplot(map_circRNA_GDP8L_2, aes(start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_2)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_2[1], circRNA_GDP8L_2[2]) +
  labs(title="Mapping of motifs on circRNA_GDP8L_2")
gg <- ggplot(map_circRNA_GDP8L_2, aes(comp_start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_2)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_2[1], circRNA_GDP8L_2[2])
plot_grid(g, gg, labels=c("A", "B"), ncol = 1, nrow = 2)


g <- ggplot(map_circRNA_GDP8L_3, aes(start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_3)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_3[1], circRNA_GDP8L_3[2]) +
  labs(title="Mapping of motifs on circRNA_GDP8L_3")
gg <- ggplot(map_circRNA_GDP8L_3, aes(comp_start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_3)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_3[1], circRNA_GDP8L_3[2])
plot_grid(g, gg, labels=c("A", "B"), ncol = 1, nrow = 2)


g <- ggplot(map_circRNA_GDP8L_4, aes(start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_4)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_4[1], circRNA_GDP8L_4[2]) +
  labs(title="Mapping of motifs on circRNA_GDP8L_4")
gg <- ggplot(map_circRNA_GDP8L_4, aes(comp_start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_4)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_4[1], circRNA_GDP8L_4[2])
plot_grid(g, gg, labels=c("A", "B"), ncol = 1, nrow = 2)


g <- ggplot(map_circRNA_GDP8L_5, aes(start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_5)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_5[1], circRNA_GDP8L_5[2]) +
  labs(title="Mapping of motifs on circRNA_GDP8L_5")
gg <- ggplot(map_circRNA_GDP8L_5, aes(comp_start)) + scale_fill_brewer(palette = "Spectral") + 
  geom_histogram(aes(fill=rownames(map_circRNA_GDP8L_5)),
                 binwidth = 10, 
                 col="black", 
                 size=0.1) + theme(legend.position = "none") + xlim(circRNA_GDP8L_5[1], circRNA_GDP8L_5[2])
plot_grid(g, gg, labels=c("A", "B"), ncol = 1, nrow = 2)


#--------------------- Distribution conserved motifs all k cov2 ---------------------#


cov2 <- read.table(file = 'cov2_allk.tsv', sep = '\t', header = TRUE)
filtered_cov2 <- cov2[!duplicated(cov2[,c('start','comp_start')]),]
conserved_cov2 <- subset(filtered_cov2, count == 7)
conserved_cov2$k <- as.factor(conserved_cov2$k)

g_cov2 <- ggplot(conserved_cov2, aes(start)) + scale_fill_brewer(palette = "Spectral")
g_cov2 + geom_histogram(aes(fill=k),
                   binwidth = 1000, 
                   col="black", 
                   size=0.1) +
  labs(title="Distribution on the genome of conserved motifs across the cov2-like sequences")

#------------------- Distribution conserved motifs k > 4 cov2 --------------------#

cov2 <- read.table(file = 'cov2_allk.tsv', sep = '\t', header = TRUE)
filtered_cov2 <- cov2[!duplicated(cov2[,c('start','comp_start')]),]
conserved_cov2 <- subset(filtered_cov2, count == 7)
cov2_ksup4 <- subset(conserved_cov2, k > 4)
cov2_ksup4$k <- as.factor(cov2_ksup4$k)

g_cov2sup4 <- ggplot(cov2_ksup4, aes(start)) + scale_fill_brewer(palette = "Spectral")
g_cov2sup4 + geom_histogram(aes(fill=k),
                        binwidth = 1000, 
                        col="black", 
                        size=0.1) +
  labs(title="Distribution on the genome of conserved motifs across the cov2-like sequences")

#------------------- Mapping on circRNA all k cov2-----------------------#

cov2 <- read.table(file = 'cov2_allk.tsv', sep = '\t', header = TRUE)
filtered_cov2 <- cov2[!duplicated(cov2[,c('start','comp_start')]),]
conserved_cov2 <- subset(filtered_cov2, count == 7)
conserved_cov2$k <- as.factor(conserved_cov2$k)

map_circRNA_GDP5L_1 <- subset(conserved_cov2, start > circRNA_GDP5L_1[1] & comp_start < circRNA_GDP5L_1[2])
map_circRNA_GDP5L_2 <- subset(conserved_cov2, start > circRNA_GDP5L_2[1] & comp_start < circRNA_GDP5L_2[2])
map_circRNA_GDP5L_3 <- subset(conserved_cov2, start > circRNA_GDP5L_3[1] & comp_start < circRNA_GDP5L_3[2])
map_circRNA_GDP7L_1 <- subset(conserved_cov2, start > circRNA_GDP7L_1[1] & comp_start < circRNA_GDP7L_1[2])
map_circRNA_GDP7L_2 <- subset(conserved_cov2, start > circRNA_GDP7L_2[1] & comp_start < circRNA_GDP7L_2[2])
map_circRNA_GDP8L_1 <- subset(conserved_cov2, start > circRNA_GDP8L_1[1] & comp_start < circRNA_GDP8L_1[2])
map_circRNA_GDP8L_2 <- subset(conserved_cov2, start > circRNA_GDP8L_2[1] & comp_start < circRNA_GDP8L_2[2])
map_circRNA_GDP8L_3 <- subset(conserved_cov2, start > circRNA_GDP8L_3[1] & comp_start < circRNA_GDP8L_3[2])
map_circRNA_GDP8L_4 <- subset(conserved_cov2, start > circRNA_GDP8L_4[1] & comp_start < circRNA_GDP8L_4[2])
map_circRNA_GDP8L_5 <- subset(conserved_cov2, start > circRNA_GDP8L_5[1] & comp_start < circRNA_GDP8L_5[2])

#-------------------- Compensate mutation cov2 ---------------#

cov2 <- read.table(file = 'cov2_allk.tsv', sep = '\t', header = TRUE)
conserved_cov2 <- subset(cov2, count == 7)
compensate_cov2 <- subset(conserved_cov2, count_kmer < 7)
filtered_cov2 <- compensate_cov2[!duplicated(compensate_cov2[,c('start','comp_start')]),]
filtered_cov2$k <- as.factor(filtered_cov2$k)
g_cov2_compensate <- ggplot(filtered_cov2, aes(start)) + scale_fill_brewer(palette = "Spectral")
g_cov2_compensate + geom_histogram(aes(fill=k),
                        binwidth = 1000, 
                        col="black", 
                        size=0.1) +
  labs(title="Distribution on the genome of compensate mutation across the cov2-like sequences")

#------------------- Mapping on circRNA compensate mutation cov2-----------------------#

cov2 <- read.table(file = 'cov2_allk.tsv', sep = '\t', header = TRUE)
conserved_cov2 <- subset(cov2, count == 7)
compensate_cov2 <- subset(conserved_cov2, count_kmer < 7)
filtered_cov2 <- compensate_cov2[!duplicated(compensate_cov2[,c('start','comp_start')]),]
filtered_cov2$k <- as.factor(filtered_cov2$k)

map_circRNA_GDP5L_1 <- subset(filtered_cov2, start > circRNA_GDP5L_1[1] & comp_start < circRNA_GDP5L_1[2])
map_circRNA_GDP5L_2 <- subset(filtered_cov2, start > circRNA_GDP5L_2[1] & comp_start < circRNA_GDP5L_2[2])
map_circRNA_GDP5L_3 <- subset(filtered_cov2, start > circRNA_GDP5L_3[1] & comp_start < circRNA_GDP5L_3[2])
map_circRNA_GDP7L_1 <- subset(filtered_cov2, start > circRNA_GDP7L_1[1] & comp_start < circRNA_GDP7L_1[2])
map_circRNA_GDP7L_2 <- subset(filtered_cov2, start > circRNA_GDP7L_2[1] & comp_start < circRNA_GDP7L_2[2])
map_circRNA_GDP8L_1 <- subset(filtered_cov2, start > circRNA_GDP8L_1[1] & comp_start < circRNA_GDP8L_1[2])
map_circRNA_GDP8L_2 <- subset(filtered_cov2, start > circRNA_GDP8L_2[1] & comp_start < circRNA_GDP8L_2[2])
map_circRNA_GDP8L_3 <- subset(filtered_cov2, start > circRNA_GDP8L_3[1] & comp_start < circRNA_GDP8L_3[2])
map_circRNA_GDP8L_4 <- subset(filtered_cov2, start > circRNA_GDP8L_4[1] & comp_start < circRNA_GDP8L_4[2])
map_circRNA_GDP8L_5 <- subset(filtered_cov2, start > circRNA_GDP8L_5[1] & comp_start < circRNA_GDP8L_5[2])


#--------------------- Distribution conserved motifs all k cov ---------------------#

cov <- read.table(file = 'cov_allk.tsv', sep = '\t', header = TRUE)
filtered_cov <- cov[!duplicated(cov[,c('start','comp_start')]),]
conserved_cov <- subset(filtered_cov, count == 44)
conserved_cov$k <- as.factor(conserved_cov$k)

g_cov <- ggplot(conserved_cov, aes(start)) + scale_fill_brewer(palette = "Spectral")
g_cov + geom_histogram(aes(fill=k),
                        binwidth = 1000, 
                        col="black", 
                        size=0.1) +
  labs(title="Distribution on the genome of conserved motifs across the cov-like sequences")


#------------------- Distribution conserved motifs k > 4 cov --------------------#

cov <- read.table(file = 'cov_allk.tsv', sep = '\t', header = TRUE)
filtered_cov <- cov[!duplicated(cov[,c('start','comp_start')]),]
conserved_cov <- subset(filtered_cov, count == 44)
cov_ksup4 <- subset(conserved_cov, k > 4)
cov_ksup4$k <- as.factor(cov_ksup4$k)

g_covsup4 <- ggplot(cov_ksup4, aes(start)) + scale_fill_brewer(palette = "Spectral")
g_covsup4 + geom_histogram(aes(fill=k),
                            binwidth = 1000, 
                            col="black", 
                            size=0.1) +
  labs(title="Distribution on the genome of conserved motifs across the cov-like sequences")



#------------------- Mapping on circRNA all k cov2-----------------------#

cov <- read.table(file = 'cov_allk.tsv', sep = '\t', header = TRUE)
filtered_cov <- cov[!duplicated(cov[,c('start','comp_start')]),]
conserved_cov <- subset(filtered_cov, count == 44)
conserved_cov$k <- as.factor(conserved_cov$k)

map_circRNA_GDP5L_1 <- subset(conserved_cov, start > circRNA_GDP5L_1[1] & comp_start < circRNA_GDP5L_1[2])
map_circRNA_GDP5L_2 <- subset(conserved_cov, start > circRNA_GDP5L_2[1] & comp_start < circRNA_GDP5L_2[2])
map_circRNA_GDP5L_3 <- subset(conserved_cov, start > circRNA_GDP5L_3[1] & comp_start < circRNA_GDP5L_3[2])
map_circRNA_GDP7L_1 <- subset(conserved_cov, start > circRNA_GDP7L_1[1] & comp_start < circRNA_GDP7L_1[2])
map_circRNA_GDP7L_2 <- subset(conserved_cov, start > circRNA_GDP7L_2[1] & comp_start < circRNA_GDP7L_2[2])
map_circRNA_GDP8L_1 <- subset(conserved_cov, start > circRNA_GDP8L_1[1] & comp_start < circRNA_GDP8L_1[2])
map_circRNA_GDP8L_2 <- subset(conserved_cov, start > circRNA_GDP8L_2[1] & comp_start < circRNA_GDP8L_2[2])
map_circRNA_GDP8L_3 <- subset(conserved_cov, start > circRNA_GDP8L_3[1] & comp_start < circRNA_GDP8L_3[2])
map_circRNA_GDP8L_4 <- subset(conserved_cov, start > circRNA_GDP8L_4[1] & comp_start < circRNA_GDP8L_4[2])
map_circRNA_GDP8L_5 <- subset(conserved_cov, start > circRNA_GDP8L_5[1] & comp_start < circRNA_GDP8L_5[2])

#-------------------- Compensate mutation cov ---------------#

cov <- read.table(file = 'cov_allk.tsv', sep = '\t', header = TRUE)
conserved_cov <- subset(cov, count == 44)
compensate_cov <- subset(conserved_cov, count_kmer < 44)
filtered_cov <- compensate_cov[!duplicated(compensate_cov[,c('start','comp_start')]),]
filtered_cov$k <- as.factor(filtered_cov$k)
g_cov_compensate <- ggplot(filtered_cov, aes(start)) + scale_fill_brewer(palette = "Spectral")
g_cov_compensate + geom_histogram(aes(fill=k),
                                   binwidth = 1000, 
                                   col="black", 
                                   size=0.1) +
  labs(title="Distribution on the genome of compensate mutation across the cov-like sequences")

#------------------- Mapping on circRNA compensate mutation cov-----------------------#

cov <- read.table(file = 'cov_allk.tsv', sep = '\t', header = TRUE)
conserved_cov <- subset(cov, count == 44)
compensate_cov <- subset(conserved_cov, count_kmer < 44)
filtered_cov <- compensate_cov[!duplicated(compensate_cov[,c('start','comp_start')]),]
filtered_cov$k <- as.factor(filtered_cov$k)

map_circRNA_GDP5L_1 <- subset(filtered_cov, start > circRNA_GDP5L_1[1] & comp_start < circRNA_GDP5L_1[2])
map_circRNA_GDP5L_2 <- subset(filtered_cov, start > circRNA_GDP5L_2[1] & comp_start < circRNA_GDP5L_2[2])
map_circRNA_GDP5L_3 <- subset(filtered_cov, start > circRNA_GDP5L_3[1] & comp_start < circRNA_GDP5L_3[2])
map_circRNA_GDP7L_1 <- subset(filtered_cov, start > circRNA_GDP7L_1[1] & comp_start < circRNA_GDP7L_1[2])
map_circRNA_GDP7L_2 <- subset(filtered_cov, start > circRNA_GDP7L_2[1] & comp_start < circRNA_GDP7L_2[2])
map_circRNA_GDP8L_1 <- subset(filtered_cov, start > circRNA_GDP8L_1[1] & comp_start < circRNA_GDP8L_1[2])
map_circRNA_GDP8L_2 <- subset(filtered_cov, start > circRNA_GDP8L_2[1] & comp_start < circRNA_GDP8L_2[2])
map_circRNA_GDP8L_3 <- subset(filtered_cov, start > circRNA_GDP8L_3[1] & comp_start < circRNA_GDP8L_3[2])
map_circRNA_GDP8L_4 <- subset(filtered_cov, start > circRNA_GDP8L_4[1] & comp_start < circRNA_GDP8L_4[2])
map_circRNA_GDP8L_5 <- subset(filtered_cov, start > circRNA_GDP8L_5[1] & comp_start < circRNA_GDP8L_5[2])