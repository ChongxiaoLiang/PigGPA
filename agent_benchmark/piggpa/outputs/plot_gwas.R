
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

gwas <- fread("/workspace/output/gcta-gwas/lmd_gcta.mlma")
setnames(gwas, old=c("Chr","bp","p"), new=c("CHR","BP","P"))
gwas <- gwas[!is.na(CHR) & !is.na(BP) & !is.na(P) & P > 0 & P <= 1]
gwas <- gwas[order(CHR, BP)]

# Cumulative position for Manhattan
chr_offsets <- gwas[, .(max_pos = max(BP)), by=CHR]
chr_offsets[, cum_offset := cumsum(as.numeric(c(0, max_pos[-.N])))]
gwas <- merge(gwas, chr_offsets[, .(CHR, cum_offset)], by="CHR")
gwas[, pos_cum := BP + cum_offset]
gwas[, chr_color := ifelse(as.integer(CHR) %% 2 == 1, "#1f77b4", "#ff7f0e")]

# Genomic inflation
chisq <- qchisq(gwas$P, df=1, lower.tail=FALSE)
lambda <- median(chisq, na.rm=TRUE) / qchisq(0.5, df=1)
cat("Lambda:", round(lambda, 4), "\n")

# Sample for plotting speed (80K points)
set.seed(42)
gwas_plot <- gwas[sample(.N, min(.N, 80000))]

# Manhattan
chr_centers <- gwas_plot[, .(center = mean(pos_cum)), by=CHR]
PALETTE_A <- c('#4285B0', '#86A7CB', '#8E3E78', '#B380A4', '#EC8528', '#8D80B3', '#EAC94D', '#1A9899')

p1 <- ggplot(gwas_plot, aes(x=pos_cum, y=-log10(P), color=chr_color)) +
  geom_point(size=0.4, alpha=0.7) + scale_color_identity() +
  geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="#E6550D", linewidth=0.5) +
  geom_hline(yintercept=-log10(1e-5), linetype="dotted", color="#2CA02C", linewidth=0.4) +
  scale_x_continuous(breaks=chr_centers$center, labels=chr_centers$CHR, expand=c(0.01,0.01)) +
  labs(title=paste0("GCTA MLMA: Loin Muscle Depth (lambda=", round(lambda,3), ")"), x="Chromosome", y=expression(-log[10](italic(P)))) +
  theme_bw(base_size=11) + theme(panel.grid=element_blank(), legend.position="none")

ggsave("/workspace/output/gcta-gwas/lmd_manhattan.png", p1, width=12, height=6, dpi=300)
ggsave("/workspace/output/gcta-gwas/lmd_manhattan.pdf", p1, width=12, height=6)

# QQ plot  
gwas[, `:=`(observed=-log10(sort(P)), expected=-log10(ppoints(.N)))]
qq_sample <- gwas[seq(1, .N, length.out=min(.N, 50000))]
p2 <- ggplot(qq_sample, aes(x=expected, y=observed)) +
  geom_point(size=0.5, alpha=0.5, color="#1f77b4") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="#E6550D", linewidth=0.7) +
  labs(title=paste0("QQ Plot (lambda=", round(lambda,3), ")"),
       x=expression(Expected~-log[10](italic(P))), y=expression(Observed~-log[10](italic(P)))) +
  theme_bw(base_size=11) + theme(panel.grid=element_blank()) + coord_fixed(ratio=1)

ggsave("/workspace/output/gcta-gwas/lmd_qq.png", p2, width=6, height=6, dpi=300)
ggsave("/workspace/output/gcta-gwas/lmd_qq.pdf", p2, width=6, height=6)

cat("Plots saved successfully!\n")

# Save suggestive SNPs
sig <- gwas[P < 1e-5]
sig <- sig[order(P)]
fwrite(sig, "/workspace/output/gcta-gwas/lmd_suggestive.csv")
cat("Suggestive SNPs saved:", nrow(sig), "\n")
for (i in 1:min(nrow(sig), 10)) {
  cat(sprintf("  %s  chr%s:%d  p=%.2e\n", sig$SNP[i], sig$CHR[i], sig$BP[i], sig$P[i]))
}
