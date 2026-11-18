# MAG x run heatmap of CoverM relative abundance, runs ordered by site category and site,
# MAGs ordered by the category where they are most abundant (then by that abundance).
# Usage: Rscript heatmap.R <outdir> <out.png>
#   <outdir>: pipeline --outdir (07_mags/coverm_genome.tsv, metadata/samples.tsv,
#             findings/mag_presence.tsv from scripts/findings/findings.py)
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
cov <- read.delim(file.path(out, "07_mags", "coverm_genome.tsv"), check.names = FALSE)
sheet <- read.delim(file.path(out, "metadata", "samples.tsv"))
pres <- read.delim(file.path(out, "findings", "mag_presence.tsv"))

cats <- c("pre_treatment", "plant", "creek_downstream", "creek_upstream", "control")
labels <- c(pre_treatment = "Pre-treatment", plant = "Plant steps", creek_downstream = "Creek downstream",
            creek_upstream = "Creek upstream", control = "Ctrl")
site_cat <- c(HOSP = "pre_treatment", RES = "pre_treatment", INF = "pre_treatment", PCI = "pre_treatment",
              PCE = "pre_treatment", ATE = "plant", FCE = "plant", UV = "plant", DS_A = "creek_downstream",
              DS_B = "creek_downstream", UP_A = "creek_upstream", UP_B = "creek_upstream", MNT_A = "control",
              MNT_B = "control", UW_A = "control", UW_B = "control")
site_order <- names(site_cat)
sheet$category <- factor(site_cat[sheet$site], levels = cats)
sheet <- sheet[order(sheet$category, match(sheet$site, site_order), sheet$location, sheet$timepoint), ]

ra <- cov[cov$Genome != "unmapped", c("Genome", grep("Relative Abundance", names(cov), value = TRUE))]
names(ra) <- c("mag", sub("_1.fastq.gz Relative Abundance \\(%\\)", "", names(ra)[-1]))
m <- as.matrix(ra[, sheet$run]); rownames(m) <- ra$mag

# MAG order: category of highest mean abundance, then that mean (descending)
means <- sapply(cats, function(c) rowMeans(m[, sheet$run[sheet$category == c], drop = FALSE]))
home <- factor(cats[max.col(means, ties.method = "first")], levels = cats)
ord <- order(home, -apply(means, 1, max))
long <- data.frame(mag = factor(rep(rownames(m), ncol(m)), levels = rev(rownames(m)[ord])),
                   run = factor(rep(colnames(m), each = nrow(m)), levels = sheet$run),
                   ra = as.vector(m))
long$category <- sheet$category[match(long$run, sheet$run)]
long$value <- log10(pmax(long$ra, 1e-3))
arg <- pres$mag[pres$resistance_genes != ""]

p <- ggplot(long, aes(run, mag, fill = value)) +
  geom_raster() +
  facet_grid(. ~ category, scales = "free_x", space = "free_x", labeller = as_labeller(labels)) +
  scale_fill_gradientn(colours = c("#F4F8F7", "#A8D5CE", "#2A9D8F", "#0B4F5C", "#1B2B30"),
                       limits = c(-3, 1.5), oob = scales::squish,
                       breaks = c(-3, -2, -1, 0, 1), labels = c("<0.001", "0.01", "0.1", "1", "10"),
                       name = "Relative\nabundance (%)") +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = sprintf("%d runs (sites in order HOSP, RES, INF, PCI, PCE | ATE, FCE, UV | DS | UP | controls)", nrow(sheet)),
       y = sprintf("%d catalog MAGs (grouped by where they are most abundant)", nrow(m))) +
  theme_minimal(base_size = 13) +
  theme(panel.spacing = unit(4, "pt"), strip.text = element_text(face = "bold", colour = "#0B4F5C"),
        panel.grid = element_blank(), legend.position = "right")
# mark the MAGs that carry resistance genes (red > at the left edge)
ticks <- data.frame(mag = factor(arg, levels = levels(long$mag)),
                    category = factor("pre_treatment", levels = cats))   # left edge, first panel only
p <- p + geom_point(data = ticks, aes(x = -Inf, y = mag), inherit.aes = FALSE, shape = 62, size = 2.2,
                    colour = "#C8553D")
ggsave(args[2], p, width = 12, height = 6.2, dpi = 200, bg = "white")
cat("wrote", args[2], "\n")
