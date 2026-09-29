# Plasmid yield against sequencing depth (stage 9 depth experiment).
# Usage: Rscript depth_plot.R <09_depth_experiment.tsv> <output prefix>
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
d <- read.delim(args[1], stringsAsFactors = FALSE)
d$million_pairs <- d$pairs / 1e6
long <- rbind(
  data.frame(run = d$run, million_pairs = d$million_pairs, level = d$level,
             measure = "plasmid contigs", value = d$plasmid_contigs),
  data.frame(run = d$run, million_pairs = d$million_pairs, level = d$level,
             measure = "plasmid bp", value = d$plasmid_bp)
)
p <- ggplot(long, aes(million_pairs, value, colour = run)) +
  geom_line() + geom_point() +
  facet_wrap(~measure, scales = "free_y") +
  labs(x = "host-depleted read pairs (millions)", y = NULL,
       title = "geNomad plasmid yield vs depth") +
  theme_bw()
ggsave(paste0(args[2], ".pdf"), p, width = 9, height = 4)
ggsave(paste0(args[2], ".png"), p, width = 9, height = 4, dpi = 150)
cat("wrote", paste0(args[2], ".{pdf,png}"), "\n")
