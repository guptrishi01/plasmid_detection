# nonpareil_summary.R -- side analysis: one row per Nonpareil curve (.npo).
#   Rscript nonpareil_summary.R <out.tsv> <file.npo> [<file.npo> ...]
# Columns: run, kappa (raw redundancy at full effort), coverage (modelled C),
# effort_bp (LR, bp sequenced), effort_bp_for_95 (LRstar), diversity (Nd), model_r.
# Needs the Nonpareil R package (pd-assembly).
suppressMessages(library(Nonpareil))

args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
files <- args[-1]
stopifnot(length(files) > 0)

rows <- lapply(files, function(f) {
  s <- summary(Nonpareil.curve(f, plot = FALSE))
  data.frame(
    run = sub("\\.npo$", "", basename(f)),
    kappa = unname(s["kappa"]),
    coverage = unname(s["C"]),
    effort_bp = unname(s["LR"]),
    effort_bp_for_95 = unname(s["LRstar"]),
    diversity = unname(s["diversity"]),
    model_r = unname(s["modelR"])
  )
})
write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE, row.names = FALSE)
