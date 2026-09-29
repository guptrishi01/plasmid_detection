# Stage 9 statistics on a resistome count matrix (gene x run).
# Usage: Rscript stats.R <stats_samples.tsv> <resistome counts tsv> <output dir>
# Parameters come from the environment (config/config.yaml via scripts/lib/common.sh):
#   PARAMS_STATS_COMPARISONS  "groupA:groupB ..."   PARAMS_STATS_METHOD  aldex2 | ancombc2
#   PARAMS_STATS_MIN_GROUP_SIZE, PARAMS_STATS_ALPHA
# Controls (site_group == "control") are a reference: summarised, never tested or subtracted.
suppressPackageStartupMessages({ library(vegan) })
args <- commandArgs(trailingOnly = TRUE)
meta <- read.delim(args[1], stringsAsFactors = FALSE)
counts <- read.delim(args[2], row.names = 1, check.names = FALSE)
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
comparisons <- strsplit(Sys.getenv("PARAMS_STATS_COMPARISONS"), " ")[[1]]
method <- Sys.getenv("PARAMS_STATS_METHOD", "aldex2")
min_n <- as.integer(Sys.getenv("PARAMS_STATS_MIN_GROUP_SIZE", "3"))
alpha <- as.numeric(Sys.getenv("PARAMS_STATS_ALPHA", "0.05"))
log_lines <- character()
note <- function(...) { msg <- paste0(...); message(msg); log_lines <<- c(log_lines, msg) }

meta <- meta[meta$run %in% colnames(counts), ]
counts <- as.matrix(counts[, meta$run, drop = FALSE])
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
note("genes with reads: ", nrow(counts), "; runs: ", ncol(counts))

# CLR (pseudocount 0.5) for description and the Aitchison PERMANOVA
clr <- apply(counts + 0.5, 2, function(x) log(x) - mean(log(x)))
write.table(data.frame(gene = rownames(clr), clr, check.names = FALSE),
            file.path(out, "clr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
grp_means <- sapply(split(meta$run, meta$site_group), function(r) rowMeans(clr[, r, drop = FALSE]))
write.table(data.frame(gene = rownames(clr), grp_means, check.names = FALSE),
            file.path(out, "clr_mean_by_site_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

test_meta <- meta[meta$site_group != "control", ]

# PERMANOVA (exploratory)
if (nrow(test_meta) >= 4 && nrow(counts) >= 2) {
  d <- dist(t(clr[, test_meta$run, drop = FALSE]))
  terms <- c("location", "timepoint", "site_group")
  terms <- terms[sapply(terms, function(t) length(unique(test_meta[[t]])) > 1)]
  if (length(terms)) {
    test_meta$timepoint <- factor(test_meta$timepoint)
    fit <- adonis2(as.formula(paste("d ~", paste(terms, collapse = " + "))),
                   data = test_meta, by = "margin", permutations = 9999)
    capture.output(fit, file = file.path(out, "permanova.txt"))
    note("PERMANOVA (exploratory) on ", nrow(test_meta), " runs: ", paste(terms, collapse = " + "))
  } else note("PERMANOVA skipped: no term varies across the non-control runs")
} else note("PERMANOVA skipped: ", nrow(test_meta), " non-control runs (need >= 4)")

run_test <- function(label, a_runs, b_runs) {
  if (length(a_runs) < min_n || length(b_runs) < min_n) {
    note(label, ": skipped (n = ", length(a_runs), " vs ", length(b_runs), "; need >= ", min_n, ")")
    return(invisible(NULL))
  }
  x <- counts[, c(a_runs, b_runs), drop = FALSE]
  cond <- c(rep("a", length(a_runs)), rep("b", length(b_runs)))
  if (method == "aldex2") {
    suppressPackageStartupMessages(library(ALDEx2))
    res <- aldex(x, cond, mc.samples = 128, test = "t", effect = TRUE, denom = "all")
    res$gene <- rownames(res)
    res$significant_BH <- res$wi.eBH < alpha   # ALDEx2's eBH columns are BH-adjusted
  } else {
    suppressPackageStartupMessages(library(ANCOMBC))
    # ANCOMBC >= 2.x takes a taxa-by-sample matrix + meta_data directly (no TSE/phyloseq)
    # two columns: a one-column data.frame collapses to a vector inside ANCOMBC's checks
    md <- data.frame(group = cond, run = colnames(x), row.names = colnames(x),
                     stringsAsFactors = FALSE)
    fit <- ancombc2(data = x, taxa_are_rows = TRUE, meta_data = md, fix_formula = "group",
                    p_adj_method = "BH", prv_cut = 0, alpha = alpha, group = "group",
                    verbose = FALSE)
    res <- fit$res
  }
  f <- file.path(out, paste0(gsub("[^A-Za-z0-9_]+", "_", label), ".", method, ".tsv"))
  write.table(res, f, sep = "\t", quote = FALSE, row.names = FALSE)
  note(label, ": ", method, " on ", length(a_runs), " vs ", length(b_runs), " -> ", basename(f))
}

for (cmp in comparisons) {
  g <- strsplit(cmp, ":")[[1]]
  for (loc in unique(test_meta$location)) {
    m <- test_meta[test_meta$location == loc, ]
    run_test(paste(cmp, "within", loc), m$run[m$site_group == g[1]], m$run[m$site_group == g[2]])
  }
  for (tp in unique(test_meta$timepoint)) {
    m <- test_meta[test_meta$timepoint == tp, ]
    run_test(paste(cmp, "within timepoint", tp), m$run[m$site_group == g[1]], m$run[m$site_group == g[2]])
  }
}
writeLines(log_lines, file.path(out, "stats_log.txt"))
