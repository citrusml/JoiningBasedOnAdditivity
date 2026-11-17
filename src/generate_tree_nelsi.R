#!/usr/bin/env Rscript
#
# generate_tree_nelsi.R - Generate synthetic phylogenetic trees with Yule topology
# and autocorrelated branch lengths using the NELSI R package.
#
# Usage:
#   Rscript bin/generate_tree_nelsi.R --taxa 32 --scale 1.0 --seed 42 \
#     --model nelsi_kishino --v 0.3 --initial_rate 1.0 \
#     --scale_mode root_to_tip --outtree tree.nwk
#
# Models:
#   - nelsi_kishino     (alias: kishino, nelsi):       NELSI::simulate.autocor.kishino (Kishino et al. 2001)
#   - nelsi_thorne      (alias: thorne):               NELSI::simulate.autocor.thorne (Thorne et al. 1998)
#   - nelsi_uncor_lnorm (alias: uncor_lnorm, ucln):    NELSI::simulate.uncor.lnorm (Uncorrelated lognormal)
#
# Output:
#   Standard Newick format without internal node labels, leaf names T1..TN,
#   compatible with INDELible, AliSim, and DendroPy.

print_usage <- function() {
  cat(paste0(
    "Usage: Rscript bin/generate_tree_nelsi.R [options] --outtree <file.nwk>\n",
    "Options:\n",
    "  --taxa <int>           Number of terminal taxa (default: 32, min: 2)\n",
    "  --scale <float>        Evolutionary distance scaling factor D (default: 1.0)\n",
    "  --seed <int>           Random seed for reproducibility (optional)\n",
    "  --model <str>          NELSI model: nelsi_kishino (default), nelsi_thorne, or nelsi_uncor_lnorm\n",
    "  --v <float>            NELSI autocorrelation variance parameter v (default: 0.3)\n",
    "  --initial_rate <float> NELSI initial root substitution rate (default: 1.0)\n",
    "  --mean_log <float>     NELSI uncor.lnorm mean.log parameter (default: 0.0)\n",
    "  --sd_log <float>       NELSI uncor.lnorm sd.log parameter (default: 0.5)\n",
    "  --scale_mode <str>     Scaling mode: root_to_tip (default), mean_branch, or raw\n",
    "  --lba_ratio <float>    Long-branch attraction ratio on T1 and TN (default: 1.0)\n",
    "  --min_length <float>   Minimum branch length floor (default: 1e-6)\n",
    "  --outtree <path>       Output Newick tree filepath (required)\n",
    "  --help, -h             Show this help message and exit\n"
  ))
}

parse_cli_args <- function(args) {
  opts <- list(
    taxa = 32L,
    scale = 1.0,
    seed = NULL,
    model = "nelsi_kishino",
    v = 0.3,
    initial_rate = 1.0,
    mean_log = 0.0,
    sd_log = 0.5,
    scale_mode = "root_to_tip",
    lba_ratio = 1.0,
    min_length = 1e-6,
    outtree = NULL,
    help = FALSE
  )

  i <- 1L
  n <- length(args)
  while (i <= n) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      opts$help <- TRUE
      return(opts)
    }

    # Support --key=value as well as --key value
    if (grepl("^--[^=]+=", arg)) {
      parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      val <- paste(parts[-1], collapse = "=")
      i <- i + 1L
    } else if (startsWith(arg, "--")) {
      key <- sub("^--", "", arg)
      if (i + 1L > n) {
        stop(sprintf("Missing value for argument --%s", key))
      }
      val <- args[[i + 1L]]
      i <- i + 2L
    } else {
      stop(sprintf("Unexpected positional argument: %s", arg))
    }

    if (key == "taxa") {
      opts$taxa <- as.integer(val)
    } else if (key == "scale") {
      opts$scale <- as.numeric(val)
    } else if (key == "seed") {
      if (tolower(val) != "none" && val != "") {
        opts$seed <- as.integer(val)
      }
    } else if (key == "model") {
      opts$model <- tolower(val)
    } else if (key %in% c("v", "nelsi_v")) {
      opts$v <- as.numeric(val)
    } else if (key %in% c("initial_rate", "nelsi_initial_rate")) {
      opts$initial_rate <- as.numeric(val)
    } else if (key %in% c("mean_log", "nelsi_mean_log")) {
      opts$mean_log <- as.numeric(val)
    } else if (key %in% c("sd_log", "nelsi_sd_log")) {
      opts$sd_log <- as.numeric(val)
    } else if (key %in% c("scale_mode", "nelsi_scale_mode")) {
      opts$scale_mode <- tolower(val)
    } else if (key == "lba_ratio") {
      opts$lba_ratio <- as.numeric(val)
    } else if (key == "min_length") {
      opts$min_length <- as.numeric(val)
    } else if (key == "outtree") {
      opts$outtree <- val
    } else {
      stop(sprintf("Unknown argument: --%s", key))
    }
  }
  return(opts)
}

args <- commandArgs(trailingOnly = TRUE)
opts <- tryCatch(
  parse_cli_args(args),
  error = function(e) {
    message(sprintf("Argument error: %s", e$message))
    print_usage()
    quit(status = 1)
  }
)

if (opts$help) {
  print_usage()
  quit(status = 0)
}

if (is.null(opts$outtree) || opts$outtree == "") {
  message("Error: --outtree is required.")
  print_usage()
  quit(status = 1)
}
if (is.na(opts$taxa) || opts$taxa < 2L) {
  stop(sprintf("Taxa count must be at least 2, got %s", as.character(opts$taxa)))
}
if (is.na(opts$scale) || opts$scale <= 0.0) {
  stop(sprintf("Scale must be strictly positive, got %s", as.character(opts$scale)))
}
if (is.na(opts$v) || opts$v < 0.0) {
  stop(sprintf("NELSI autocorrelation parameter v must be non-negative, got %s", as.character(opts$v)))
}
if (is.na(opts$initial_rate) || opts$initial_rate <= 0.0) {
  stop(sprintf("NELSI initial_rate must be strictly positive, got %s", as.character(opts$initial_rate)))
}
if (is.na(opts$mean_log)) {
  stop("NELSI mean_log must be a valid finite number.")
}
if (is.na(opts$sd_log) || opts$sd_log < 0.0) {
  stop(sprintf("NELSI sd_log must be non-negative, got %s", as.character(opts$sd_log)))
}
if (!(opts$scale_mode %in% c("root_to_tip", "mean_branch", "raw"))) {
  stop(sprintf("Invalid --scale_mode '%s'. Must be root_to_tip, mean_branch, or raw.", opts$scale_mode))
}

# Check required R packages strictly without fallback
if (!requireNamespace("ape", quietly = TRUE)) {
  message("Error: R package 'ape' is required but not installed.")
  message("Install in R with: install.packages('ape', repos='https://cloud.r-project.org')")
  quit(status = 1)
}
if (!requireNamespace("NELSI", quietly = TRUE)) {
  message("Error: R package 'NELSI' is required but not installed.")
  message("Install dependencies and NELSI in R with:")
  message("  install.packages(c('ape', 'castor', 'digest', 'remotes'), repos='https://cloud.r-project.org')")
  message("  remotes::install_github('sebastianduchene/NELSI')")
  quit(status = 1)
}

suppressPackageStartupMessages({
  library(ape)
  library(NELSI)
})

if (!is.null(opts$seed) && !is.na(opts$seed)) {
  set.seed(opts$seed)
}

# 1. Generate rooted Yule chronogram (pure-birth process, death = 0)
chronogram <- ape::rphylo(n = opts$taxa, birth = 1.0, death = 0.0)
chronogram <- ape::reorder.phylo(chronogram, order = "cladewise")
chronogram$tip.label <- paste0("T", seq_len(opts$taxa))
chronogram$node.label <- NULL

# Normalize time tree root height to 1.0
root_height <- max(ape::node.depth.edgelength(chronogram))
if (!is.finite(root_height) || root_height <= 0.0) {
  stop("Generated Yule chronogram has non-positive root height.")
}
chronogram$edge.length <- chronogram$edge.length / root_height

# 2. Simulate rates and branch lengths with NELSI
nelsi_params <- list(initial.rate = opts$initial_rate, v = opts$v)

if (opts$model %in% c("nelsi_kishino", "kishino", "nelsi", "nelsi_yule")) {
  sim_res <- NELSI::simulate.rate(chronogram, NELSI::simulate.autocor.kishino, params = nelsi_params)
  canonical_model <- "nelsi_kishino"
} else if (opts$model %in% c("nelsi_thorne", "thorne")) {
  sim_res <- NELSI::simulate.rate(chronogram, NELSI::simulate.autocor.thorne, params = nelsi_params)
  canonical_model <- "nelsi_thorne"
} else if (opts$model %in% c("nelsi_uncor_lnorm", "uncor_lnorm", "ucln")) {
  uncor_params <- list(mean.log = opts$mean_log, sd.log = opts$sd_log)
  sim_res <- NELSI::simulate.rate(chronogram, NELSI::simulate.uncor.lnorm, params = uncor_params)
  canonical_model <- "nelsi_uncor_lnorm"
} else {
  stop(sprintf("Unsupported NELSI model '%s'. Must be 'nelsi_kishino', 'nelsi_thorne', or 'nelsi_uncor_lnorm'.", opts$model))
}

phylogram <- sim_res$phylogram
phylogram$node.label <- NULL
phylogram$root.edge <- NULL

# Guard against any zero/non-finite edge lengths prior to scaling
if (any(!is.finite(phylogram$edge.length))) {
  stop("NELSI produced non-finite branch lengths.")
}
phylogram$edge.length <- pmax(phylogram$edge.length, opts$min_length)

# 3. Scale branch lengths to target evolutionary distance D (--scale)
if (opts$scale_mode == "root_to_tip") {
  tip_depths <- ape::node.depth.edgelength(phylogram)[seq_len(opts$taxa)]
  mean_rtt <- mean(tip_depths)
  if (!is.finite(mean_rtt) || mean_rtt <= 0.0) {
    stop("Mean root-to-tip distance is non-positive; cannot rescale.")
  }
  phylogram$edge.length <- phylogram$edge.length * (opts$scale / mean_rtt)
} else if (opts$scale_mode == "mean_branch") {
  mean_bl <- mean(phylogram$edge.length)
  if (!is.finite(mean_bl) || mean_bl <= 0.0) {
    stop("Mean branch length is non-positive; cannot rescale.")
  }
  phylogram$edge.length <- phylogram$edge.length * (opts$scale / mean_bl)
} else {
  phylogram$edge.length <- phylogram$edge.length * opts$scale
}

# 4. Optional LBA ratio on terminal branches T1 and TN
if (opts$lba_ratio > 1.0) {
  lba_tips <- c("T1", paste0("T", opts$taxa))
  lba_tip_indices <- which(phylogram$tip.label %in% lba_tips)
  lba_edges <- which(phylogram$edge[, 2] %in% lba_tip_indices)
  phylogram$edge.length[lba_edges] <- phylogram$edge.length[lba_edges] * opts$lba_ratio
}

# 5. Final minimum branch length floor for INDELible / AliSim compatibility
phylogram$edge.length <- pmax(phylogram$edge.length, opts$min_length)

# 6. Ensure destination directory exists and write Newick file
#    Convert any scientific notation branch lengths (e.g., 1.25e-05) to fixed-point
#    decimals because INDELible V1.03 rejects 'e' in guide tree branch lengths.
out_dir <- dirname(opts$outtree)
if (out_dir != "." && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
nwk_str <- ape::write.tree(phylogram, file = "")
sci_matches <- gregexpr(":[0-9]+(\\.[0-9]+)?[eE][+-]?[0-9]+", nwk_str, perl = TRUE)
matched_tokens <- regmatches(nwk_str, sci_matches)[[1]]
if (length(matched_tokens) > 0L) {
  vals <- as.numeric(sub("^:", "", matched_tokens))
  regmatches(nwk_str, sci_matches) <- list(sprintf(":%.10f", vals))
}
writeLines(nwk_str, con = opts$outtree)

cat(sprintf(
  "Generated %s tree [NELSI] (%d taxa, scale=%.6g, v=%.6g, initial_rate=%.6g, scale_mode=%s, seed=%s) -> %s\n",
  canonical_model,
  opts$taxa,
  opts$scale,
  opts$v,
  opts$initial_rate,
  opts$scale_mode,
  ifelse(is.null(opts$seed) || is.na(opts$seed), "None", as.character(opts$seed)),
  opts$outtree
))

