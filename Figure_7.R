##############################################################################
## Code to generate Figure 7 of the manuscript:
# "Scale- and genotype-dependent spatial concordance of Mycobacterium bovis
#  in cattle and badgers"
#
# Authors: Assel Akhmetova, Liliana C. M. Salvador, Katarina Oravcova,
# Eleanor Presho, Carl McCormick, Suzan Thompson, Lorraine Wright,
# Fraser Menzies, Nigel Trimble, Roland Harwood, Rowland R. Kao,
# Adrian Allen, Theo Pepler, Robin Skuce
#
# Code developers: Assel Akhmetova, Theo Pepler, Liliana C. M. Salvador.
#
# Description: For each genotype combination and each KDE % level:
#   1. For each genotype in the combo, build a KDE from that genotype's
#      combined badger+cattle points
#   2. Keep only points (of either species) falling inside that genotype's
#      X% contour.
#   3. Recombine genotypes' filtered points into badger/cattle subsets.
#   4. if any genotype has <=2 samples in either species after contour
#     filtering, that combo x % cell is skipped
#
# Note: working directory setup was deleted from original version
#
# Input:  Supplementary_Table_2.xls / see Data Availability Statement]
# Output: Figure_7.pdf
# Requirements: R 4.4.2; packages: ks, readxl
# Contact: lilianasalvador@arizona.edu
# Last updated: October 1, 2026KDA across multiple KDE % contour thresholds, per shared genotype combo
##############################################################################

install.packages("ks")
install.packages("readxl")

library(ks)
library(readxl)

set.seed(12345)

TVRdata <- read_excel("Supplementary_Table_2.xlsx", sheet = 1,
                      col_names = TRUE, na = c("", "NA"))

# Match read.csv output: plain data.frame with read.csv-style column names
TVRdata <- as.data.frame(TVRdata)
names(TVRdata) <- make.names(names(TVRdata), unique = TRUE)

TVRdata$DVO <- ifelse(TVRdata$DVO == "Badger", 'B', 'C')
TVRdata$species <- TVRdata$DVO
TVRdata <- subset(TVRdata, !is.na(Genotype))

## KDE percentage levels to test
kde_pcts <- c(25, 50, 75, 95)

## Genotype combinations (same as analysis_comb.R)
genotype_combos <- list(
  c("4", "297"),
  c("4", "6"),
  c("6", "297"),
  c("4", "6", "297")
)

## Minimum samples per genotype per species required to run KDA safely.
min_n_per_genotype <- 3   

##############################################################################
## Contour membership: which points of a genotype's combined badger+cattle
## locations fall inside its own X% KDE contour
##############################################################################
get_contour_membership <- function(xy, pct) {
  # xy: data.frame/matrix with columns x, y for one genotype's combined
  #     badger+cattle points
  # pct: contour percentage (e.g. 25, 50, 80, 95)
  # returns: logical vector, TRUE = point falls inside the pct% contour

  if (nrow(xy) < 5) {
    # KDE is unstable/undefined with too few points; mark all as outside
    # so the downstream sample-size check flags and skips this cell
    return(rep(FALSE, nrow(xy)))
  }

  kde_fit <- ks::kde(x = as.matrix(xy))
  dens_at_points <- ks::kde(x = as.matrix(xy), eval.points = as.matrix(xy))$estimate
  thresh <- ks::contourLevels(kde_fit, cont = pct)
  dens_at_points >= thresh
}

##############################################################################
## Subset badger/cattle data to within-contour points, per genotype, for
## one combo and one KDE %
##############################################################################
filter_by_kde_pct <- function(ub, uc, geno_combo, pct) {
  ub_keep <- ub[0, ]
  uc_keep <- uc[0, ]

  for (g in geno_combo) {
    gb <- ub[ub$Genotype == g, ]
    gc <- uc[uc$Genotype == g, ]

    combined <- rbind(
      data.frame(x = gb$X.coord, y = gb$y.coord, species = "B"),
      data.frame(x = gc$X.coord, y = gc$y.coord, species = "C")
    )

    if (nrow(combined) < 5) next  # leaves gb/gc fully excluded; safety
                                    # check downstream will catch this

    inside <- get_contour_membership(combined[, c("x", "y")], pct)

    # combined was built as B-rows-then-C-rows in the same order as
    # gb/gc, so the species-split masks align directly with gb/gc rows
    gb_inside <- gb[inside[combined$species == "B"], , drop = FALSE]
    gc_inside <- gc[inside[combined$species == "C"], , drop = FALSE]

    ub_keep <- rbind(ub_keep, gb_inside)
    uc_keep <- rbind(uc_keep, gc_inside)
  }

  list(ub = ub_keep, uc = uc_keep)
}

##############################################################################
## Sample-size safety check: returns TRUE (safe to run) or FALSE (skip),
## plus a message explaining which genotype/species failed
##############################################################################
check_sample_sizes <- function(ub, uc, geno_combo) {
  problems <- character(0)
  for (g in geno_combo) {
    n_b <- sum(ub$Genotype == g)
    n_c <- sum(uc$Genotype == g)
    if (n_b < min_n_per_genotype) {
      problems <- c(problems, sprintf("badger n=%d for genotype %s", n_b, g))
    }
    if (n_c < min_n_per_genotype) {
      problems <- c(problems, sprintf("cattle n=%d for genotype %s", n_c, g))
    }
  }
  list(ok = length(problems) == 0, message = paste(problems, collapse = "; "))
}

##############################################################################
## Core KDA computation 
##############################################################################
calc_spatial_genotype <- function(badgerdata, cattledata) {
  badgerdata$X.coord <- as.numeric(badgerdata$X.coord)
  badgerdata$y.coord <- as.numeric(badgerdata$y.coord)
  cattledata$X.coord <- as.numeric(cattledata$X.coord)
  cattledata$y.coord <- as.numeric(cattledata$y.coord)

  alldata <- data.frame(rbind(badgerdata, cattledata))

  badger.EN <- subset(alldata, species == 'B', select = c("X.coord", "y.coord"))
  badger.EN <- as.matrix(badger.EN)
  badger.genotype <- as.numeric(alldata$Genotype[alldata$species == 'B'])

  cattle.EN <- subset(alldata, species == 'C', select = c("X.coord", "y.coord"))
  cattle.EN <- as.matrix(cattle.EN)
  cattle.genotype <- as.numeric(alldata$Genotype[alldata$species == 'C'])

  kda.out1 <- ks::kda(x = badger.EN, x.group = badger.genotype, eval.points = cattle.EN)
  pred.genotype1 <- kda.out1$x.group.estimate
  badger.cattle.rate <- sum(cattle.genotype != pred.genotype1) / length(cattle.genotype)

  kda.out2 <- ks::kda(x = cattle.EN, x.group = cattle.genotype, eval.points = badger.EN)
  pred.genotype2 <- kda.out2$x.group.estimate
  cattle.badger.rate <- sum(badger.genotype != pred.genotype2) / length(badger.genotype)

  list(badger.cattle.rate = badger.cattle.rate, cattle.badger.rate = cattle.badger.rate)
}


##############################################################################
## One-sided (left-tailed) permutation p-value 
##############################################################################
calc_pvalue <- function(stat, nulldist) {
  (sum(nulldist <= stat) + 1) / (length(nulldist) + 1)
}


##############################################################################
## Run one genotype combo at one KDE %, with the safety check gating
## whether the permutation analysis runs at all
##############################################################################
run_combo_at_pct <- function(data, geno_combo, pct, reps = 1000, seed = 12345) {

  combo_label <- paste(geno_combo, collapse = "_")

  d <- subset(data, Genotype %in% geno_combo)

  ub <- subset(d, DVO == 'B' & !duplicated(Tno),
               select = c("species", "Genotype", "X.coord", "y.coord"))
  uc <- subset(d, DVO == 'C' & !duplicated(paste(Genotype, Herd.ID)),
               select = c("species", "Genotype", "X.coord", "y.coord"))

  ub$X.coord <- as.numeric(ub$X.coord); ub$y.coord <- as.numeric(ub$y.coord)
  uc$X.coord <- as.numeric(uc$X.coord); uc$y.coord <- as.numeric(uc$y.coord)

  ## Contour membership computed ONCE on observed (true-label) data
  filtered <- filter_by_kde_pct(ub, uc, geno_combo, pct)
  ub_f <- filtered$ub
  uc_f <- filtered$uc

  ## Sample-size safety check on the CONTOUR-FILTERED data
  check <- check_sample_sizes(ub_f, uc_f, geno_combo)

  n_b <- table(factor(ub_f$Genotype, levels = geno_combo))
  n_c <- table(factor(uc_f$Genotype, levels = geno_combo))

  base_row <- data.frame(
    combination           = paste(geno_combo, collapse = ","),
    kde_pct                = pct,
    n_badgers_total        = nrow(ub_f),
    n_cattle_total         = nrow(uc_f),
    n_badgers_by_genotype  = paste(names(n_b), "=", as.integer(n_b), collapse = "; "),
    n_cattle_by_genotype   = paste(names(n_c), "=", as.integer(n_c), collapse = "; "),
    stringsAsFactors = FALSE
  )

  if (!check$ok) {
    message(sprintf("SKIPPED combo=%s pct=%d%%: %s", combo_label, pct, check$message))
    skipped_row <- cbind(base_row, reason = check$message, stringsAsFactors = FALSE)
    return(list(summary = NULL, skipped = skipped_row, panel_data = NULL))
  }

  ## Jitter coordinates to avoid singular matrices 
  set.seed(seed)
  ub_f$X.coord <- jitter(ub_f$X.coord, factor = 5)
  ub_f$y.coord <- jitter(ub_f$y.coord, factor = 5)
  uc_f$X.coord <- jitter(uc_f$X.coord, factor = 5)
  uc_f$y.coord <- jitter(uc_f$y.coord, factor = 5)

  obs <- calc_spatial_genotype(ub_f, uc_f)

  ## Permutation: reshuffle genotype labels WITHIN the fixed, already
  ## contour-filtered point sets (see note at top of script)
  perm_bc <- numeric(reps)
  perm_cb <- numeric(reps)
  for (i in seq_len(reps)) {
    pb <- ub_f; pc <- uc_f
    pb$Genotype <- sample(ub_f$Genotype)
    pc$Genotype <- sample(uc_f$Genotype)
    out <- calc_spatial_genotype(pb, pc)
    perm_bc[i] <- out$badger.cattle.rate
    perm_cb[i] <- out$cattle.badger.rate
  }

  p_bc <- calc_pvalue(obs$badger.cattle.rate, perm_bc)
  p_cb <- calc_pvalue(obs$cattle.badger.rate, perm_cb)

  summary_row <- cbind(base_row, data.frame(
    rate_badger_to_cattle   = obs$badger.cattle.rate,
    pvalue_badger_to_cattle = p_bc,
    rate_cattle_to_badger   = obs$cattle.badger.rate,
    pvalue_cattle_to_badger = p_cb,
    stringsAsFactors = FALSE
  ))

  list(
    summary = summary_row,
    skipped = NULL,
    panel_data = list(
      geno_combo = geno_combo, pct = pct,
      perm_bc = perm_bc, perm_cb = perm_cb,
      obs_bc = obs$badger.cattle.rate, obs_cb = obs$cattle.badger.rate,
      p_bc = p_bc, p_cb = p_cb
    )
  )
}


##############################################################################
## Run everything: all combos x all KDE %
##############################################################################
all_runs <- list()
idx <- 1
for (pct in kde_pcts) {
  for (g in genotype_combos) {
    message(sprintf("Running combination %s at KDE %d%%", paste(g, collapse=","), pct))
    all_runs[[idx]] <- run_combo_at_pct(TVRdata, g, pct, reps = 10)
    idx <- idx + 1
  }
}

saveRDS(all_runs, file.path(project_dir, "all_runs.rds"))     
message("Saved all_runs.rds -- reload this to redo figures without rerunning the analysis")

## Collect results and skipped cells separately
results <- do.call(rbind, lapply(all_runs, function(x) x$summary))
skipped <- do.call(rbind, lapply(all_runs, function(x) x$skipped))

print(results)
write.csv(results, "spatial_genotype_results_by_kde_pct.csv", row.names = FALSE)

if (!is.null(skipped)) {
  cat("\n--- SKIPPED combinations (insufficient sample size after contour filtering) ---\n")
  print(skipped)
  write.csv(skipped, "spatial_genotype_SKIPPED_by_kde_pct.csv", row.names = FALSE)
} else {
  cat("\nNo combinations were skipped.\n")
}


########################################################################
########################################################################
# Code to generate figure 7
# This can be run separately if all_runs.rds has already been generated

# Load the all_runs.rds file
all_runs <- readRDS("all_runs.rds")

message("Loaded all_runs.rds")

# Define KDE percentages and genotype combos
kde_pcts <- c(25, 50, 75, 95)

genotype_combos <- list(
  c("4", "297"),
  c("4", "6"),
  c("6", "297"),
  c("4", "6", "297")
)

# Color assignments
combo_colors <- c()
combo_colors["4,297"]    <- "#1F77B4"
combo_colors["4,6"]      <- "#FFB81C"
combo_colors["6,297"]    <- "#17BECF"
combo_colors["4,6,297"]  <- "#D7573B"

################################################################################
# Helper function: get_panel_data
################################################################################
get_panel_data <- function(geno_combo, kde_pct, direction) {
  matching_run <- NULL
  
  for (run_idx in seq_along(all_runs)) {
    run <- all_runs[[run_idx]]
    
    if (is.null(run$panel_data)) next
    if (run$panel_data$pct != kde_pct) next
    if (!identical(run$panel_data$geno_combo, geno_combo)) next
    
    matching_run <- run
    break
  }
  
  if (is.null(matching_run)) return(NULL)
  
  pd <- matching_run$panel_data
  combo_str <- paste(pd$geno_combo, collapse = ",")
  
  if (direction == "bc") {
    return(list(
      combo_str = combo_str,
      perm_rates = pd$perm_bc,
      obs_rate = pd$obs_bc,
      pval = pd$p_bc
    ))
  }
  
  if (direction == "cb") {
    return(list(
      combo_str = combo_str,
      perm_rates = pd$perm_cb,
      obs_rate = pd$obs_cb,
      pval = pd$p_cb
    ))
  }
  
  return(NULL)
}

################################################################################
# Function to create ONE panel
################################################################################
create_panel <- function(combo_data_list, kde_pct, y_max, break_points = NULL, show_title = TRUE) {
  valid_combos <- c()
  for (i in seq_along(combo_data_list)) {
    if (!is.null(combo_data_list[[i]])) {
      valid_combos <- c(valid_combos, i)
    }
  }
  
  if (length(valid_combos) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data", size = 3) +
             theme_void())
  }
  
  plot_data <- data.frame()
  obs_data <- data.frame()
  
  for (idx_num in seq_along(valid_combos)) {
    idx <- valid_combos[idx_num]
    combo <- combo_data_list[[idx]]
    
    dens <- density(combo$perm_rates, adjust = 2)
    
    combo_dens_df <- data.frame(
      x = dens$x,
      y = dens$y,
      combo = combo$combo_str,
      pval = combo$pval,
      obs_rate = combo$obs_rate
    )
    
    plot_data <- rbind(plot_data, combo_dens_df)
    
    obs_row <- data.frame(
      obs_rate = combo$obs_rate,
      combo = combo$combo_str,
      pval = combo$pval
    )
    
    obs_data <- rbind(obs_data, obs_row)
  }
  
  title_text <- if (show_title) sprintf("KDE %d%%", kde_pct) else ""
  
  p <- ggplot(plot_data, aes(x = x, y = y, color = combo, fill = combo)) +
    geom_ribbon(aes(ymin = 0, ymax = y, fill = combo), 
                alpha = 0.3, color = NA) +
    geom_line(linewidth = 0.6) +
    geom_vline(data = obs_data, aes(xintercept = obs_rate, color = combo), 
               linetype = 2, linewidth = 0.6, show.legend = FALSE) +
    scale_color_manual(values = combo_colors) +
    scale_fill_manual(values = combo_colors) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
    labs(x = "Misclassification rate",
         y = "Density",
         title = title_text,
         color = "Genotype",
         fill = "Genotype") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.3, "cm")
    )
  
  # Add axis breaks if provided
  if (!is.null(break_points)) {
    for (bp in break_points) {
      p <- p + scale_y_break(breaks = bp, scales = 0.3, space = 0.2)
    }
    # Remove the duplicated right y-axis that ggbreak adds
    p <- p + theme(
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      axis.title.y.right = element_blank()
    )
  }
  
  return(p)
}

################################################################################
# Break points for Cattle->Badger panels
################################################################################
breaks_per_cb_panel <- list(
  list(c(100, 280)),                    # 25%
  list(c(100, 280)),                    # 50%
  list(c(100, 280)),                    # 75%
  list(c(50, 170), c(220, 280))         # 95%
)

################################################################################
# Generate all 8 panels as single-page PDFs
################################################################################
bc_files <- c()
cb_files <- c()

# Badger -> Cattle panels (top row)
for (i in 1:4) {
  pct <- kde_pcts[i]
  combo_list <- list()
  
  for (j in seq_along(genotype_combos)) {
    geno <- genotype_combos[[j]]
    data <- get_panel_data(geno, pct, "bc")
    combo_list[[j]] <- data
  }
  
  p <- create_panel(combo_list, pct, y_max = 30, break_points = NULL, show_title = TRUE)
  
  filename <- sprintf("temp_bc_kde_%d.pdf", pct)
  pdf(filename, width = 4, height = 4, onefile = FALSE)
  print(p)
  dev.off()
  
  bc_files <- c(bc_files, filename)
  message(sprintf("Saved: %s", filename))
}

# Cattle -> Badger panels (bottom row)
for (i in 1:4) {
  pct <- kde_pcts[i]
  break_points <- breaks_per_cb_panel[[i]]
  combo_list <- list()
  
  for (j in seq_along(genotype_combos)) {
    geno <- genotype_combos[[j]]
    data <- get_panel_data(geno, pct, "cb")
    combo_list[[j]] <- data
  }
  
  p <- create_panel(combo_list, pct, y_max = 320, break_points = break_points, show_title = FALSE)
  
 filename <- sprintf("temp_cb_kde_%d.pdf", pct)
 pdf(filename, width = 4, height = 4, onefile = FALSE)
 print(p)
 dev.off()
  
  cb_files <- c(cb_files, filename)
  message(sprintf("Saved: %s", filename))
}

################################################################################
# Combine all 8 panels into a 2x4 grid using magick
################################################################################

# Read all PDFs as images
bc_images <- lapply(bc_files, function(f) image_read_pdf(f, density = 300))
cb_images <- lapply(cb_files, function(f) image_read_pdf(f, density = 300))

# Combine top row (Badger -> Cattle)
top_row <- image_append(do.call(c, bc_images))

# Combine bottom row (Cattle -> Badger)
bottom_row <- image_append(do.call(c, cb_images))

# Get dimensions
top_info <- image_info(top_row)
bottom_info <- image_info(bottom_row)

# Create row label images (blank with rotated text) - to the LEFT of the plots
label_width <- 120  # pixels of extra space for row labels

# Top row label: Badger -> Cattle
top_label <- image_blank(width = label_width, height = top_info$height, color = "white")
top_label <- image_annotate(top_label, "Badger -> Cattle",
                            size = 40,
                            gravity = "center",
                            degrees = 270,
                            weight = 700,
                            color = "black")

# Bottom row label: Cattle -> Badger
bottom_label <- image_blank(width = label_width, height = bottom_info$height, color = "white")
bottom_label <- image_annotate(bottom_label, "Cattle -> Badger",
                               size = 40,
                               gravity = "center",
                               degrees = 270,
                               weight = 700,
                               color = "black")

# Combine label + row for both directions
top_row_labeled <- image_append(c(top_label, top_row))
bottom_row_labeled <- image_append(c(bottom_label, bottom_row))

# Stack rows vertically
combined <- image_append(c(top_row_labeled, bottom_row_labeled), stack = TRUE)

# Save final combined figure
image_write(combined, "Figure_7.pdf", format = "pdf")

message("Done! Combined figure saved: Figure_7.pdf")

################################################################################
# Optional: Clean up temporary files
################################################################################
# Uncomment to remove temporary panel PDFs after combining:
# file.remove(c(bc_files, cb_files))