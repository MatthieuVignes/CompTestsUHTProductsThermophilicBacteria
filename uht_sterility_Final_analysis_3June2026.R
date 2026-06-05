##################################################
# UHT Sterility Testing - all 5 Products
# Matthieu Vignes, for Julie Warren (Fonterra) - 3 June 2026
##################################################
# Products : UHT Milk, Cream, Whipping Cream, Beverage, Chocolate Milk
# Tests    : Plate Count (reference), Charm (RLU/test), Attune Diluted (AFU/mL)
# Bacteria : 8 common + 4 Chocolate-specific (BL 1, BL 2, AF 5 Hi, AF 8 Hi)
# Design   : Blank + 5 replicates (A–E) × 8/12 bacteria × 3–5 time-points
##################################################
# OUTPUT FILES (written to ./output/)
#   Plots (PDF)  : one per product/test combination
#   Tables (CSV) : sweep results per product × bacteria × time × threshold
#   Summary CSV  : best thresholds per product × test × time-point
##################################################
library(tidyverse)   # data wrangling + ggplot2
library(patchwork)   # multi-panel figures
library(scales)      # axis formatting
library(kableExtra)  # (optional) pretty tables if knitting later

setwd("~/Documents/Work-notMassey/FONTERRA/Sterility_Testing_Julie/") # enter yours!
dir.create("output", showWarnings = FALSE)

##################################################
# SECTION 1 - CONSTANTS & HELPERS
##################################################
PLATE_POS   <- 10    # CFU/mL: plate count > 10 is considered positive
PLATE_CENS  <- 300   # values above this are TNTC; plotted as ">300", stored as 305
PLATE_STORE <- 305   # numeric stand-in for >300 in all calculations

CHARM_MFR   <- 300   # manufacturer recommended Charm threshold (RLU/test)
BG_N_SD     <- 2     # blank band = mean +/- BG_N_SD * sd (log10 scale)

# Canonical time-point order (covers all products; subset used per product)
TIME_ORDER <- c("6 Hours","12 Hours","24 Hours","48 Hours","72 Hours")

# Colour palette for bacteria (up to 12)
BACT_COLS <- setNames(
  c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
    "#FF7F00","#A65628","#F781BF","#999999",
    "#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),
  c("Geo 1","Geo 2","Geo 3","Geo 4",
    "AF 5","AF 6","AF 7","AF 8",
    "BL 1","BL 2","AF 5 Hi","AF 8 Hi")
)

# Parse plate count: strip ">", coerce to numeric, cap at PLATE_STORE
parse_plate <- function(x) {
  v <- as.numeric(gsub(">\\s*", "", as.character(x)))
  ifelse(v > PLATE_CENS, PLATE_STORE, v)
}

# Standardise bacteria names (trim whitespace, unify "AF 8 Hi" variants)
clean_bact <- function(x) {
  x <- trimws(x)
  x <- gsub("AF\\s+8\\s+Hi.*",  "AF 8 Hi",  x)
  x <- gsub("AF\\s+5\\s+Hi.*",  "AF 5 Hi",  x)
  x
}

# Standardise time labels -> numeric hours value stored in TimeNum,
# plus an ordered factor Time used for colour/legend ordering.
# Numeric hours allow proportional spacing on the x-axis.
clean_time <- function(x) {
  x <- trimws(x)
  x <- sub("(?i)^(\\d+)\\s+hours?$", "\\1 Hours", x, perl = TRUE)
  factor(x, levels = TIME_ORDER)
}

# Extract the numeric part of "X Hours" -> integer (6, 12, 24, 48, 72)
time_to_num <- function(f) {
  as.integer(sub(" Hours", "", as.character(f)))
}

# Safe log10: zeros -> NA (handled explicitly later)
log10s <- function(x) ifelse(x > 0, log10(x), NA_real_)

# Performance metrics from a 2×2 confusion table
# Returns a named numeric vector; NA where denominator = 0
metrics <- function(TP, FP, FN, TN) {
  N   <- TP + FP + FN + TN
  Sen <- if ((TP + FN) > 0) TP / (TP + FN)    else NA_real_
  Spe <- if ((TN + FP) > 0) TN / (TN + FP)    else NA_real_
  PPV <- if ((TP + FP) > 0) TP / (TP + FP)    else NA_real_
  NPV <- if ((TN + FN) > 0) TN / (TN + FN)    else NA_real_
  Prev <- (TP + FN) / N
  Acc  <- (TP + TN) / N
  F1   <- if (!is.na(Sen) & !is.na(PPV) & (Sen + PPV) > 0)
    2 * Sen * PPV / (Sen + PPV) else NA_real_
  Jacc <- if ((TP + FP + FN) > 0) TP / (TP + FP + FN) else NA_real_
  denom_MCC <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC  <- if (denom_MCC > 0) (TP*TN - FP*FN) / denom_MCC else NA_real_
  YouJ <- if (!is.na(Sen) & !is.na(Spe)) Sen + Spe - 1 else NA_real_
  c(TP=TP, FP=FP, FN=FN, TN=TN,
    Sensitivity=Sen, Specificity=Spe, PPV=PPV, NPV=NPV,
    Prevalence=Prev, Accuracy=Acc, F1=F1, Jaccard=Jacc, MCC=MCC,
    Youden_J=YouJ)
}

##################################################
# SECTION 2 - DATA LOADING
##################################################
# File paths and official product names
PRODUCTS <- list(
  Milk          = list(path = "FinalMilk.csv",          name = "UHT Milk"),
  Cream         = list(path = "FinalCream.csv",          name = "Cream"),
  WhippingCream = list(path = "FinalWhippingCream.csv",  name = "Whipping Cream"),
  MedBeverage   = list(path = "FinalMedBeverage.csv",    name = "Med Beverage"),
  Chocolate     = list(path = "FinalChocolate.csv",      name = "Chocolate")
)

read_product <- function(path, prod_name) {
  raw <- read.csv(path, skip = 2, header = FALSE, stringsAsFactors = FALSE)
  colnames(raw) <- c("Bact","Reps","Time","Plate","Charm","AttuneRaw","AttuneDilut")
  raw |>
    mutate(
      Product    = prod_name,
      Bact       = factor(clean_bact(Bact), levels = names(BACT_COLS)),
      Reps       = trimws(Reps),
      Time        = clean_time(Time),
      TimeNum     = time_to_num(clean_time(Time)),  # numeric hours for x-axis spacing
      Plate       = parse_plate(Plate),
      Charm       = as.numeric(Charm),
      AttuneDilut = as.numeric(AttuneDilut),
      # Positivity flag (for analysis): >10 CFU/mL = positive
      PlatePos    = as.integer(Plate > PLATE_POS),
      # For zeros in Charm/Attune: replace with 0.5 so log10 is finite
      # and the value is always below any operationally meaningful threshold
      Charm_adj       = if_else(Charm       == 0 | is.na(Charm),       0.5, Charm),
      AttuneDilut_adj = if_else(AttuneDilut == 0 | is.na(AttuneDilut), 0.5, AttuneDilut)
    ) |>
    # Drop Attune Raw per instructions; also drop any rows where Bact is NA
    # (catches stray header/comment rows, e.g. Chocolate row 641 in inoculated)
    filter(!is.na(Bact)) |>
    select(-AttuneRaw)
}

all_data <- imap_dfr(PRODUCTS, ~ read_product(.x$path, .x$name)) |>
  mutate(
    is_blank    = (Reps == "Blank"),
    Product     = factor(Product, levels = map_chr(PRODUCTS, "name"))
  )

blanks     <- filter(all_data,  is_blank)
inoculated <- filter(all_data, !is_blank)

# Pooled blank background per product (all bacteria combined).
# Stats are on log10 scale to match the signal plots.
# Used to draw a shaded noise-floor band on Charm and Attune signal plots.
blank_bg <- blanks |>
  group_by(Product) |>
  summarise(
    charm_bg_mean  = mean(log10s(Charm_adj),       na.rm = TRUE),
    charm_bg_sd    = sd(log10s(Charm_adj),         na.rm = TRUE),
    attune_bg_mean = mean(log10s(AttuneDilut_adj), na.rm = TRUE),
    attune_bg_sd   = sd(log10s(AttuneDilut_adj),   na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    charm_bg_lo  = charm_bg_mean  - BG_N_SD * charm_bg_sd,
    charm_bg_hi  = charm_bg_mean  + BG_N_SD * charm_bg_sd,
    attune_bg_lo = attune_bg_mean - BG_N_SD * attune_bg_sd,
    attune_bg_hi = attune_bg_mean + BG_N_SD * attune_bg_sd
  )

cat("Data loaded:\n")
inoculated |> count(Product, Time)

##################################################
# SECTION 3 - PLOT FUNCTIONS
##################################################
# Colour map for time-points
time_pal <- c("6 Hours"  = "#1B9E77",
              "12 Hours" = "#D95F02",
              "24 Hours" = "#7570B3",
              "48 Hours" = "#E7298A",
              "72 Hours" = "#66A61E")

# 3A: Plate Count plot
plot_plate <- function(df, prod_name) {
  time_breaks <- sort(unique(df$TimeNum))
  df |>
    mutate(PlotPlate = pmin(Plate, PLATE_STORE)) |>
    ggplot(aes(x = TimeNum, y = PlotPlate, colour = Bact,
               group = interaction(Bact, Reps))) +
    geom_line(alpha = 0.55, linewidth = 0.5) +
    geom_point(size = 1.5, alpha = 0.75,
               position = position_jitter(width = 0.3, seed = 42)) +
    geom_hline(yintercept = PLATE_POS, linetype = "dashed",
               colour = "firebrick", linewidth = 0.7) +
    annotate("text", x = max(time_breaks), y = PLATE_POS,
             label = paste0(">", PLATE_POS, " CFU/mL = positive"),
             hjust = 1, vjust = -0.4, colour = "firebrick", size = 2.8) +
    scale_x_continuous(breaks = time_breaks,
                       labels = paste0(time_breaks, "h")) +
    scale_y_continuous(
      breaks = c(0, 50, 100, 150, 200, 250, 300, PLATE_STORE),
      labels = c("0","50","100","150","200","250","300",">300")
    ) +
    scale_colour_manual(values = BACT_COLS, drop = TRUE) +
    labs(title = paste(prod_name, "- Plate Count"),
         x = "Time (h)", y = "Plate Count (CFU/mL)", colour = "Bacteria") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
}

# 3B: Generic log10 signal plot (Charm or Attune)
# bg is a named list with elements $lo and $hi (both on log10 scale),
# derived from pooled blank stats: mean +/- BG_N_SD * sd.
# A grey shaded band is drawn first so threshold lines sit on top.
plot_signal <- function(df, test_col, ylabel, hlines = NULL, hline_labels = NULL,
                        prod_name = "", bg = NULL) {
  time_breaks <- sort(unique(df$TimeNum))
  p <- df |>
    mutate(log_val = log10s(.data[[test_col]])) |>
    filter(!is.na(log_val)) |>
    ggplot(aes(x = TimeNum, y = log_val, colour = Bact,
               group = interaction(Bact, Reps))) +
    geom_line(alpha = 0.55, linewidth = 0.5) +
    geom_point(size = 1.5, alpha = 0.75) +
    scale_x_continuous(breaks = time_breaks,
                       labels = paste0(time_breaks, "h")) +
    scale_colour_manual(values = BACT_COLS, drop = TRUE) +
    labs(title = paste(prod_name, "-", ylabel),
         x = "Time (h)", y = paste0("log10(", ylabel, ")"), colour = "Bacteria") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
  # Blank noise-floor band: drawn before threshold lines so lines sit on top.
  # Band spans full x range; label placed at left edge just above the band.
  if (!is.null(bg) && !is.na(bg$lo) && !is.na(bg$hi)) {
    p <- p +
      annotate("rect",
               xmin = -Inf, xmax = Inf,
               ymin = bg$lo, ymax = bg$hi,
               fill = "grey70", alpha = 0.35) +
      annotate("text",
               x = min(time_breaks), y = bg$hi,
               label = paste0("blank mean ± ", BG_N_SD, " SD"),
               hjust = 0, vjust = -0.4, colour = "grey40", size = 2.4)
  }
  if (!is.null(hlines)) {
    for (i in seq_along(hlines)) {
      p <- p +
        geom_hline(yintercept = log10(hlines[i]), linetype = "dashed",
                   colour = "firebrick", linewidth = 0.6) +
        annotate("text",
                 x     = max(time_breaks),
                 y     = log10(hlines[i]),
                 label = hline_labels[i],
                 hjust = 1, vjust = -0.4, colour = "firebrick", size = 2.5)
    }
  }
  p
}

# 3C: Per-bacteria boxplot vs Plate Count
plot_vs_plate <- function(df, test_col, ylabel, prod_name) {
  df |>
    mutate(
      log_val  = log10s(.data[[test_col]]),
      # Use a fixed two-level factor so labels are always valid even when
      # only one plate-count class is present (e.g. Chocolate, all positive).
      PlateGrp = factor(PlatePos,
                        levels = c(0L, 1L),
                        labels = c("Plate >= 10", "Plate >10"))
    ) |>
    filter(!is.na(log_val)) |>
    ggplot(aes(x = TimeNum, y = log_val, fill = PlateGrp, colour = PlateGrp,
               group = interaction(TimeNum, PlateGrp))) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.4, size = 1, alpha = 0.65) +
    facet_wrap(~ Bact, nrow = 2) +
    scale_x_continuous(
      breaks = sort(unique(df$TimeNum)),
      labels = paste0(sort(unique(df$TimeNum)), "h")
    ) +
    scale_fill_manual(values  = c("steelblue","tomato")) +
    scale_colour_manual(values = c("steelblue","tomato")) +
    labs(title  = paste(prod_name, "-", ylabel, "vs Plate Count"),
         x = "Time (h)", y = paste0("log10(", ylabel, ")"),
         fill = "Plate count", colour = "Plate count") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7),
          legend.position = "top")
}

# 3D: Accuracy vs threshold (one curve per bacteria)
plot_accuracy_curves <- function(sweep_df, time_val, prod_name, test_label) {
  sweep_df |>
    filter(Time == time_val, !is.na(Accuracy)) |>
    ggplot(aes(x = log10_thresh, y = Accuracy, colour = Bact)) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    scale_colour_manual(values = BACT_COLS, drop = TRUE) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(1)) +
    labs(title   = paste(prod_name, test_label, time_val),
         x       = paste0("Threshold (log10)"),
         y       = "Accuracy",
         colour  = "Bacteria") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
}

# 3E: Mean Sensitivity & Specificity vs threshold
# For time-points where all samples are plate-positive, Specificity is
# NA for every threshold and every bacterium. mean(NA, na.rm=TRUE) = NaN,
# which ggplot silently plots as 0 - a misleading flat line. We keep the
# Sensitivity line for such time-points (it is valid) but suppress the
# Specificity line and flag the time-point in the subtitle.
plot_mean_sensspe <- function(sweep_df, prod_name, test_label) {
  # Identify degenerate time-points (no valid Specificity anywhere)
  degen_times <- sweep_df |>
    group_by(Time) |>
    summarise(all_spe_na = all(is.na(Specificity)), .groups = "drop") |>
    filter(all_spe_na) |>
    pull(Time) |>
    as.character()
  
  ss_data <- sweep_df |>
    group_by(Time, log10_thresh) |>
    summarise(
      Sen = mean(Sensitivity, na.rm = TRUE),
      # Set Spe to NA for degenerate time-points so the line is suppressed
      Spe = if_else(as.character(Time[1]) %in% degen_times,
                    NA_real_,
                    mean(Specificity, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    pivot_longer(c(Sen, Spe), names_to = "Metric", values_to = "Value") |>
    mutate(Metric = recode(Metric, Sen = "Sensitivity", Spe = "Specificity")) |>
    filter(!is.na(Value))
  
  subtitle_txt <- if (length(degen_times) > 0)
    paste0("Specificity not shown for: ", paste(degen_times, collapse = ", "),
           " (all samples plate-positive)")
  else ""
  
  ggplot(ss_data, aes(x = log10_thresh, y = Value,
                      colour = Time, linetype = Metric)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_colour_manual(values = time_pal, drop = TRUE) +
    scale_linetype_manual(values = c(Sensitivity = "solid", Specificity = "dashed")) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(1)) +
    labs(title    = paste(prod_name, "-", test_label,
                          "mean Sensitivity & Specificity"),
         subtitle = subtitle_txt,
         x        = "Threshold (log10)",
         y        = NULL,
         colour   = "Time", linetype = NULL) +
    theme_bw(base_size = 10)
}

# 3F: Mean ROC (one curve per time-point)
# A degenerate time-point (all samples plate-positive → no true negatives)
# makes Specificity undefined for every threshold. mean(NA, na.rm=TRUE)
# returns NaN, which ggplot maps to FPR=0, producing a spurious vertical
# line from (0,0) to (0,1). We detect and exclude such time-points and
# note them in the subtitle so the reader knows they are missing.
plot_mean_roc <- function(sweep_df, prod_name, test_label) {
  roc_data <- sweep_df |>
    group_by(Time, log10_thresh) |>
    summarise(
      # Count how many bacteria had a valid (non-NA) Specificity value
      n_valid_spe = sum(!is.na(Specificity)),
      FPR = mean(1 - Specificity, na.rm = TRUE),
      TPR = mean(Sensitivity,     na.rm = TRUE),
      .groups = "drop"
    ) |>
    # Exclude time-points where NO bacterium had a computable Specificity
    # (n_valid_spe == 0 means every bacterium had only positives → FPR is
    # meaningless and would produce a degenerate vertical line at FPR=0)
    filter(n_valid_spe > 0, !is.na(TPR)) |>
    # Sort by DECREASING threshold to trace the ROC correctly from the
    # top-left corner (high threshold = few positives called) to the
    # bottom-right (low threshold = everything called positive).
    # Sorting by FPR instead causes ties at FPR=0 (many thresholds with
    # no false positives yet) to be connected in arbitrary TPR order,
    # producing a spurious vertical line on the left edge that appears
    # as a second incompatible curve on the plot.
    arrange(Time, desc(log10_thresh))
  
  # Identify excluded time-points for the subtitle
  all_times     <- levels(droplevels(sweep_df$Time))
  plotted_times <- as.character(unique(roc_data$Time))
  excluded      <- setdiff(all_times, plotted_times)
  subtitle_txt  <- if (length(excluded) > 0)
    paste0("Mean over all bacteria; one curve per time-point
",
           "Note: ", paste(excluded, collapse = ", "),
           " excluded (all samples plate-positive → Specificity undefined)")
  else
    "Mean over all bacteria; one curve per time-point"
  
  ggplot(roc_data, aes(x = FPR, y = TPR, colour = Time)) +
    geom_line(linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = time_pal, drop = TRUE) +
    scale_x_continuous(labels = percent_format(1), limits = c(0, 1)) +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
    labs(title    = paste(prod_name, "-", test_label, "mean ROC"),
         subtitle = subtitle_txt,
         x        = "1 - Specificity (FPR)",
         y        = "Sensitivity (TPR)",
         colour   = "Time") +
    theme_bw(base_size = 10)
}

##################################################
# SECTION 4 - THRESHOLD SWEEP ANALYSIS
##################################################
# Build a grid of thresholds for a given test column (log10 scale)
make_thres_grid <- function(df, col, step = 0.1) {
  vals <- df[[col]][df[[col]] > 0 & !is.na(df[[col]])]
  lo   <- floor(log10(min(vals)) * 10) / 10
  hi   <- ceiling(log10(max(vals)) * 10) / 10
  seq(lo, hi, by = step)
}

# Runs the full sweep for one test column on one product.
# blanks_df : blank rows for this product (PlatePos always 0).
#   For each Bact × Time cell, blank rows with a matching (Bact, Time) in
#   the inoculated data are appended before computing TP/FP/FN/TN.
#   This gives real FP counts (a blank above the threshold = false positive)
#   and makes Specificity meaningful even at time-points where all inoculated
#   replicates are plate-positive.
#   Blanks with no matching inoculated cell are silently ignored.
# The threshold grid is derived from inoculated data only so that very low
#   blank signals do not drag the grid floor down unnecessarily.
# Returns a tidy data frame with one row per (Bact × Time × threshold).
run_sweep <- function(df, test_col, blanks_df = NULL) {
  thres_grid <- make_thres_grid(df, test_col)
  
  map_dfr(levels(df$Bact)[levels(df$Bact) %in% unique(df$Bact)], function(b) {
    map_dfr(levels(df$Time)[levels(df$Time) %in% unique(df$Time)], function(t) {
      sub_inoc <- df[df$Bact == b & df$Time == t, ]
      if (nrow(sub_inoc) == 0) return(NULL)
      
      # Append matching blank rows (same Bact and Time) if available.
      # Blanks are always PlatePos = 0 (plate count = 0, well below 10).
      sub_blank <- if (!is.null(blanks_df))
        blanks_df[blanks_df$Bact == b & blanks_df$Time == t, ]
      else
        sub_inoc[0, ]   # empty frame with same columns
      
      sub <- bind_rows(sub_inoc, sub_blank)
      
      plate_pos <- sub$PlatePos            # 0/1
      log_test  <- log10s(sub[[test_col]]) # NA where original was 0
      
      map_dfr(thres_grid, function(thr) {
        test_pos <- !is.na(log_test) & (log_test > thr)
        # Rows where log_test is NA are excluded from all cells
        valid    <- !is.na(log_test)
        TP <- sum( plate_pos[valid] == 1 &  test_pos[valid])
        FP <- sum( plate_pos[valid] == 0 &  test_pos[valid])
        FN <- sum( plate_pos[valid] == 1 & !test_pos[valid])
        TN <- sum( plate_pos[valid] == 0 & !test_pos[valid])
        m  <- metrics(TP, FP, FN, TN)
        as_tibble_row(c(
          list(Bact = b, Time = t, log10_thresh = thr),
          as.list(m)
        ))
      })
    })
  }) |>
    mutate(
      Bact = factor(Bact, levels = names(BACT_COLS)),
      Time = factor(Time, levels = TIME_ORDER)
    )
}

##################################################
# SECTION 5 - TWO-STAGE THRESHOLD SELECTION
##################################################
# NOTE: not what is exactly used in the paper, see method description. Was used as a guide. 
# Youden, "best" Sensitivity=1 and above balk/bg threshold was explored
# Stage 1 - Hard background floor
#   Only thresholds strictly above bg$hi (mean + BG_N_SD*sd of blanks on
#   log10 scale) are considered. This eliminates thresholds that cannot
#   distinguish signal from instrument noise, regardless of classification
#   performance.
#
# Stage 2 - Three nested optima within the above-floor candidates
#   (a) Youden's J = Sensitivity + Specificity - 1  (balanced optimum,
#       standard in diagnostic test literature, Youden 1950)
#   (b) Sensitivity-first range: all thresholds achieving the highest
#       attainable Sensitivity (= 1 if possible), then among those the
#       range that also maximises Specificity
#   Both are computed per bacterium then intersected across all bacteria
#   at each time-point to give a common operating range.
#   NA/NaN values are excluded at every step and never used to restrict
#   the range.
##################################################
# intersection of per-bacteria ranges
# per_bact : data frame with (Time, Bact, lower, upper, achieved) - one row
#            per bacterium per time-point, produced by youden_optimum() or
#            sensitivity_first(). Contains only bacteria that survived the
#            Sensitivity filter (i.e. had at least one plate-positive replicate
#            above the background floor).
# bact_counts : named integer vector Time -> n_total_bact, giving the number
#            of bacteria present in the raw inoculated data for each time-point
#            (before any filtering). Used to report n_above_floor / n_total.
#
# Three bacteria counts are reported per time-point:
#   n_total        : all bacteria present in the raw data at this time-point
#   n_above_floor  : bacteria that had >= 1 threshold above the background
#                    floor AND had at least one plate-positive replicate
#                    (i.e. rows reaching intersect_ranges)
#   n_in_intersect : bacteria whose individual range overlaps the common
#                    intersection (n_above_floor minus problematic ones)
intersect_ranges <- function(per_bact, bact_counts) {
  per_bact |>
    group_by(Time) |>
    summarise(
      common_lower   = {
        vals <- as.numeric(lower[!is.na(lower)])
        if (length(vals)) max(vals) else NA_real_
      },
      common_upper   = {
        vals <- as.numeric(upper[!is.na(upper)])
        if (length(vals)) min(vals) else NA_real_
      },
      # Bacteria that reached this function (survived Sensitivity filter
      # and had at least one above-floor threshold)
      n_above_floor  = n(),
      # Bacteria with no valid individual range (lower is NA - optimisation
      # found no threshold meeting the criterion for this bacterium)
      no_range        = {
        b <- Bact[is.na(lower)]
        if (length(b)) paste(b, collapse = ", ") else NA_character_
      },
      # Bacteria whose individual range exists but does not overlap the
      # common intersection - they are the binding constraint
      problematic     = {
        lo_valid <- as.numeric(lower[!is.na(lower)])
        hi_valid <- as.numeric(upper[!is.na(upper)])
        if (!length(lo_valid) || !length(hi_valid)) {
          NA_character_
        } else {
          b <- Bact[!is.na(lower) & !is.na(upper) &
                      (as.numeric(lower) > min(hi_valid) |
                         as.numeric(upper) < max(lo_valid))]
          if (length(b)) paste(b, collapse = ", ") else NA_character_
        }
      },
      # Best achieved metric value for bacteria with no valid range
      best_achieved   = {
        idx <- is.na(lower)
        if (!any(idx)) NA_character_
        else paste(paste0(Bact[idx], "=", round(as.numeric(achieved[idx]), 2)),
                   collapse = ", ")
      },
      .groups = "drop"
    ) |>
    mutate(
      # Total bacteria in raw data for this time-point (from bact_counts)
      n_total        = as.integer(bact_counts[as.character(Time)]),
      # Bacteria in intersection = above-floor minus those with no range
      # or whose range doesn't overlap
      n_problematic  = lengths(strsplit(
        if_else(is.na(problematic), "", problematic), ",\\s*"
      )) * (!is.na(problematic)),
      n_no_range     = lengths(strsplit(
        if_else(is.na(no_range), "", no_range), ",\\s*"
      )) * (!is.na(no_range)),
      n_in_intersect = n_above_floor - n_no_range - n_problematic,
      valid          = !is.na(common_lower) & !is.na(common_upper) &
        common_lower <= common_upper,
      common_lower   = if_else(valid, common_lower, NA_real_),
      common_upper   = if_else(valid, common_upper, NA_real_)
    ) |>
    select(-n_problematic, -n_no_range)
}

# Stage 1: apply background floor to a sweep data frame
# bg_hi : upper edge of blank noise band (log10 scale) for this product/test.
# Returns the sweep filtered to thresholds > bg_hi only.
# Rows with NA log10_thresh are dropped; NaN is treated as NA.
apply_bg_floor <- function(sweep_df, bg_hi) {
  if (is.null(bg_hi) || is.na(bg_hi) || is.nan(bg_hi)) {
    warning("Background upper limit is NA/NaN - floor not applied.")
    return(sweep_df)
  }
  sweep_df |> filter(!is.na(log10_thresh), log10_thresh > bg_hi)
}

# Stage 2a: Youden's J optimum
# For each (Bact, Time): find the threshold maximising J = Sen + Spe - 1
# among above-floor candidates. Sen or Spe = NA at a threshold means that
# threshold is skipped (J cannot be computed). Returns per-bact and common.
youden_optimum <- function(sweep_df, bact_counts) {
  per_bact <- sweep_df |>
    filter(!is.na(Sensitivity), !is.na(Specificity)) |>
    mutate(J = Sensitivity + Specificity - 1) |>
    group_by(Time, Bact) |>
    summarise(
      best_J  = {
        jv <- J[!is.na(J)]
        if (length(jv)) max(jv) else NA_real_
      },
      # Range of thresholds achieving best_J (ties included)
      lower   = {
        jv  <- as.numeric(J[!is.na(J)])
        if (!length(jv)) NA_real_
        else { bj <- max(jv); thr <- as.numeric(log10_thresh[!is.na(J) & abs(J-bj)<1e-9]); if (length(thr)) min(thr) else NA_real_ }
      },
      upper   = {
        jv  <- as.numeric(J[!is.na(J)])
        if (!length(jv)) NA_real_
        else { bj <- max(jv); thr <- as.numeric(log10_thresh[!is.na(J) & abs(J-bj)<1e-9]); if (length(thr)) max(thr) else NA_real_ }
      },
      best_Sen = {
        jv <- as.numeric(J[!is.na(J)])
        if (!length(jv)) NA_real_
        else { bj <- max(jv); mean(as.numeric(Sensitivity[!is.na(J) & abs(J-bj)<1e-9]), na.rm=TRUE) }
      },
      best_Spe = {
        jv <- as.numeric(J[!is.na(J)])
        if (!length(jv)) NA_real_
        else { bj <- max(jv); mean(as.numeric(Specificity[!is.na(J) & abs(J-bj)<1e-9]), na.rm=TRUE) }
      },
      achieved = if_else(is.na(best_J), NA_real_, best_J),
      .groups  = "drop"
    )
  list(per_bact = per_bact, common = intersect_ranges(per_bact, bact_counts))
}

# Stage 2b: Sensitivity-first range
# Step 1: find the highest attainable Sensitivity per (Bact, Time) among
#         above-floor thresholds (excluding NA Sensitivity rows).
# Step 2: within thresholds achieving that max Sensitivity, find the range
#         that also maximises Specificity (again excluding NA Specificity).
# Step 3: intersect across bacteria.
sensitivity_first <- function(sweep_df, bact_counts) {
  per_bact <- sweep_df |>
    filter(!is.na(Sensitivity)) |>
    group_by(Time, Bact) |>
    summarise(
      max_sen  = max(Sensitivity, na.rm = TRUE),
      # Among thresholds achieving max Sen, find those also maximising Spe
      # (drop rows where Spe is NA before looking for max Spe)
      best_Spe = {
        ms   <- max(Sensitivity, na.rm = TRUE)
        sub  <- Specificity[abs(Sensitivity - ms) < 1e-9 & !is.na(Specificity)]
        if (length(sub)) max(sub) else NA_real_
      },
      lower    = {
        ms  <- max(as.numeric(Sensitivity), na.rm = TRUE)
        bs  <- { sub <- as.numeric(Specificity[abs(Sensitivity-ms)<1e-9 & !is.na(Specificity)]); if (length(sub)) max(sub) else NA_real_ }
        thr <- if (is.na(bs)) as.numeric(log10_thresh[abs(Sensitivity-ms)<1e-9])
        else           as.numeric(log10_thresh[abs(Sensitivity-ms)<1e-9 & !is.na(Specificity) & abs(Specificity-bs)<1e-9])
        if (length(thr)) min(thr) else NA_real_
      },
      upper    = {
        ms  <- max(as.numeric(Sensitivity), na.rm = TRUE)
        bs  <- { sub <- as.numeric(Specificity[abs(Sensitivity-ms)<1e-9 & !is.na(Specificity)]); if (length(sub)) max(sub) else NA_real_ }
        thr <- if (is.na(bs)) as.numeric(log10_thresh[abs(Sensitivity-ms)<1e-9])
        else           as.numeric(log10_thresh[abs(Sensitivity-ms)<1e-9 & !is.na(Specificity) & abs(Specificity-bs)<1e-9])
        if (length(thr)) max(thr) else NA_real_
      },
      achieved = max(Sensitivity, na.rm = TRUE),
      .groups  = "drop"
    )
  list(per_bact = per_bact, common = intersect_ranges(per_bact, bact_counts))
}

# Stage 2c: Best Accuracy
# Finds the threshold(s) maximising Accuracy = (TP+TN)/(TP+FP+FN+TN) per
# (Bact, Time) among above-floor candidates. Where Accuracy ties, the range
# of tied thresholds is reported. achieved = best Accuracy value.
accuracy_best <- function(sweep_df, bact_counts) {
  per_bact <- sweep_df |>
    filter(!is.na(Accuracy)) |>
    group_by(Time, Bact) |>
    summarise(
      best_Acc = max(as.numeric(Accuracy), na.rm = TRUE),
      lower    = {
        ba  <- max(as.numeric(Accuracy), na.rm = TRUE)
        thr <- as.numeric(log10_thresh[!is.na(Accuracy) & abs(Accuracy - ba) < 1e-9])
        if (length(thr)) min(thr) else NA_real_
      },
      upper    = {
        ba  <- max(as.numeric(Accuracy), na.rm = TRUE)
        thr <- as.numeric(log10_thresh[!is.na(Accuracy) & abs(Accuracy - ba) < 1e-9])
        if (length(thr)) max(thr) else NA_real_
      },
      achieved = max(as.numeric(Accuracy), na.rm = TRUE),
      .groups  = "drop"
    )
  list(per_bact = per_bact, common = intersect_ranges(per_bact, bact_counts))
}

# Main summary function (replaces summarise_best)
# Returns a tidy data frame with one row per (Time, Method, above_floor)
# combining all three Stage-2 approaches run on both the floored and the
# full sweep, so the caller can compare ranges above and below the floor.
# above_floor = TRUE  : range derived from thresholds > bg_hi only
# above_floor = FALSE : range derived from all thresholds (may include
#                       sub-floor values); report these for manual inspection
summarise_best <- function(sweep_df, prod_name, test_label, bg_hi) {
  # Stage 1: filter above background floor
  sweep_floor <- apply_bg_floor(sweep_df, bg_hi)
  n_below_floor <- nrow(sweep_df) - nrow(sweep_floor)
  if (n_below_floor > 0)
    cat(sprintf("    [%s %s] %d threshold steps below background floor removed\n",
                prod_name, test_label, n_below_floor))
  
  # n_total per time-point from the ORIGINAL sweep (before floor).
  # A bacterium is absent from a time-point only if it had no plate-positive
  # replicates there (Sensitivity undefined = not yet growing above 10 CFU/mL).
  bact_counts <- sweep_df |>
    group_by(Time) |>
    summarise(n = n_distinct(Bact), .groups = "drop") |>
    { function(d) setNames(d$n, as.character(d$Time)) }()
  
  # Rund all three optimisations on a given sweep data frame and
  # return a bound data frame tagged with above_floor = af_flag.
  run_all_methods <- function(sw, af_flag) {
    youd <- youden_optimum(sw, bact_counts)
    senf <- sensitivity_first(sw, bact_counts)
    accu <- accuracy_best(sw, bact_counts)
    
    # Helper: mean achieved value per time-point from per_bact
    mean_achieved <- function(pb)
      pb |>
      group_by(Time) |>
      summarise(v = mean(as.numeric(achieved), na.rm = TRUE),
                .groups = "drop") |>
      { function(d) setNames(d$v, as.character(d$Time)) }()
    
    bind_rows(
      youd$common |> mutate(Method = "Youden J",
                            criterion_value =
                              mean_achieved(youd$per_bact)[as.character(Time)]),
      senf$common |> mutate(Method = "Sensitivity-first",
                            criterion_value =
                              mean_achieved(senf$per_bact)[as.character(Time)]),
      accu$common |> mutate(Method = "Best Accuracy",
                            criterion_value =
                              mean_achieved(accu$per_bact)[as.character(Time)])
    ) |>
      mutate(above_floor = af_flag)
  }
  
  fmt_range <- function(lo, hi, valid, af) {
    case_when(
      !valid | is.na(lo) ~ if_else(af, "no range above floor", "no range"),
      lo == hi           ~ sprintf("%.2f", lo),
      TRUE               ~ sprintf("%.2f – %.2f", lo, hi)
    )
  }
  
  bind_rows(
    run_all_methods(sweep_floor, TRUE),   # above-floor results
    run_all_methods(sweep_df,    FALSE)   # full-sweep results (may be sub-floor)
  ) |>
    mutate(
      Product         = prod_name,
      Test            = test_label,
      bg_floor_log10  = round(bg_hi, 2),
      criterion_value = round(as.numeric(criterion_value), 4),
      Range_log10     = fmt_range(common_lower, common_upper, valid, above_floor),
      Range_units     = case_when(
        !valid | is.na(common_lower) ~ "-",
        common_lower == common_upper ~
          formatC(10^common_lower, format = "e", digits = 1),
        TRUE ~ paste0(formatC(10^common_lower, format = "e", digits = 1),
                      " – ",
                      formatC(10^common_upper, format = "e", digits = 1))
      )
    ) |>
    select(Product, Test, Method, above_floor, Time, bg_floor_log10,
           criterion_value, Range_log10, Range_units,
           n_in_intersect, n_above_floor, n_total,
           no_range, problematic, best_achieved, valid)
}

##################################################
# SECTION 6 - MAIN LOOP OVER ALL PRODUCTS
##################################################

all_sweep_charm   <- list()
all_sweep_attune  <- list()
all_best_summary  <- list()
all_avg_tables    <- list()

for (prod_key in names(PRODUCTS)) {
  
  prod_name <- PRODUCTS[[prod_key]]$name
  cat("\n\n", strrep("=", 60), "\n")
  cat("Processing:", prod_name, "\n")
  cat(strrep("=", 60), "\n")
  
  df <- filter(inoculated, Product == prod_name)
  times_here <- levels(droplevels(df$Time))
  
  # Extract pooled blank background stats for this product
  bg_row     <- filter(blank_bg, Product == prod_name)
  bg_charm   <- list(lo = bg_row$charm_bg_lo,  hi = bg_row$charm_bg_hi)
  bg_attune  <- list(lo = bg_row$attune_bg_lo, hi = bg_row$attune_bg_hi)
  
  # 6.1 SIGNAL PLOTS
  
  pdf(file.path("output", paste0(prod_key, "_01_plate_count.pdf")),
      width = 8, height = 5)
  print(plot_plate(df, prod_name))
  dev.off()
  
  # Charm: manufacturer threshold shown now; optimal added after sweep
  pdf(file.path("output", paste0(prod_key, "_02_charm.pdf")),
      width = 8, height = 5)
  print(plot_signal(df, "Charm_adj", "Charm (RLU/test)",
                    hlines       = CHARM_MFR,
                    hline_labels = "300 (mfr.)",
                    prod_name    = prod_name,
                    bg           = bg_charm))
  dev.off()
  
  # Attune Diluted (optimal threshold added after sweep below)
  pdf(file.path("output", paste0(prod_key, "_03_attune_dilut.pdf")),
      width = 8, height = 5)
  print(plot_signal(df, "AttuneDilut_adj", "Attune Diluted (AFU/mL)",
                    prod_name = prod_name,
                    bg        = bg_attune))
  dev.off()
  
  # Per-bacteria vs Plate Count
  pdf(file.path("output", paste0(prod_key, "_04_charm_vs_plate.pdf")),
      width = 11, height = 7)
  print(plot_vs_plate(df, "Charm_adj", "Charm (RLU/test)", prod_name))
  dev.off()
  
  pdf(file.path("output", paste0(prod_key, "_05_attune_vs_plate.pdf")),
      width = 11, height = 7)
  print(plot_vs_plate(df, "AttuneDilut_adj", "Attune Diluted (AFU/mL)", prod_name))
  dev.off()
  
  # 6.2 THRESHOLD SWEEP
  
  # Extract blanks for this product to include in sweep (matched by Bact+Time)
  df_blanks <- filter(blanks, Product == prod_name)
  
  cat("  Running Charm sweep...\n")
  sw_charm  <- run_sweep(df, "Charm_adj",       blanks_df = df_blanks)
  cat("  Running Attune sweep...\n")
  sw_attune <- run_sweep(df, "AttuneDilut_adj",  blanks_df = df_blanks)
  
  all_sweep_charm[[prod_key]]  <- sw_charm  |> mutate(Product = prod_name)
  all_sweep_attune[[prod_key]] <- sw_attune |> mutate(Product = prod_name)
  
  # Write full sweep tables to CSV
  write_csv(sw_charm  |> mutate(Product = prod_name, Test = "Charm"),
            file.path("output", paste0(prod_key, "_sweep_charm.csv")))
  write_csv(sw_attune |> mutate(Product = prod_name, Test = "AttuneDilut"),
            file.path("output", paste0(prod_key, "_sweep_attune.csv")))
  
  # 6.3 AVERAGE-OVER-BACTERIA TABLE (per product × time × threshold)
  
  avg_charm <- sw_charm |>
    group_by(Time, log10_thresh) |>
    summarise(across(c(Sensitivity, Specificity, PPV, NPV,
                       Prevalence, Accuracy, F1, Jaccard, MCC, Youden_J),
                     ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") |>
    mutate(Product = prod_name, Test = "Charm")
  
  avg_attune <- sw_attune |>
    group_by(Time, log10_thresh) |>
    summarise(across(c(Sensitivity, Specificity, PPV, NPV,
                       Prevalence, Accuracy, F1, Jaccard, MCC, Youden_J),
                     ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") |>
    mutate(Product = prod_name, Test = "AttuneDilut")
  
  all_avg_tables[[prod_key]] <- bind_rows(avg_charm, avg_attune)
  
  # Write per-product averaged performance to CSV immediately (no screen print)
  write_csv(
    bind_rows(avg_charm, avg_attune) |>
      select(Product, Test, Time, log10_thresh, everything()) |>
      mutate(across(where(is.numeric), ~ round(.x, 4))),
    file.path("output", paste0(prod_key, "_avg_performance.csv"))
  )
  cat("  Avg performance CSV written for", prod_name, "\n")
  
  # 6.4 BEST THRESHOLD SUMMARIES
  # Pass the background floor (bg$hi, log10 scale) so Stage 1 can filter
  # out thresholds that cannot be distinguished from blank noise.
  
  best_charm  <- summarise_best(sw_charm,  prod_name, "Charm",
                                bg_hi = bg_charm$hi)
  best_attune <- summarise_best(sw_attune, prod_name, "AttuneDilut",
                                bg_hi = bg_attune$hi)
  all_best_summary[[prod_key]] <- bind_rows(best_charm, best_attune)
  
  # Print to console
  cat("\n  --- Best Charm thresholds (two-stage) ---\n")
  best_charm |>
    select(Method, above_floor, Time, criterion_value, Range_log10, Range_units,
           n_in_intersect, n_above_floor, n_total,
           no_range, problematic) |>
    arrange(Method, above_floor, Time) |>
    print(n = Inf)
  
  cat("\n  --- Best Attune thresholds (two-stage) ---\n")
  best_attune |>
    select(Method, above_floor, Time, criterion_value, Range_log10, Range_units,
           n_in_intersect, n_above_floor, n_total,
           no_range, problematic) |>
    arrange(Method, above_floor, Time) |>
    print(n = Inf)
  
  # 6.5 ACCURACY PLOTS (per time-point, one curve per bacterium)
  
  pdf(file.path("output", paste0(prod_key, "_06_accuracy_charm.pdf")),
      width = 9, height = 5)
  for (t in times_here) {
    print(plot_accuracy_curves(sw_charm, t, prod_name, "Charm"))
  }
  dev.off()
  
  pdf(file.path("output", paste0(prod_key, "_07_accuracy_attune.pdf")),
      width = 9, height = 5)
  for (t in times_here) {
    print(plot_accuracy_curves(sw_attune, t, prod_name, "Attune Diluted"))
  }
  dev.off()
  
  # 6.6 MEAN SENSITIVITY & SPECIFICITY PLOTS
  
  pdf(file.path("output", paste0(prod_key, "_08_sensspe_charm.pdf")),
      width = 8, height = 5)
  print(plot_mean_sensspe(sw_charm,  prod_name, "Charm"))
  dev.off()
  
  pdf(file.path("output", paste0(prod_key, "_09_sensspe_attune.pdf")),
      width = 8, height = 5)
  print(plot_mean_sensspe(sw_attune, prod_name, "Attune Diluted"))
  dev.off()
  
  # ── 6.7 ROC CURVES ───────────────────────────────────────────────────────────
  
  pdf(file.path("output", paste0(prod_key, "_10_roc_charm.pdf")),
      width = 6, height = 6)
  print(plot_mean_roc(sw_charm,  prod_name, "Charm"))
  dev.off()
  
  pdf(file.path("output", paste0(prod_key, "_11_roc_attune.pdf")),
      width = 6, height = 6)
  print(plot_mean_roc(sw_attune, prod_name, "Attune Diluted"))
  dev.off()
  
  cat("  Plots and CSVs written for", prod_name, "\n")
}

# 6B - UPDATED SIGNAL PLOTS WITH OPTIMAL THRESHOLDS
# Thresholds are on the log10 scale, chosen manually after inspecting the
# sweep outputs. Each product gets two PDFs: one for Charm, one for Attune.
# The manufacturer reference line (300 RLU) is also shown on Charm plots.
# The background noise band (mean ± BG_N_SD SD of blanks) is shown on both.

# Named vectors indexed by prod_key for easy lookup inside the loop
MANUAL_THRESH <- list(
  charm  = c(Milk         = 1.6,
             Cream        = 0.6,
             WhippingCream = 0.9,
             MedBeverage  = 1.2,
             Chocolate    = 1.3),
  attune = c(Milk         = 5.4,
             Cream        = 4.2,
             WhippingCream = 3.0,
             MedBeverage  = 4.0,
             Chocolate    = 4.7)
)

for (prod_key in names(PRODUCTS)) {
  prod_name <- PRODUCTS[[prod_key]]$name
  df        <- filter(inoculated, Product == prod_name)
  bg_row    <- filter(blank_bg,   Product == prod_name)
  bg_charm  <- list(lo = bg_row$charm_bg_lo,  hi = bg_row$charm_bg_hi)
  bg_attune <- list(lo = bg_row$attune_bg_lo, hi = bg_row$attune_bg_hi)
  
  thr_ch <- MANUAL_THRESH$charm [prod_key]
  thr_at <- MANUAL_THRESH$attune[prod_key]
  
  # Charm: manufacturer reference + manual optimal threshold
  pdf(file.path("output", paste0(prod_key, "_12_charm_manual_threshold.pdf")),
      width = 8, height = 5)
  print(plot_signal(
    df, "Charm_adj", "Charm (RLU/test)",
    hlines       = c(CHARM_MFR,          10^thr_ch),
    hline_labels = c("300 (mfr.)",
                     sprintf("optimal log10=%.1f", thr_ch)),
    prod_name    = prod_name,
    bg           = bg_charm
  ))
  dev.off()
  
  # Attune Diluted: manual optimal threshold only
  pdf(file.path("output", paste0(prod_key, "_13_attune_manual_threshold.pdf")),
      width = 8, height = 5)
  print(plot_signal(
    df, "AttuneDilut_adj", "Attune Diluted (AFU/mL)",
    hlines       = 10^thr_at,
    hline_labels = sprintf("optimal log10=%.1f", thr_at),
    prod_name    = prod_name,
    bg           = bg_attune
  ))
  dev.off()
  
  cat("  Manual-threshold plots written for", prod_name, "\n")
}

##################################################
# 6B2 - CHARM VS PLATE COUNT WITH MANUAL THRESHOLD
# Same per-bacteria boxplot as _04_charm_vs_plate.pdf but with the
# manually chosen log10 threshold overlaid as a horizontal dashed line.
##################################################

# for Charm
for (prod_key in names(PRODUCTS)) {
  prod_name <- PRODUCTS[[prod_key]]$name
  df        <- filter(inoculated, Product == prod_name)
  thr_ch    <- MANUAL_THRESH$charm[prod_key]
  
  # Build the standard vs-plate plot then add the threshold line
  p <- plot_vs_plate(df, "Charm_adj", "Charm (RLU/test)", prod_name) +
    geom_hline(yintercept = thr_ch,
               linetype   = "dashed",
               colour     = "black",
               linewidth  = 0.7) +
    annotate("text",
             x     = max(sort(unique(df$TimeNum))),
             y     = thr_ch,
             label = sprintf("threshold log10=%.1f", thr_ch),
             hjust = 1, vjust = -0.4,
             colour = "black", size = 2.5)
  
  pdf(file.path("output",
                paste0(prod_key, "_15_charm_vs_plate_threshold.pdf")),
      width = 11, height = 7)
  print(p)
  dev.off()
}

# for Attune Dilut
for (prod_key in names(PRODUCTS)) {
  prod_name <- PRODUCTS[[prod_key]]$name
  df        <- filter(inoculated, Product == prod_name)
  thr_ch    <- MANUAL_THRESH$attune[prod_key]
  
  # Build the standard vs-plate plot then add the threshold line
  p <- plot_vs_plate(df, "AttuneDilut_adj", "Attune Diluted (AFU/mL)", prod_name) +
    geom_hline(yintercept = thr_ch,
               linetype   = "dashed",
               colour     = "black",
               linewidth  = 0.7) +
    annotate("text",
             x     = max(sort(unique(df$TimeNum))),
             y     = thr_ch,
             label = sprintf("threshold log10=%.1f", thr_ch),
             hjust = 1, vjust = -0.4,
             colour = "black", size = 2.5)
  
  pdf(file.path("output",
                paste0(prod_key, "_16_attune_vs_plate_threshold.pdf")),
      width = 11, height = 7)
  print(p)
  dev.off()
}

##################################################
# 6C - POD ANALYSIS (Probability of Detection) and Add-on plot with three-test together
# Structure mirrors AOAC PTM Table 4 / document Tables 1-2:
#   Each row = one product × bacteria-group × spike-level combination at 24h.
#
# Three spike level groups (defined by INOCULATED LEVEL, a fixed design
# property, NOT by the observed 24h plate count):
#   ">10"  : bacteria inoculated at > 10 CFU/mL (high positives)
#   "0-10" : bacteria inoculated at 0–10 CFU/mL  (fractional positives)
#   "=0"   : blank replicates (never inoculated; PlatePos always 0)
#
# Design table - which bacteria were inoculated at ">10" per product:
#   UHT Milk       : Geo 1 only
#   Cream          : Geo 3 only
#   Whipping Cream : Geo 1 only
#   Med Beverage   : Geo 1 only
#   Chocolate      : Geo 3, AF 5 Hi, AF 8 Hi
#   All others in every product are "0-10".
#
# N per group (inoculated at 24h, 5 reps per bacterium):
#   ">10"  = 5 × (number of ">10" bacteria in that group for that product)
#   "0-10" = 5 × (number of "0-10" bacteria in that group for that product)
#   "=0"   = number of blank replicates for that group in that product
#            (1 blank per bacterium in the group, pooled across time-points)
#
# X_ref  = number of ">10" or "0-10" inoculated replicates at 24h where
#          plate count > 10 CFU/mL  (i.e. sum of PlatePos at 24h)
# X_test = number of those same replicates where the test signal
#          (Charm or Attune) exceeds the log10 threshold at 24h
#
# CI method:
#   Individual PODs: Clopper-Pearson exact binomial (standard for AOAC PTM).
#   For N <= 5 (small samples at high-positive spike level), the dPOD CI is
#   derived by propagating uncertainty from the individual POD CI limits
#   (Wehling-Wilson approach, avoiding variance collapse at boundary values).
#   For N > 5: Newcombe Stat Med (1998) score-based CI for the difference.
#
# A dPOD CI not containing 0 indicates a statistically significant
# difference at the 5% level (highlighted in the document tables in red).
##################################################

# Clopper-Pearson exact CI for a single proportion
cp_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(pod = NA_real_, lo = NA_real_, hi = NA_real_))
  pod <- x / n
  lo  <- if (x == 0) 0    else qbeta((1 - conf) / 2,     x,     n - x + 1)
  hi  <- if (x == n) 1    else qbeta((1 + conf) / 2, x + 1,     n - x)
  c(pod = pod, lo = lo, hi = hi)
}

# dPOD CI: Newcombe (N>5) or Wehling-Wilson propagation (N<=5)
dpod_ci <- function(x1, n1, x2, n2, conf = 0.95) {
  # Always return a plain named numeric(3) - no nested names, no attributes.
  # as.numeric() strips names from cp_ci() sub-results so arithmetic never
  # produces a named intermediate that breaks c(dPOD=, lo=, hi=) assembly.
  if (n1 == 0 || n2 == 0)
    return(c(dPOD = NA_real_, lo = NA_real_, hi = NA_real_))
  dPOD <- as.numeric(x1/n1 - x2/n2)
  z    <- qnorm((1 + conf) / 2)
  
  if (n1 <= 5 || n2 <= 5) {
    # Wehling-Wilson: propagate uncertainty from individual CP CIs.
    # Extract scalars explicitly with as.numeric() to avoid named-vector
    # arithmetic producing length-1 named results inside lo/hi.
    ci1 <- cp_ci(x1, n1, conf); ci2 <- cp_ci(x2, n2, conf)
    p1  <- as.numeric(ci1["pod"]); l1 <- as.numeric(ci1["lo"]); u1 <- as.numeric(ci1["hi"])
    p2  <- as.numeric(ci2["pod"]); l2 <- as.numeric(ci2["lo"]); u2 <- as.numeric(ci2["hi"])
    lo  <- as.numeric(dPOD - sqrt((p1 - l1)^2 + (u2 - p2)^2))
    hi  <- as.numeric(dPOD + sqrt((u1 - p1)^2 + (p2 - l2)^2))
  } else {
    # Newcombe score-based CI
    w_lo <- function(x, n) {
      p <- x/n
      as.numeric((2*n*p + z^2 - z*sqrt(z^2 + 4*n*p*(1-p))) / (2*(n + z^2)))
    }
    w_hi <- function(x, n) {
      p <- x/n
      as.numeric((2*n*p + z^2 + z*sqrt(z^2 + 4*n*p*(1-p))) / (2*(n + z^2)))
    }
    l1 <- w_lo(x1,n1); u1 <- w_hi(x1,n1)
    l2 <- w_lo(x2,n2); u2 <- w_hi(x2,n2)
    lo  <- as.numeric(dPOD - z * sqrt(l1*(1-l1)/n1 + u2*(1-u2)/n2))
    hi  <- as.numeric(dPOD + z * sqrt(u1*(1-u1)/n1 + l2*(1-l2)/n2))
  }
  c(dPOD = dPOD, lo = lo, hi = hi)
}

# Bacteria group lookup
bact_group <- function(bact) {
  dplyr::case_when(
    bact %in% c("Geo 1","Geo 2","Geo 3","Geo 4")           ~ "Geobacillus",
    bact %in% c("AF 5","AF 6","AF 7","AF 8",
                "AF 5 Hi","AF 8 Hi")                       ~ "Anoxybacillus",
    bact %in% c("BL 1","BL 2")                             ~ "B. licheniformis",
    TRUE                                                   ~ as.character(bact)
  )
}

# Design-time lookup: which bacteria were inoculated at ">10" CFU/mL.
# All other inoculated bacteria in the same product are "0-10".
# This is a FIXED experimental design property, independent of observed data.
HIGH_INOC <- list(
  "UHT Milk"      = c("Geo 1"),
  "Cream"         = c("Geo 3"),
  "Whipping Cream"= c("Geo 1"),
  "Med Beverage"  = c("Geo 1"),
  "Chocolate"     = c("Geo 3","AF 5 Hi","AF 8 Hi")
)

# Assign spike_grp from the inoculated-level design table.
# Returns ">10" for the designated high-inoculum bacteria,
# "0-10" for all other inoculated bacteria (never blank rows).
assign_spike_grp <- function(bact_vec, prod_name) {
  high_bacts <- HIGH_INOC[[prod_name]]
  dplyr::if_else(as.character(bact_vec) %in% high_bacts, ">10", "0-10")
}

# Main POD table: one table per test (Charm or Attune), at 24h.
# Produces rows matching the AOAC document table structure:
#   Product | Isolate group | Spike level | N | X_test | POD_test | CI |
#   X_ref   | POD_ref | CI | dPOD | CI | sig
#
# N is determined by the design (5 reps × number of bacteria in that
# spike_grp × bact_grp combination), NOT from the data row count, ensuring
# N is correct even when some 24h measurements are missing.
# Blanks contribute one row per bacterium per product (pooled across
# time-points), so N_blank = number of distinct Bact × Reps combinations
# in the blank data for that product and bacteria group.
compute_pod_table_24h <- function(inoculated_df, blanks_df, prod_name,
                                  test_col, log_thr) {
  time_val <- "24 Hours"
  
  # Inoculated rows at 24h
  # spike_grp assigned from design table, NOT from observed PlatePos.
  inoc_24 <- inoculated_df |>
    dplyr::filter(Time == time_val) |>
    dplyr::mutate(
      bact_grp  = bact_group(as.character(Bact)),
      spike_grp = assign_spike_grp(Bact, prod_name),
      test_pos  = !is.na(log10s(.data[[test_col]])) &
        log10s(.data[[test_col]]) > log_thr
    )
  
  # Blank rows
  #  Blanks must be evaluated at the SAME time-point as inoculated rows (24h),
  # both for X_test (signal threshold) and for N consistency.
  # N_blank = one row per bacterium in the group at 24h
  # (e.g. N=4 for Geobacillus in Milk: Geo 1, Geo 2, Geo 3, Geo 4 × 1 blank rep).
  # Do NOT pool or deduplicate across time-points: distinct(Bact, Reps) keeps an
  # arbitrary time-point's signal row, giving wrong X_test counts.
  bl <- blanks_df |>
    dplyr::filter(Time == time_val) |>
    dplyr::mutate(
      bact_grp  = bact_group(as.character(Bact)),
      spike_grp = "=0",
      test_pos  = !is.na(log10s(.data[[test_col]])) &
        log10s(.data[[test_col]]) > log_thr,
      PlatePos  = 0L
    )
  
  combined <- dplyr::bind_rows(inoc_24, bl)
  
  # Counts
  # X_ref  = number of replicates in that group with plate count > 10 at 24h.
  #          For ">10" bacteria this should approach N (all positive);
  #          for "0-10" bacteria some may be positive, some not.
  #          For "=0" (blanks) this is always 0 by definition.
  # X_test = number of replicates in that group with test signal > threshold.
  counts <- combined |>
    dplyr::group_by(bact_grp, spike_grp) |>
    dplyr::summarise(
      N      = dplyr::n(),
      X_test = sum(test_pos,       na.rm = TRUE),
      X_ref  = sum(PlatePos == 1L, na.rm = TRUE),
      .groups = "drop"
    )
  
  # CI computation
  # pmap_dfr avoids rowwise() list-column issues across dplyr versions.
  ci_rows <- pmap_dfr(counts, function(bact_grp, spike_grp, N, X_test, X_ref) {
    ct <- cp_ci(X_test, N)
    cr <- cp_ci(X_ref,  N)
    cd <- dpod_ci(X_test, N, X_ref, N)
    tibble(
      bact_grp   = bact_grp,
      spike_grp  = spike_grp,
      N          = N,
      X_test     = X_test,
      X_ref      = X_ref,
      POD_test   = ct[["pod"]],
      CI_test_lo = ct[["lo"]],
      CI_test_hi = ct[["hi"]],
      POD_ref    = cr[["pod"]],
      CI_ref_lo  = cr[["lo"]],
      CI_ref_hi  = cr[["hi"]],
      dPOD       = cd[["dPOD"]],
      dPOD_lo    = cd[["lo"]],
      dPOD_hi    = cd[["hi"]],
      sig        = !is.na(cd[["lo"]]) & !is.na(cd[["hi"]]) &
        (cd[["lo"]] > 0 | cd[["hi"]] < 0)
    )
  })
  
  ci_rows |>
    dplyr::mutate(
      Product         = prod_name,
      Test            = test_col,
      Threshold_log10 = log_thr,
      dplyr::across(c(POD_test, CI_test_lo, CI_test_hi,
                      POD_ref,  CI_ref_lo,  CI_ref_hi,
                      dPOD, dPOD_lo, dPOD_hi), ~ round(.x, 2))
    ) |>
    dplyr::select(Product, Test, Threshold_log10,
                  Isolate     = bact_grp,
                  Spike_level = spike_grp,
                  N, X_test, POD_test, CI_test_lo, CI_test_hi,
                  X_ref,  POD_ref,  CI_ref_lo,  CI_ref_hi,
                  dPOD, dPOD_lo, dPOD_hi, sig) |>
    dplyr::arrange(
      Product,
      factor(Isolate,      levels = c("Geobacillus","Anoxybacillus",
                                      "B. licheniformis")),
      factor(Spike_level,  levels = c(">10","0-10","=0"))
    )
}

# Three-test simultaneous comparison plot
# Shows plate count, Charm and Attune for one bacterium in one product,
# all time-points on the x-axis, log10 scale on the right axes.
# Layout: left y-axis = plate count (CFU/mL), right y-axis = log10 signal.
plot_fig1_style <- function(df_inoc, bact_name, prod_name,
                            thr_charm_log10, thr_attune_log10,
                            charm_thr_range  = NULL,  # c(lo, hi) optional range
                            attune_thr_range = NULL) {
  sub <- df_inoc |>
    filter(as.character(Bact) == bact_name) |>
    mutate(
      TimeNum    = as.numeric(sub(" Hours", "", as.character(Time))),
      log_charm  = log10s(Charm_adj),
      log_attune = log10s(AttuneDilut_adj),
      PlotPlate  = pmin(Plate, PLATE_STORE)
    )
  if (nrow(sub) == 0) {
    message("No data for ", bact_name, " in ", prod_name); return(NULL)
  }
  
  time_breaks <- sort(unique(sub$TimeNum))
  
  # Scaling factor: map log10 signal range onto plate count axis
  plate_max  <- PLATE_STORE
  log_lo     <- 0; log_hi <- 8
  scale_fac  <- plate_max / log_hi
  
  # Jitter width proportional to smallest time gap so points do not
  # visually collapse adjacent time-points on the continuous axis.
  min_gap    <- min(diff(time_breaks))
  jit_width  <- min_gap * 0.06   # 6% of the smallest gap
  
  p <- sub |>
    ggplot(aes(x = TimeNum)) +
    # Connecting lines per replicate (drawn first, behind points)
    geom_line(aes(y = PlotPlate, group = Reps),
              colour = "black", alpha = 0.25, linewidth = 0.4) +
    geom_line(aes(y = log_charm  * scale_fac, group = Reps),
              colour = "#377EB8", alpha = 0.25, linewidth = 0.4) +
    geom_line(aes(y = log_attune * scale_fac, group = Reps),
              colour = "#E41A1C", alpha = 0.25, linewidth = 0.4) +
    # Plate count (left axis) - black
    geom_jitter(aes(y = PlotPlate, shape = "Plate count"),
                width = jit_width, size = 2, colour = "black", alpha = 0.80) +
    geom_hline(yintercept = PLATE_POS, linetype = "dotted",
               colour = "black", linewidth = 0.7) +
    # Charm (right axis, scaled) - blue
    geom_jitter(aes(y = log_charm * scale_fac, shape = "ATP Bioluminescence"),
                width = jit_width, size = 2, colour = "#377EB8", alpha = 0.80) +
    geom_hline(yintercept = thr_charm_log10 * scale_fac,
               linetype = "dotted", colour = "#377EB8", linewidth = 0.7) +
    # Attune (right axis, scaled) - red
    geom_jitter(aes(y = log_attune * scale_fac, shape = "Flow Cytometry"),
                width = jit_width, size = 2, colour = "#E41A1C", alpha = 0.80) +
    geom_hline(yintercept = thr_attune_log10 * scale_fac,
               linetype = "dotted", colour = "#E41A1C", linewidth = 0.7) +
    scale_x_continuous(breaks = time_breaks,
                       labels = paste0(time_breaks, "h"),
                       # Expand slightly so edge points are not clipped
                       expand = expansion(mult = 0.05)) +
    scale_y_continuous(
      name   = "Plate count (CFU/mL)",
      breaks = seq(0, plate_max, 50),
      labels = c(seq(0, 250, 50), "≥300"),
      limits = c(0, plate_max),
      sec.axis = sec_axis(
        transform = ~ . / scale_fac,
        name   = "log10(signal)",
        breaks = seq(0, log_hi, 1)
      )
    ) +
    scale_shape_manual(values = c("Plate count" = 16,
                                  "ATP Bioluminescence" = 17,
                                  "Flow Cytometry" = 15)) +
    labs(title  = paste0(prod_name, " - ", bact_name),
         x      = "Pre-incubation time (h)",
         shape  = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position  = "bottom",
          axis.title.y.right = element_text(colour = "grey40"),
          axis.text.y.right  = element_text(colour = "grey40"))
  
  # Optional threshold ranges (shaded bands)
  if (!is.null(charm_thr_range))
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = charm_thr_range[1] * scale_fac,
                      ymax = charm_thr_range[2] * scale_fac,
                      fill = "#377EB8", alpha = 0.08)
  if (!is.null(attune_thr_range))
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = attune_thr_range[1] * scale_fac,
                      ymax = attune_thr_range[2] * scale_fac,
                      fill = "#E41A1C", alpha = 0.08)
  p
}

# Run POD + 3 tests together plots for all products
all_pod_charm  <- list()
all_pod_attune <- list()

for (prod_key in names(PRODUCTS)) {
  prod_name  <- PRODUCTS[[prod_key]]$name
  df         <- filter(inoculated, Product == prod_name)
  df_blanks  <- filter(blanks,     Product == prod_name)
  thr_ch     <- MANUAL_THRESH$charm [prod_key]
  thr_at     <- MANUAL_THRESH$attune[prod_key]
  
  # POD tables (one per test)
  pod_ch <- compute_pod_table_24h(df, df_blanks, prod_name,
                                  "Charm_adj",       thr_ch)
  pod_at <- compute_pod_table_24h(df, df_blanks, prod_name,
                                  "AttuneDilut_adj",  thr_at)
  
  all_pod_charm [[prod_key]] <- pod_ch
  all_pod_attune[[prod_key]] <- pod_at
  
  write_csv(pod_ch, file.path("output",
                              paste0(prod_key, "_pod_charm_24h.csv")))
  write_csv(pod_at, file.path("output",
                              paste0(prod_key, "_pod_attune_24h.csv")))
  cat("  POD CSVs written for", prod_name, "\n")
  
  # 3-test plots: one PDF per product, one page per bacterium
  bacts_here <- levels(droplevels(df$Bact))
  pdf(file.path("output", paste0(prod_key, "_14_fig1_comparison.pdf")),
      width = 7, height = 5)
  for (b in bacts_here) {
    p <- plot_fig1_style(df, b, prod_name,
                         thr_charm_log10  = thr_ch,
                         thr_attune_log10 = thr_at)
    if (!is.null(p)) print(p)
  }
  dev.off()
  cat("  Fig-1 comparison plots written for", prod_name, "\n")
}

# Master POD tables
write_csv(bind_rows(all_pod_charm),
          "output/master_pod_charm_24h.csv")
write_csv(bind_rows(all_pod_attune),
          "output/master_pod_attune_24h.csv")

##################################################
# 6D - ADDITIONAL CHARM POD TABLES AT MANUFACTURER THRESHOLDS
# Charm 150 RLU  (log10 = log10(150) ≈ 2.18, lower "suspect" boundary)
# Charm 300 RLU  (log10 = log10(300) ≈ 2.48, upper "fail" / manufacturer threshold)
# Same compute_pod_table_24h() function, only the threshold changes.
##################################################

CHARM_THR_150 <- log10(150)   # ≈ 2.176
CHARM_THR_300 <- log10(300)   # ≈ 2.477

all_pod_charm_150 <- list()
all_pod_charm_300 <- list()

for (prod_key in names(PRODUCTS)) {
  prod_name <- PRODUCTS[[prod_key]]$name
  df        <- filter(inoculated, Product == prod_name)
  df_blanks <- filter(blanks,     Product == prod_name)
  
  pod_150 <- compute_pod_table_24h(df, df_blanks, prod_name,
                                   "Charm_adj", CHARM_THR_150)
  pod_300 <- compute_pod_table_24h(df, df_blanks, prod_name,
                                   "Charm_adj", CHARM_THR_300)
  
  all_pod_charm_150[[prod_key]] <- pod_150
  all_pod_charm_300[[prod_key]] <- pod_300
  
  write_csv(pod_150, file.path("output",
                               paste0(prod_key, "_pod_charm_150RLU_24h.csv")))
  write_csv(pod_300, file.path("output",
                               paste0(prod_key, "_pod_charm_300RLU_24h.csv")))
  cat("  Charm 150/300 RLU POD CSVs written for", prod_name, "\n")
}

write_csv(bind_rows(all_pod_charm_150), "output/master_pod_charm_150RLU_24h.csv")
write_csv(bind_rows(all_pod_charm_300), "output/master_pod_charm_300RLU_24h.csv")

##################################################
# SECTION 7 - COMBINED OUTPUTS ACROSS ALL PRODUCTS
##################################################

# 7.1 Write master best-threshold table
# summarise_best() already formats Range_log10 and Range_units; just bind
# all products and select the columns we want.

master_best <- bind_rows(all_best_summary) |>
  select(Product, Test, Method, above_floor, Time, bg_floor_log10,
         criterion_value, Range_log10, Range_units,
         n_in_intersect, n_above_floor, n_total,
         no_range, problematic, best_achieved, valid)

write_csv(master_best, "output/master_best_thresholds.csv")

cat("\n\n", strrep("=", 60), "\n")
cat("MASTER BEST THRESHOLD TABLE (Youden J + Sensitivity-first)\n")
cat(strrep("=", 60), "\n")
print(master_best, n = Inf)

# 7.2 Write master average-performance table

# Write master best-threshold summary (all products, all methods) to CSV
master_best_all <- bind_rows(all_best_summary)
write_csv(master_best_all, "output/master_best_thresholds_all.csv")

# master_avg combines all products; individual per-product CSVs were already
# written inside the loop. This master file is for cross-product queries.
master_avg <- bind_rows(all_avg_tables) |>
  mutate(across(where(is.numeric), ~ round(.x, 4))) |>
  select(Product, Test, Time, log10_thresh, everything())

write_csv(master_avg, "output/master_avg_performance.csv")

# 7.3 Conclusions printed to console
# For each product and test, report:
#   - Background floor applied
#   - Youden J optimum per time-point
#   - Sensitivity-first range per time-point
#   - Any bacteria with no valid range (flagged)

cat("\n\n", strrep("=", 60), "\n")
cat("CONCLUSIONS\n")
cat(strrep("=", 60), "\n\n")

for (prod_key in names(PRODUCTS)) {
  prod_name <- PRODUCTS[[prod_key]]$name
  cat("--", prod_name, "--\n")
  
  for (test in c("Charm", "AttuneDilut")) {
    cat("  Test:", test, "\n")
    sub <- filter(bind_rows(all_best_summary[[prod_key]]), Test == test)
    if (nrow(sub) == 0) { cat("    No results\n"); next }
    
    cat(sprintf("    Background floor (log10): %.2f\n",
                sub$bg_floor_log10[1]))
    
    for (meth in c("Youden J", "Sensitivity-first", "Best Accuracy")) {
      cat("   ", meth, ":\n")
      # Show above-floor rows first, then below-floor rows (for manual comparison)
      for (af in c(TRUE, FALSE)) {
        rows <- filter(sub, Method == meth, above_floor == af)
        if (nrow(rows) == 0) next
        cat(sprintf("      [%s background floor]\n",
                    if (af) "ABOVE" else "BELOW"))
        for (i in seq_len(nrow(rows))) {
          r <- rows[i, ]
          cat(sprintf("        %s: %s  (criterion=%.4f)  [%d in intersect / %d above floor / %d total bacteria]",
                      r$Time, r$Range_log10, r$criterion_value,
                      r$n_in_intersect, r$n_above_floor, r$n_total))
          if (!is.na(r$no_range) && nchar(r$no_range) > 0)
            cat(sprintf("  NO RANGE: %s", r$no_range))
          if (!is.na(r$problematic) && nchar(r$problematic) > 0)
            cat(sprintf("  PROBLEMATIC: %s", r$problematic))
          if (!is.na(r$best_achieved) && nchar(r$best_achieved) > 0)
            cat(sprintf("  best achieved: %s", r$best_achieved))
          cat("\n")
        }
      }
    }
  }
  cat("\n")
}
