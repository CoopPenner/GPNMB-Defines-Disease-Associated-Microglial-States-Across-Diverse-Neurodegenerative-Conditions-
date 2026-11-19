suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(ggplot2); library(patchwork)
  library(MASS); library(lme4); library(lmerTest); library(ggseg); library(ggnewscale); library(data.table)
})
select <- dplyr::select

# ---------------------------------------------------------- CONFIG ------------------------------------------------------------
box <- c(file.path(Sys.getenv("USERPROFILE"), "Box"), "~/Library/CloudStorage/Box-Box", "~/Box")
box <- box[dir.exists(box)][1]

RUN_FIG7 <- Sys.getenv("RUN_FIG7", "TRUE") == "TRUE"
RUN_FIG8 <- Sys.getenv("RUN_FIG8", "TRUE") == "TRUE"

# input / output locations: everything lives in Box > Cooper + Marc > GPNMB microglia paper > Data > Clinical
clinical    <- file.path(box, "Cooper + Marc", "GPNMB microglia paper", "Data", "Clinical")
in_seaad    <- file.path(clinical, "inputs", "SEA-AD")                 # Cognitive_Scores.xlsx
in_indd     <- file.path(clinical, "inputs", "INDD")                   # INDD autopsy + regional pathology workbooks
built       <- file.path(clinical, "outputs")                          # SEAAD/AMPPD_dataset.xlsx and *_cells.csv.gz from the two *_Clinical_Build.R scripts
fig7_tables <- file.path(clinical, "outputs", "Figure7_tables")        # stats / donor tables made by Make_Figure7_tables.R (run it first)
out_root <- Sys.getenv("FINAL_FIG_OUT", unset = file.path(clinical, "outputs"))

# required inputs; a part whose inputs are missing is skipped with a message
inputs <- list(
  fig7 = c(file.path(fig7_tables, c("MTL_vs_cortex_metric_stats.csv", "AMPPD_region_groups_metric_stats.csv", "SEAAD_donor_metrics_MTG_DFC_MEC.csv",
                                    "SEAAD_metric_vs_pathology_stats.csv", "SEAAD_metric_vs_cognition_stats.csv", "AMPPD_donor_region_metrics_z.csv",
                                    "AMPPD_metric_vs_pathology_per_region_stats.csv")),
           file.path(in_seaad, "Cognitive_Scores.xlsx")),
  fig8_microglia = file.path(built, c("SEAAD_dataset.xlsx", "AMPPD_dataset.xlsx", "SEAAD_cells.csv.gz", "AMPPD_cells.csv.gz")),
  fig8_indd = file.path(in_indd, c("allPathSpreadData.xlsx", "allPathData_Cerad_Ammended.xlsx")))
missing <- lapply(inputs, function(f) f[!file.exists(f)])
for (k in names(missing)) if (length(missing[[k]])) message("MISSING inputs (", k, " will be skipped):
  ", paste(missing[[k]], collapse = "
  "))
has <- sapply(missing, length) == 0
dir.create(file.path(out_root, "Figure7"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "Figure8"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "Supplement"), recursive = TRUE, showWarnings = FALSE)
out7 <- file.path(out_root, "Figure7"); out8 <- file.path(out_root, "Figure8"); outS <- file.path(out_root, "Supplement")

# -------------------------------------------------------- SHARED HELPERS ------------------------------------------------------
num <- function(v) suppressWarnings(as.numeric(as.character(v)))
inddid <- function(v) sub("\\.?0+$", "", sprintf("%.2f", as.numeric(v)))       # INQuery ID; the decimal suffix marks different people
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), signif(p, 2))
to_dose <- c(TT = 2, CT = 1, CC = 0)
dose_of <- function(x) unname(to_dose[sub(".*(TT|CT|CC).*", "\\1", x)])
geno_levels <- c("GG", "AG", "AA")
geno_colours <- c(GG = "#DEEBF7", AG = "#9ECAE1", AA = "#4F7CAC")

score_map <- c("0" = 0, "none" = 0, "rare" = 0.5, "1" = 1, "1+" = 1, "mild" = 1, "2" = 2, "2+" = 2, "moderate" = 2, "3" = 3, "3+" = 3, "severe" = 3)
score_path <- function(x) unname(score_map[tolower(trimws(as.character(x)))])

panel_theme <- function(base = 7.5) theme_classic(base_size = base) +
  theme(axis.text = element_text(colour = "black", size = base - 1), axis.title = element_text(size = base - 0.5),
        plot.title = element_text(hjust = 0.5, size = base, face = "bold"), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), strip.background = element_blank(), strip.text = element_text(face = "bold", size = base))
save_pdf <- function(plot, file, w, h) ggsave(file, plot, width = w, height = h, units = "mm", device = cairo_pdf)


# =====================================================================================================================
# FIGURE 7   
# =====================================================================================================================
make_figure7 <- function() {
  set.seed(1)   # fixed jitter of the dots, so the panels are identical on every run
  blue <- "#2a78d6"; grey <- "#f0efec"; red <- "#c93c3c"
  dementia_colours <- c("No dementia" = "#5B8DC9", "Dementia" = "#D9534F")
  pd_colours <- c("Control" = "#5B8DC9", "PD" = "#D9534F")
  small_theme <- panel_theme(7.5)

  # ---- 7 (heatmap): AD (SEA-AD) and PD (AMP-PD) shared-OR heatmap, MTL / cortex and brainstem / basal ganglia / cortex ----
  ad_stats <- read.csv(file.path(fig7_tables, "MTL_vs_cortex_metric_stats.csv")) %>%
    filter(outcome %in% c("Braak", "Thal", "CERAD", "Memory", "Language")) %>%
    mutate(block = paste(ifelse(group == "MTL", "MTL", "Cortex"), ifelse(family == "pathology", "Pathology", "Cognition"), sep = "\n"),
           block = factor(block, c("MTL\nPathology", "MTL\nCognition", "Cortex\nPathology", "Cortex\nCognition")),
           outcome = factor(outcome, c("Braak", "Thal", "CERAD", "Memory", "Language")),
           row = factor(metric, c("HOMEO", "TREM2", "GPNMB", "GNDAM"), c("Homeostatic 1", "TREM2", "GPNMB", "GNDAM")),
           label = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p <= 0.1, sprintf("%.2f", p), "ns")))),
           bold = ifelse(p < 0.05, "bold", "plain"))
  pd_stats <- read.csv(file.path(fig7_tables, "AMPPD_region_groups_metric_stats.csv")) %>%
    filter(outcome %in% c("BraakLB", "BraakNFT")) %>%
    mutate(block = factor(group, c("Brainstem", "Basal ganglia", "Cortex"), c("Brainstem", "Basal\nganglia", "Cortex")),
           outcome = factor(outcome, c("BraakLB", "BraakNFT"), c("Braak LB", "Braak NFT")),
           row = factor(metric, c("HOMEO", "TREM2", "GPNMB", "GNDAM"), c("Homeostatic 1", "TREM2", "GPNMB", "GNDAM")),
           label = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p <= 0.1, sprintf("%.2f", p), "ns")))),
           bold = ifelse(p < 0.05, "bold", "plain"))
  heatmap_theme <- theme_minimal(base_size = 6.5) +
    theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(face = "bold", size = 6, lineheight = 0.9), plot.title = element_text(face = "bold", size = 6.5, hjust = 0.5),
          legend.title = element_text(size = 5.5), legend.text = element_text(size = 5), legend.key.height = unit(2, "mm"),
          legend.key.width = unit(6, "mm"), legend.position = "right")
  or_scale <- scale_fill_gradient2(low = blue, mid = grey, high = red, midpoint = 0, limits = c(-2, 2), breaks = c(-2, 0, 2), labels = c("0.25", "1", "4"),
                                   oob = scales::squish, name = "OR per SD", guide = guide_colourbar(direction = "horizontal", title.position = "top", order = 1))
  heatmap_ad <- ggplot(ad_stats, aes(outcome, row)) +
    geom_tile(data = filter(ad_stats, family == "pathology"), aes(fill = log2(est)), colour = "white", linewidth = 0.6) + or_scale + new_scale_fill() +
    geom_tile(data = filter(ad_stats, family == "cognition"), aes(fill = -est), colour = "white", linewidth = 0.6) +
    scale_fill_gradient2(low = blue, mid = grey, high = red, midpoint = 0, limits = c(-0.025, 0.025), breaks = c(-0.025, 0, 0.025),
                         labels = c("slower", "0", "faster"), oob = scales::squish, name = "Cognitive decline",
                         guide = guide_colourbar(direction = "horizontal", title.position = "top", order = 2)) +
    geom_text(aes(label = label, fontface = bold), size = 2) + facet_grid(. ~ block, scales = "free_x", space = "free_x") +
    scale_x_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL, title = "Alzheimer's disease (SEA-AD)") + heatmap_theme   # (title fixed: the draft said AMP-PD)
  heatmap_pd <- ggplot(pd_stats, aes(outcome, row)) + geom_tile(aes(fill = log2(OR)), colour = "white", linewidth = 0.6) + or_scale +
    geom_text(aes(label = label, fontface = bold), size = 2) + facet_grid(. ~ block, scales = "free_x", space = "free_x") +
    scale_x_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL, title = "Parkinson's disease (AMP-PD)") + heatmap_theme + theme(axis.text.y = element_blank())
  save_pdf(heatmap_ad + heatmap_pd + plot_layout(widths = c(10, 6), guides = "collect"), file.path(out7, "Panel_heatmap_AD_PD_sharedOR.pdf"), 170, 55)

  # ---- 7 (SEA-AD): microglial GPNMB and TREM2 vs Braak / Thal / CERAD and vs memory / language decline ----
  seaad <- read.csv(file.path(fig7_tables, "SEAAD_donor_metrics_MTG_DFC_MEC.csv"))
  slopes <- read_excel(file.path(in_seaad, "Cognitive_Scores.xlsx")) %>% group_by(donor = `Donor ID`) %>%
    summarise(Memory = na.omit(slope_zmem0)[1], Language = na.omit(slope_zlan0)[1], .groups = "drop")
  seaad <- seaad %>% left_join(slopes, by = "donor") %>%
    mutate(Braak = factor(Braak, c("0", "I", "II", "III", "IV", "V", "VI")), Thal = factor(Thal, c("0", "1", "2", "3", "4", "5")),
           CERAD = factor(CERAD, c("Absent", "Sparse", "Moderate", "Frequent")), dementia = factor(dementia, c("No dementia", "Dementia")))
  seaad_path_stats <- read.csv(file.path(fig7_tables, "SEAAD_metric_vs_pathology_stats.csv"))
  seaad_cog_stats <- read.csv(file.path(fig7_tables, "SEAAD_metric_vs_cognition_stats.csv"))
  stage_titles <- c(Braak = "Tau (Braak stage)", Thal = "Amyloid (Thal phase)", CERAD = "Plaques (CERAD)")
  stage_colours <- list(Braak = c("#EAF1F8", "#08306B"), Thal = c("#FCE3D6", "#99000D"), CERAD = c("#EDF8E9", "#006D2C"))
  seaad_plots <- list()
  for (metric in c("GPNMB", "TREM2")) {
    for (stage in c("Braak", "Thal", "CERAD")) {
      stats <- filter(seaad_path_stats, metric == !!metric, outcome == stage)
      plot_data <- seaad %>% filter(!is.na(.data[[stage]])) %>% mutate(stage_level = droplevels(.data[[stage]]), value = .data[[paste0(metric, "_z")]])
      box_colours <- colorRampPalette(stage_colours[[stage]])(nlevels(plot_data$stage_level))
      p <- ggplot(plot_data, aes(stage_level, value)) +
        geom_boxplot(aes(fill = stage_level), outlier.shape = NA, width = 0.7, linewidth = 0.3, show.legend = FALSE) +
        geom_jitter(aes(colour = dementia), width = 0.18, height = 0, size = 0.7) +
        scale_fill_manual(values = box_colours) + scale_colour_manual(values = dementia_colours, name = NULL) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.25, size = 2.3, fontface = "italic", label = paste0("OR ", sprintf("%.2f", stats$OR), ", FDR ", signif(stats$p_fdr, 2))) +
        coord_cartesian(ylim = c(-2.6, 3.8)) +
        labs(x = NULL, y = if (stage == "Braak") paste("Microglial", metric, "(z)") else NULL, title = if (metric == "GPNMB") stage_titles[[stage]] else NULL) + small_theme
      if (stage == "CERAD") p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1))
      seaad_plots[[paste(metric, stage)]] <- p
    }
    for (domain in c("Memory", "Language")) {
      stats <- filter(seaad_cog_stats, metric == !!metric, domain == !!domain)
      plot_data <- seaad %>% filter(!is.na(.data[[domain]])) %>% mutate(value = .data[[paste0(metric, "_z")]], slope = .data[[domain]])
      p <- ggplot(plot_data, aes(value, slope)) + geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
        geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey80", linewidth = 0.45, show.legend = FALSE) +
        geom_point(aes(colour = dementia), size = 0.8) + scale_colour_manual(values = dementia_colours, name = NULL) +
        guides(colour = guide_legend(override.aes = list(size = 2))) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        annotate("text", x = Inf, y = Inf, hjust = 1.03, vjust = 1.25, size = 2.3, fontface = "italic", label = paste0("beta ", sprintf("%.3f", stats$beta), ", FDR ", signif(stats$p_fdr, 2))) +
        labs(x = paste("Microglial", metric, "(z)"), y = if (domain == "Memory") "Annual slope (z units)" else NULL, title = if (metric == "GPNMB") domain else NULL) + small_theme
      seaad_plots[[paste(metric, domain)]] <- p
    }
  }
  save_pdf(wrap_plots(seaad_plots, nrow = 2, widths = c(7, 6, 4.2, 5.2, 5.2), guides = "collect") & theme(legend.position = "bottom"),
           file.path(out7, "Panel_SEAAD_GPNMB_TREM2_pathology_cognition.pdf"), 210, 118)

  # ---- 7 (AMP-PD): microglial GPNMB and TREM2 vs Braak LB and Braak NFT in brainstem (DMNX) and primary visual cortex (PVC) ----
  amppd <- read.csv(file.path(fig7_tables, "AMPPD_donor_region_metrics_z.csv")) %>% filter(region %in% c("DMNX", "PVC")) %>% mutate(dx = factor(dx, c("Control", "PD")))
  amppd_stats <- read.csv(file.path(fig7_tables, "AMPPD_metric_vs_pathology_per_region_stats.csv"))
  region_titles <- c(DMNX = "Brainstem (DMNX)", PVC = "Primary visual cortex (PVC)")
  amppd_plots <- list()
  for (metric in c("GPNMB", "TREM2")) for (region in c("DMNX", "PVC")) for (stage in c("BraakLB", "BraakNFT")) {
    stats <- filter(amppd_stats, region == !!region, metric == !!metric, outcome == stage)
    plot_data <- amppd %>% filter(region == !!region, !is.na(.data[[stage]])) %>% mutate(value = pmin(.data[[paste0(metric, "_z")]], 3.5))
    if (stage == "BraakLB") { plot_data$stage_level <- factor(plot_data$BraakLB, 0:6, as.character(0:6)); box_colours <- colorRampPalette(c("#FFF7BC", "#8C2D04"))(7); title <- "Lewy bodies (Braak LB)"
    } else { plot_data$stage_level <- factor(plot_data$BraakNFT, 0:3, c("0", "I", "II", "III")); box_colours <- colorRampPalette(c("#EAF1F8", "#08306B"))(4); title <- "Tau (Braak NFT)" }
    amppd_plots[[paste(metric, region, stage)]] <- ggplot(plot_data, aes(stage_level, value)) +
      geom_boxplot(aes(fill = stage_level), outlier.shape = NA, width = 0.7, linewidth = 0.3, show.legend = FALSE) +
      geom_jitter(aes(colour = dx), width = 0.18, height = 0, size = 0.7) + scale_fill_manual(values = box_colours, drop = FALSE) + scale_x_discrete(drop = FALSE) +
      scale_colour_manual(values = pd_colours, name = NULL) + guides(colour = guide_legend(override.aes = list(size = 2))) +
      annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.25, size = 2.3, fontface = "italic", label = paste0("OR ", sprintf("%.2f", stats$OR), ", FDR ", signif(stats$p_fdr, 2))) +
      coord_cartesian(ylim = c(-3.3, 4.4)) +
      labs(x = NULL, y = if (region == "DMNX" && stage == "BraakLB") paste("Microglial", metric, "(z)") else NULL,
           title = if (metric == "GPNMB") paste(region_titles[[region]], title, sep = "\n") else NULL) + small_theme
  }
  save_pdf(wrap_plots(amppd_plots, nrow = 2, widths = c(7, 4.5, 7, 4.5), guides = "collect") & theme(legend.position = "bottom"),
           file.path(out7, "Panel_AMPPD_GPNMB_TREM2_pathology_regions_DMNX_vs_PVC.pdf"), 210, 118)
  message("Figure 7 done -> ", out7)
}

# =====================================================================================================================
# FIGURE 8
# =====================================================================================================================
# ---- data shared by the INDD panels ----
load_indd <- function() {
  spread <- read_excel(file.path(in_indd, "allPathSpreadData.xlsx"), guess_max = 5000) %>%
    mutate(INDDID = inddid(INDDID)) %>% distinct(INDDID, .keep_all = TRUE)
  wide <- read_excel(file.path(in_indd, "allPathData_Cerad_Ammended.xlsx"), guess_max = 5000)
  list(spread = spread, wide = wide)
}

# ---- 8: Genotype_Expression_SEAAD / Genotype_Expression_AMPPD ------------------------------------------------------------------------
# left: microglial GPNMB (z) by genotype (donors with >= 10 nuclei);  right: % GPNMB+ nuclei by cell state (GNDAM = brown module; Homeostatic = yellow + black)
fig8_expression <- function() {
  read_donors <- function(f) read_excel(file.path(built, f), sheet = "donor") %>% filter(!is.na(rs199347), n_nuclei >= 10) %>%
    mutate(sex = factor(sex), GPNMB_z = as.numeric(scale(GPNMB)))
  cohorts <- list("SEA-AD" = read_donors("SEAAD_dataset.xlsx"), "AMP-PD" = read_donors("AMPPD_dataset.xlsx"))
  covariates <- c("SEA-AD" = "age + sex + pmi + apoe_e4", "AMP-PD" = "age + sex + PD")
  # per-nucleus module assignment (all regions) and GPNMB expression from the build outputs
  read_cells <- function(f, module_col) fread(file.path(built, f), select = c("donor", "region", module_col, "GPNMB")) %>% as.data.frame() %>% rename(module = all_of(module_col))
  cells <- bind_rows(read_cells("SEAAD_cells.csv.gz", "module_allregions") %>% mutate(cohort = "SEA-AD"),
                     read_cells("AMPPD_cells.csv.gz", "module") %>% mutate(cohort = "AMP-PD"))
  cells$state <- ifelse(cells$module == "brown", "GNDAM", ifelse(cells$module %in% c("yellow", "black"), "Homeostatic", "Other"))
  state_colours <- c(GNDAM = "#8B5A2B", Homeostatic = "#E8B83A", Other = "#9E9E9E")
  files <- c("SEA-AD" = "Genotype_Expression_SEAAD", "AMP-PD" = "Genotype_Expression_AMPPD")
  theme_e <- panel_theme(7) + theme(plot.subtitle = element_text(hjust = 0.5, size = 6), legend.key.size = unit(2.5, "mm"),
                                    legend.text = element_text(size = 5.5), legend.title = element_blank())
  stats_out <- data.frame()
  for (cohort in names(cohorts)) {
    d <- cohorts[[cohort]]
    fit <- summary(lm(as.formula(paste("GPNMB_z ~ rs199347_A_dosage +", covariates[[cohort]])), data = d))$coefficients["rs199347_A_dosage", ]
    d$genotype <- factor(d$rs199347, geno_levels)
    stats_out <- rbind(stats_out, data.frame(cohort = cohort, measure = "GPNMB z-score", state = "All microglia", beta = fit[[1]], fold_AA_vs_GG = NA, p = fit[[4]]))
    left <- ggplot(d, aes(genotype, GPNMB_z, fill = genotype)) + geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
      geom_jitter(width = 0.15, height = 0, size = 0.7, alpha = 0.7) + scale_fill_manual(values = geno_colours, guide = "none") +
      labs(x = "rs199347 genotype", y = "Microglial GPNMB (z-score)", title = cohort, subtitle = paste0("\u03b2 = ", sprintf("%.2f", fit[[1]]), " per A, P = ", fmt_p(fit[[4]]))) + theme_e
    pos <- cells %>% filter(cohort == !!cohort, !is.na(GPNMB)) %>% group_by(donor, state) %>% summarise(pct = 100 * mean(GPNMB > 0), n = n(), .groups = "drop") %>%
      filter(n >= 20) %>% inner_join(d %>% select(donor, rs199347, rs199347_A_dosage, age, sex, any_of(c("pmi", "apoe_e4", "PD"))), by = "donor") %>%
      mutate(genotype = factor(rs199347, geno_levels))
    labels <- c()
    for (s in names(state_colours)) {
      ds <- filter(pos, state == s)
      f <- summary(lm(as.formula(paste("pct ~ rs199347_A_dosage +", covariates[[cohort]])), data = ds))$coefficients["rs199347_A_dosage", ]
      fold <- mean(ds$pct[ds$genotype == "AA"]) / mean(ds$pct[ds$genotype == "GG"])
      stats_out <- rbind(stats_out, data.frame(cohort = cohort, measure = "GPNMB+ microglia (%)", state = s, beta = f[[1]], fold_AA_vs_GG = fold, p = f[[4]]))
      labels[s] <- sprintf("%s  \u00d7%.1f, P %s", s, fold, fmt_p(f[[4]]))
    }
    means <- pos %>% group_by(state, genotype) %>% summarise(mean = mean(pct), se = sd(pct) / sqrt(n()), .groups = "drop")
    right <- ggplot(means, aes(genotype, mean, colour = state, group = state)) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12, linewidth = 0.35) + geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
      scale_colour_manual(values = state_colours, labels = labels) + labs(x = "rs199347 genotype", y = "GPNMB+ microglia (%)", title = cohort) +
      theme_e + theme(legend.position = "bottom", legend.direction = "vertical")
    save_pdf(left | right, file.path(out8, paste0(files[[cohort]], ".pdf")), 180, 90)
  }
  write.csv(stats_out, file.path(out8, "Genotype_Expression_stats.csv"), row.names = FALSE)
}

# ---- AD pathology by genotype: INDD primary AD (main) and SEA-AD donors with AD pathology (supplement) -------------------------------------
# ordinal logistic regression, stage ~ A + age + sex + PMI + APOE e4 (number of alleles); amyloid = A score grouped from Thal (0 | 1-2 | 3 | 4-5)
fig8_ad_pathology <- function(indd) {
  outcomes <- list(
    Braak   = list(title = "Tau Progression\n(Braak Stage)", labels = c("0", "I", "II", "III", "IV", "V", "VI"), pal = c("#EAF1F8", "#08306B")),
    Amyloid = list(title = "Amyloid\n(A score, grouped Thal)", labels = c("A0", "A1", "A2", "A3"), pal = c("#FCE3D6", "#99000D")),
    Plaques = list(title = "Plaque Density\n(CERAD Score)", labels = c("Absent", "Sparse", "Moderate", "Frequent"), pal = c("#EDF8E9", "#006D2C")))
  seaad <- read_excel(file.path(built, "SEAAD_dataset.xlsx"), sheet = "donor") %>% filter(!is.na(rs199347), adnc != "Not AD") %>%        # drops the 9 donors without AD pathology
    mutate(genotype = factor(rs199347, geno_levels), A = rs199347_A_dosage, sex = factor(sex), apoe = apoe_e4_dose,
           Braak = match(as.character(braak), c("0", "I", "II", "III", "IV", "V", "VI")) - 1, Amyloid = c(0, 1, 1, 2, 3, 3)[num(thal) + 1],
           Plaques = match(cerad, c("Absent", "Sparse", "Moderate", "Frequent")) - 1)
  indd_ad <- indd$spread %>% mutate(A = dose_of(rs199347)) %>% filter(!is.na(A), NPDx1 == "Alzheimer's disease") %>%
    mutate(genotype = factor(geno_levels[A + 1], geno_levels), age = num(AgeatDeath), sex = factor(Sex), pmi = num(PMI), apoe_c = gsub(" ", "", APOE),
           apoe = ifelse(grepl("^E[234]/E[234]$", apoe_c), lengths(regmatches(apoe_c, gregexpr("4", apoe_c))), NA),
           Braak = num(Braak06), Amyloid = num(ABeta), Plaques = num(CERAD))
  ordinal <- function(data, o, cohort) {
    d <- data %>% filter(!is.na(.data[[o]]), !is.na(age), !is.na(sex), !is.na(pmi), !is.na(apoe)) %>% mutate(y = factor(.data[[o]], ordered = TRUE))
    s <- coef(summary(polr(y ~ A + age + sex + pmi + apoe, data = d, Hess = TRUE)))["A", ]
    data.frame(cohort = cohort, outcome = o, OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]]), p = 2 * pnorm(-abs(s[[3]])), n = nrow(d))
  }
  cohorts <- list("INDD" = indd_ad, "SEA-AD" = seaad)
  stats <- do.call(rbind, lapply(names(cohorts), function(k) do.call(rbind, lapply(names(outcomes), function(o) ordinal(cohorts[[k]], o, k)))))
  stats$p_BH <- ave(stats$p, stats$cohort, FUN = function(p) p.adjust(p, "BH"))
  write.csv(stats, file.path(out8, "AD_pathology_stats.csv"), row.names = FALSE)
  theme_p <- theme_classic(base_size = 9.5) + theme(axis.text = element_text(colour = "black", size = 8.5), axis.text.x = element_text(size = 8), axis.title = element_text(size = 9.5),
        plot.title = element_text(hjust = 0.5, size = 10, lineheight = 0.95), plot.subtitle = element_text(hjust = 0.5, size = 9, lineheight = 1.05), axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35), legend.position = "bottom", legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 8.5),
        legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1.5, "mm"), plot.margin = margin(2, 3, 2, 2))
  header <- function(text) ggplot() + annotate("text", x = 0, y = 0, label = text, fontface = "bold", size = 3.2) + ylim(-1, 1) + theme_void() + theme(plot.margin = margin(1, 0, 1, 0))
  draw <- function(data, cohort, label, file) {
    panels <- list()
    for (i in seq_along(outcomes)) {
      o <- names(outcomes)[i]; lab <- outcomes[[o]]$labels
      d <- data %>% filter(!is.na(.data[[o]])); present <- lab[sort(unique(d[[o]])) + 1]
      cols <- setNames(colorRampPalette(outcomes[[o]]$pal)(length(lab)), lab)[present]
      d$stage <- factor(lab[d[[o]] + 1], rev(present)); counts <- table(d$genotype); s <- stats[stats$cohort == cohort & stats$outcome == o, ]
      panels[[o]] <- ggplot(d, aes(genotype, fill = stage)) + geom_bar(position = "fill", colour = "black", linewidth = 0.2, width = 0.63) +
        scale_fill_manual(values = cols, breaks = present) + scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
        scale_x_discrete(labels = paste0(names(counts), "\n(n = ", counts, ")")) +
        labs(x = "rs199347 genotype", y = if (i == 1) "Donors (%)" else NULL, title = outcomes[[o]]$title, subtitle = sprintf("OR %.2f per A allele
P = %s", s$OR, fmt_p(s$p))) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + theme_p
    }
    save_pdf(wrap_plots(panels, nrow = 1) & theme(plot.margin = margin(2, 5, 2, 2)), file, 178, 110)
  }
  draw(indd_ad, "INDD", "INDD: primary AD autopsies", file.path(out8, "AD_Pathology_INDD.pdf"))
  draw(seaad, "SEA-AD", "SEA-AD: donors with AD neuropathology", file.path(outS, "AD_Pathology_SEAAD_supplement.pdf"))
}

# ---- regional amyloid: brain cartoons (beta per A allele) and Thal-phase panel -------------------------------------------------------------
# Model (adapted from Cooper's in GlobalAveragePathPerRegionRework.R): score ~ A:region + region + age + sex + APOE e4 + (1 | donor), regions with >= 100 patients,
thal_phase <- c(MF = 1, SMT = 1, TS = 1, Ang = 1, OC = 1, Cing = 1, MC = 1, OFC = 1, Amyg = 2, CS = 2, DG = 2, EC = 2, CP = 3, GP = 3, MB = 3, Med = 4, SN = 4, Pons = 5, LC = 5, CB = 5)
build_regional_long <- function(wide, grp_pattern, regions, suffix, min_patients = 100) {
  apoe <- gsub(" ", "", wide$APOE)
  info <- data.frame(donor = wide$INDDID, A = dose_of(wide$rs199347), age = num(wide$AgeatDeath), sex = wide$Sex,
                     e4 = ifelse(grepl("^E[234]/E[234]$", apoe), lengths(regmatches(apoe, gregexpr("4", apoe))), NA))
  long <- do.call(rbind, lapply(regions, function(r) data.frame(info, region = r, score = score_path(wide[[paste0(r, suffix)]]), in_group = grepl(grp_pattern, wide$NPDx, fixed = TRUE))))
  long <- long[long$in_group & complete.cases(long) & long$sex %in% c("Male", "Female"), ]
  n_patients <- tapply(long$donor, long$region, function(x) length(unique(x)))
  long <- long[long$region %in% names(n_patients)[n_patients >= min_patients], ]
  long$phase <- paste("Phase", thal_phase[long$region])
  long
}

fig8_regional_amyloid <- function(indd) {
  wide <- indd$wide; regions <- names(thal_phase)
  groups <- c(LBD = "Lewy", AD = "Alzheimer")
  data_by_disease <- lapply(groups, function(p) build_regional_long(wide, p, regions, "AntibodyPlaques"))
  # betas per region and per phase
  region_stats <- data.frame(); phase_stats <- data.frame(); interaction_p <- c()
  for (g in names(groups)) {
    long <- data_by_disease[[g]]
    fit_common <- lmer(score ~ A + region + age + sex + e4 + (1 | donor), data = long)
    fit_phase  <- lmer(score ~ A:phase + region + age + sex + e4 + (1 | donor), data = long)
    fit_region <- lmer(score ~ A:region + region + age + sex + e4 + (1 | donor), data = long)
    cf <- summary(fit_phase)$coefficients; cf <- cf[grep("^A:phase", rownames(cf)), , drop = FALSE]
    phase_stats <- rbind(phase_stats, data.frame(disease = g, phase = sub("^A:phase", "", rownames(cf)), beta = cf[, "Estimate"], p = cf[, "Pr(>|t|)"], row.names = NULL))
    interaction_p[g] <- anova(fit_common, fit_phase)$`Pr(>Chisq)`[2]
    cf <- summary(fit_region)$coefficients; cf <- cf[grep("^A:region", rownames(cf)), , drop = FALSE]; reg <- sub("^A:region", "", rownames(cf))
    region_stats <- rbind(region_stats, data.frame(disease = g, region = reg, phase = paste("Phase", thal_phase[reg]), beta = cf[, "Estimate"], p = cf[, "Pr(>|t|)"],
                                                   n_patients = sapply(reg, function(r) length(unique(long$donor[long$region == r]))), row.names = NULL))
  }
  write.csv(phase_stats,  file.path(out8, "ThalPhases_phase_betas.csv"), row.names = FALSE)
  write.csv(region_stats, file.path(out8, "ThalPhases_region_betas.csv"), row.names = FALSE)

  # ---- brain cartoons: LBD (orange) and AD (teal), beta per A allele per region ----
  dk_map <- list(MF = c("rostral middle frontal", "caudal middle frontal"), SMT = c("superior temporal", "middle temporal"), TS = "banks of superior temporal sulcus",
                 Ang = "inferior parietal", OFC = c("lateral orbitofrontal", "medial orbitofrontal"), MC = "precentral",
                 OC = c("lateral occipital", "cuneus", "pericalcarine", "lingual"),
                 Cing = c("rostral anterior cingulate", "caudal anterior cingulate", "posterior cingulate", "isthmus cingulate"), EC = "entorhinal")
  aseg_map <- list(CS = "Hippocampus", DG = "Hippocampus", Amyg = "Amygdala", CP = c("Caudate", "Putamen"), GP = "Pallidum", MB = "ventraldc",
                   SN = "Brain Stem", Med = "Brain Stem", Pons = "Brain Stem", LC = "Brain Stem", CB = "Cerebellum")
  key_colours <- c("Frontal" = "#8DD3C7", "Temporal" = "#FFD92F", "Parietal (angular)" = "#BEBADA", "Occipital" = "#FB8072", "Limbic" = "#80B1D3",
                   "Brainstem" = "#B3DE69", "Cerebellum" = "#FCCDE5", "Basal ganglia / thalamus" = "#BC80BD")
  key_cortex <- c("superior frontal" = "Frontal", "rostral middle frontal" = "Frontal", "caudal middle frontal" = "Frontal", "pars opercularis" = "Frontal", "pars triangularis" = "Frontal",
                  "pars orbitalis" = "Frontal", "lateral orbitofrontal" = "Frontal", "medial orbitofrontal" = "Frontal", "frontal pole" = "Frontal", "superior temporal" = "Temporal",
                  "middle temporal" = "Temporal", "inferior temporal" = "Temporal", "banks of superior temporal sulcus" = "Temporal", "transverse temporal" = "Temporal",
                  "temporal pole" = "Temporal", "fusiform" = "Temporal", "inferior parietal" = "Parietal (angular)", "supramarginal" = "Parietal (angular)",
                  "superior parietal" = "Parietal (angular)", "precuneus" = "Parietal (angular)", "lateral occipital" = "Occipital", "lingual" = "Occipital", "cuneus" = "Occipital",
                  "pericalcarine" = "Occipital", "entorhinal" = "Limbic", "parahippocampal" = "Limbic", "caudal anterior cingulate" = "Limbic", "rostral anterior cingulate" = "Limbic",
                  "posterior cingulate" = "Limbic", "isthmus cingulate" = "Limbic")
  key_subcortex <- c("Brain Stem" = "Brainstem", "Cerebellum" = "Cerebellum", "Thalamus" = "Basal ganglia / thalamus", "Thalamus Proper" = "Basal ganglia / thalamus",
                     "Caudate" = "Basal ganglia / thalamus", "Putamen" = "Basal ganglia / thalamus", "Pallidum" = "Basal ganglia / thalamus", "Hippocampus" = "Limbic", "Amygdala" = "Limbic")
  dk_atlas <- dk(); aseg_atlas <- aseg()
  all_cortex <- unique(dk_atlas$core$region); all_subcort <- setdiff(unique(aseg_atlas$core$region), c("vessel", "choroid plexus", "Optic Chiasm"))
  cortex_values <- function(b) { d <- do.call(rbind, lapply(intersect(b$region, names(dk_map)), function(r) data.frame(region = dk_map[[r]], beta = b$beta[b$region == r])))
    merge(data.frame(region = all_cortex), d, by = "region", all.x = TRUE) }
  subcortex_values <- function(b) { d <- do.call(rbind, lapply(intersect(b$region, names(aseg_map)), function(r) data.frame(region = aseg_map[[r]], beta = b$beta[b$region == r])))
    merge(data.frame(region = all_subcort), aggregate(beta ~ region, d, mean), by = "region", all.x = TRUE) }
  brain_theme <- theme_void(base_size = 8) + theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 7),
                                                   legend.key.height = unit(3, "mm"), plot.margin = margin(0, 2, 0, 2))
  cortex_layer <- function(d, aes_fill) geom_brain(atlas = dk_atlas, data = d, hemi = "left", view = c("lateral", "medial"), position = position_brain(hemi ~ view), mapping = aes_fill, colour = "grey40", linewidth = 0.15)
  sagittal_layer <- function(d, aes_fill) geom_brain(atlas = aseg_atlas, data = d, view = "sagittal", mapping = aes_fill, colour = "grey40", linewidth = 0.15)
  kc <- data.frame(region = all_cortex, area = key_cortex[all_cortex]); ks <- data.frame(region = all_subcort, area = key_subcortex[all_subcort])
  key_scale <- scale_fill_manual(values = key_colours, na.value = "white", breaks = names(key_colours), name = NULL, drop = FALSE, guide = guide_legend(nrow = 2, byrow = TRUE))
  key_block <- wrap_plots(ggplot() + cortex_layer(kc, aes(fill = area)) + key_scale + brain_theme, ggplot() + sagittal_layer(ks, aes(fill = area)) + key_scale + brain_theme, widths = c(2.2, 1)) /
    guide_area() + plot_layout(guides = "collect", heights = c(1, 0.3))
  scales_by_disease <- list(LBD = list(colours = c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"), limits = c(0, 0.35), breaks = c(0, 0.1, 0.2, 0.3)),
                            AD  = list(colours = c("#EAF6F5", "#B4DEDB", "#5FB6B2", "#028C8C", "#015757"), limits = c(0, 0.20), breaks = c(0, 0.05, 0.10, 0.15, 0.20)))
  for (g in names(groups)) {
    b <- region_stats[region_stats$disease == g, ]; sc <- scales_by_disease[[g]]
    value_scale <- scale_fill_gradientn(colours = sc$colours, limits = sc$limits, breaks = sc$breaks, oob = scales::squish, na.value = "white",
                                        name = "Change in mean A\u03b2 plaque grade per A allele (\u03b2)", guide = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 10, barheight = 0.8))
    value_block <- wrap_plots(ggplot() + cortex_layer(cortex_values(b), aes(fill = beta)) + value_scale + brain_theme,
                              ggplot() + sagittal_layer(subcortex_values(b), aes(fill = beta)) + value_scale + brain_theme, widths = c(2.2, 1)) /
      guide_area() + plot_layout(guides = "collect", heights = c(1, 0.3))
    figure <- wrap_plots(key_block, value_block, ncol = 1) + plot_annotation(title = paste0("Neuropath ", g, ": A\u03b2 plaques"), theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5)))
    save_pdf(figure, file.path(out8, paste0("Cartoon_Abeta_beta_per_allele_", g, ".pdf")), 115, 125)
  }
  # ---- merged cartoon 
  merged_theme <- theme_void(base_size = 6.5) + theme(legend.position = "right", legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
                                                      legend.key.size = unit(2.6, "mm"), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(0, 0, 0, 0))
  row_colour <- c(LBD = "#eb6834", AD = "#028C8C")   # same as disease_colour of the Thal-phase panel
  row_label <- function(text, colour) ggplot() + annotate("text", x = 0, y = 0, label = text, fontface = "bold", size = 7.5 / .pt, angle = 90, colour = colour) + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  view_layer <- function(d, aes_fill, view) if (view == "sagittal") sagittal_layer(d, aes_fill) else
    geom_brain(atlas = dk_atlas, data = d, hemi = "left", view = view, mapping = aes_fill, colour = "grey40", linewidth = 0.15)
  views <- c("lateral", "medial", "sagittal")
  # a single view keeps the canvas of the full atlas, so crop each view to its own drawing and use its width / height for the column widths
  view_box <- lapply(setNames(views, views), function(v) { ld <- layer_data(ggplot() + view_layer(if (v == "sagittal") ks else kc, aes(fill = area), v)); list(x = range(ld$x), y = range(ld$y)) })
  aspect <- sapply(view_box, function(b) diff(b$x) / diff(b$y))
  crop <- function(v) coord_fixed(xlim = view_box[[v]]$x, ylim = view_box[[v]]$y, expand = TRUE, clip = "off")
  brain_row <- function(d_cortex, d_sub, fill_scale, label = NULL) {
    plots <- lapply(views, function(v) ggplot() + view_layer(if (v == "sagittal") d_sub else d_cortex, aes(fill = .data[[if ("area" %in% names(d_cortex)) "area" else "beta"]]), v) +
                      crop(v) + fill_scale + merged_theme + theme(legend.position = "none"))
    legend_src <- if ("area" %in% names(d_cortex)) data.frame(v = factor(names(key_colours), names(key_colours))) else data.frame(v = range(fill_scale$limits))   # every category / full range
    legend <- cowplot::get_plot_component(ggplot(legend_src, aes(0, 0, fill = v)) + geom_tile() + fill_scale + merged_theme, "guide-box-right")
    list(if (is.null(label)) plot_spacer() else row_label(label, row_colour[[label]]), plots[[1]], plot_spacer(), plots[[2]], plot_spacer(), plots[[3]], plot_spacer(), wrap_elements(full = legend))
  }
  key_scale_col <- scale_fill_manual(values = key_colours, na.value = "white", breaks = names(key_colours), name = NULL, drop = FALSE, guide = guide_legend(ncol = 1))
  beta_scale <- function(g) { sc <- scales_by_disease[[g]]
    scale_fill_gradientn(colours = sc$colours, limits = sc$limits, breaks = sc$breaks, oob = scales::squish, na.value = "white", name = "β per
A allele",
                         guide = guide_colourbar(title.position = "top", barwidth = unit(2.2, "mm"), barheight = unit(16, "mm"))) }
  gap <- 0.08
  cells <- c(brain_row(kc, ks, key_scale_col),
             brain_row(cortex_values(region_stats[region_stats$disease == "AD", ]), subcortex_values(region_stats[region_stats$disease == "AD", ]), beta_scale("AD"), "AD"),
             brain_row(cortex_values(region_stats[region_stats$disease == "LBD", ]), subcortex_values(region_stats[region_stats$disease == "LBD", ]), beta_scale("LBD"), "LBD"))
  merged <- wrap_plots(cells, ncol = 8, widths = c(0.12, aspect[["lateral"]], gap, aspect[["medial"]], gap, aspect[["sagittal"]], gap, 1.2), heights = c(1, 1, 1))
  save_pdf(merged, file.path(out8, "Cartoon_Abeta_beta_per_allele_merged.pdf"), 138, 80)

  # ---- Thal-phase panel 
  roman <- c("I", "II", "III", "IV", "V"); grades <- c("0", "Rare", "1+", "2+", "3+")
  disease_colour <- c(LBD = "#eb6834", AD = "#028C8C")
  heat_ramp <- list(LBD = c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"), AD = c("#EAF6F5", "#B4DEDB", "#5FB6B2", "#028C8C", "#015757"))
  heat_limit <- list(LBD = c(0, 0.35), AD = c(0, 0.20))
  for (g in names(data_by_disease)) {
    x <- data_by_disease[[g]]; x$genotype <- factor(geno_levels[x$A + 1], geno_levels); x$grade <- factor(grades[match(x$score, c(0, 0.5, 1, 2, 3))], grades); data_by_disease[[g]] <- x
  }
  theme_h <- theme_classic(base_size = 5.5) + theme(axis.text = element_text(colour = "black", size = 5), axis.title = element_text(size = 5.5), axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25), axis.ticks.length = unit(0.6, "mm"), plot.title = element_text(face = "bold", hjust = 0.5, size = 5.5, margin = margin(0, 0, 0.6, 0)),
        plot.margin = margin(1, 1, 1, 1), legend.title = element_blank(), legend.text = element_text(size = 5), legend.key.size = unit(2.2, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"), legend.spacing.y = unit(0.4, "mm"))
  stack_panel <- function(g, phase, show_y) {
    d <- data_by_disease[[g]] %>% filter(phase == !!phase)
    fills <- setNames(colorRampPalette(c("#f3f2ef", disease_colour[[g]]))(5), grades)
    shares <- d %>% count(genotype, grade, .drop = FALSE) %>% group_by(genotype) %>% mutate(share = n / sum(n)) %>% ungroup() %>% filter(n > 0)
    st <- phase_stats[phase_stats$disease == g & phase_stats$phase == phase, ]
    ggplot(shares, aes(genotype, share, fill = grade)) + geom_col(position = "stack", width = 0.72, colour = "white", linewidth = 0.15) +
      scale_fill_manual(values = fills, breaks = rev(grades), drop = FALSE) +
      scale_y_continuous(labels = c("0", "50", "100"), breaks = c(0, 0.5, 1), expand = c(0, 0)) + coord_cartesian(ylim = c(0, 1), clip = "off") +
      labs(x = NULL, y = if (show_y) "Regional scores (%)" else NULL) +
      labs(title = paste0("Phase ", roman[as.integer(sub("Phase ", "", phase))]),
           subtitle = paste0("β ", sprintf("%.2f", st$beta), ", P ", sub("e-0", "e-", fmt_p(st$p)))) + theme_h +
      theme(plot.subtitle = element_text(hjust = 0.5, size = 4.6, margin = margin(0, 0, 1.2, 0))) +
      theme(axis.text.y = if (show_y) element_text(size = 5) else element_blank(), axis.ticks.y = if (show_y) element_line(linewidth = 0.25) else element_blank())
  }
  label_row <- function(g) wrap_elements(full = (ggplot() + annotate("text", x = 0, y = 0, label = g, fontface = "bold", size = 6.5 / .pt, angle = 90, colour = disease_colour[[g]]) + theme_void()))
  rows <- lapply(c("AD", "LBD"), function(g) wrap_plots(lapply(1:5, function(i) stack_panel(g, paste("Phase", i), i == 1)), nrow = 1, guides = "collect") & theme(legend.position = "right"))
  v1 <- wrap_plots(label_row("AD"), rows[[1]], label_row("LBD"), rows[[2]], ncol = 2, widths = c(0.035, 1))   # AD on top
  message(sprintf("Thal phases: genotype x phase interaction LBD P = %s, AD P = %s", fmt_p(interaction_p[["LBD"]]), fmt_p(interaction_p[["AD"]])))
  save_pdf(v1, file.path(out8, "ThalPhases_stacked_all_phases.pdf"), 300 * 25.4 / 72, 155 * 25.4 / 72)
}

# ---- three-disease forest plot (alpha-synuclein: brainstem regions only) --------------------------------------------------
# one mixed model per disease x pathology: z ~ A + age + sex + APOE e4 + region + (1 | donor), z = regional score z-scored inside each region
fig8_forest <- function(indd) {
  wide <- indd$wide; apoe <- gsub(" ", "", wide$APOE)
  donor_info <- data.frame(donor = wide$INDDID, A = dose_of(wide$rs199347), age = num(wide$AgeatDeath), sex = wide$Sex,
                           e4 = ifelse(grepl("^E[234]/E[234]$", apoe), lengths(regmatches(apoe, gregexpr("4", apoe))), NA))
  suffix <- c(Abeta = "AntibodyPlaques", Tau = "Tau", aSyn = "aSyn", TDP43 = "TDP43")
  pathology_label <- c(Abeta = "A\u03b2 plaques", Tau = "Tau", aSyn = "\u03b1-synuclein", TDP43 = "pTDP-43")
  pathology_colour <- setNames(c("#eb6834", "#2a78d6", "#4a3aa7", "#2e8b57"), pathology_label)
  diseases <- c("Neuropath AD" = "Alzheimer", "Neuropath LBD" = "Lewy", "Neuropath ALS" = "Amyotrophic")
  abeta_regions <- c("MF", "SMT", "TS", "Ang", "OC", "Cing"); tau_regions <- c("EC", "CS", "Amyg"); brainstem <- c("Med", "LC", "SN"); motor_brainstem <- c("Pons", "Med", "MC", "SC")
  make_sets <- function(asyn) list("Neuropath AD" = list(Abeta = abeta_regions, Tau = tau_regions, aSyn = asyn, TDP43 = NULL),
                                   "Neuropath LBD" = list(Abeta = abeta_regions, Tau = tau_regions, aSyn = asyn, TDP43 = NULL),
                                   "Neuropath ALS" = list(Abeta = abeta_regions, Tau = tau_regions, aSyn = NULL, TDP43 = motor_brainstem))   # NULL = not tested (pathology absent)
  forest_row <- function(disease, p, regs) {
    if (is.null(regs)) return(NULL)
    long <- do.call(rbind, lapply(regs, function(r) data.frame(donor_info, region = r, score = score_path(wide[[paste0(r, suffix[[p]])]]), in_group = grepl(diseases[[disease]], wide$NPDx, fixed = TRUE))))
    long <- long[long$in_group & complete.cases(long) & long$sex %in% c("Male", "Female"), ]
    long$z <- ave(long$score, long$region, FUN = function(x) (x - mean(x)) / sd(x)); long <- long[is.finite(long$z), ]
    if (length(unique(long$donor)) < 30) return(NULL)
    f <- if (length(unique(long$region)) > 1) z ~ A + age + sex + e4 + region + (1 | donor) else z ~ A + age + sex + e4 + (1 | donor)
    s <- summary(lmer(f, data = long))$coefficients["A", ]
    data.frame(disease = disease, pathology = p, beta = s[["Estimate"]], lo = s[["Estimate"]] - 1.96 * s[["Std. Error"]], hi = s[["Estimate"]] + 1.96 * s[["Std. Error"]],
               p = s[["Pr(>|t|)"]], n = length(unique(long$donor)), regions = paste(unique(long$region), collapse = ","))
  }
  run_forest <- function(sets) { f <- do.call(rbind, lapply(names(diseases), function(d) do.call(rbind, lapply(names(suffix), function(p) forest_row(d, p, sets[[d]][[p]])))))
    f$p_BH <- ave(f$p, f$disease, FUN = function(p) p.adjust(p, "BH")); f }      # BH within each disease (csv only; the plot shows nominal P)
  draw_forest <- function(forest, sets, file) {
    plot_data <- forest %>% mutate(pathology = factor(pathology_label[pathology], rev(unname(pathology_label))), disease = factor(disease, names(diseases)),
                                   p_label = paste0("P = ", signif(p, 2)), p_face = ifelse(p < 0.05, "bold", "plain"), n_label = paste0("n = ", n))
    skipped <- do.call(rbind, lapply(names(diseases), function(d) { m <- names(suffix)[sapply(names(suffix), function(p) is.null(sets[[d]][[p]]))]
      if (length(m)) data.frame(disease = factor(d, names(diseases)), pathology = factor(pathology_label[m], rev(unname(pathology_label)))) }))
    figure <- ggplot(plot_data, aes(beta, pathology, colour = pathology)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.25) +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", linewidth = 0.4) + geom_point(size = 1.6) +
      geom_text(aes(x = 0.45, label = p_label, fontface = p_face), hjust = 0, nudge_y = 0.15, size = 2.1, colour = "black") +
      geom_text(aes(x = 0.45, label = n_label), hjust = 0, nudge_y = -0.2, size = 1.8, colour = "grey35") +
      facet_wrap(~ disease, nrow = 1) + scale_colour_manual(values = pathology_colour, guide = "none") + scale_y_discrete(drop = FALSE) +
      scale_x_continuous(limits = c(-0.3, 0.8), breaks = c(-0.2, 0, 0.2, 0.4)) + labs(x = "SD per A allele (95% CI)", y = NULL) + panel_theme(8) + theme(panel.spacing = unit(4, "mm"))
    if (!is.null(skipped)) figure <- figure + geom_text(data = skipped, aes(x = 0.1, y = pathology, label = "not tested"), inherit.aes = FALSE, colour = "grey45", size = 2, fontface = "italic")
    save_pdf(figure, file, 185, 62)
  }
  fa <- run_forest(make_sets(brainstem)); write.csv(fa, file.path(out8, "Forest_stats_aSyn_brainstem.csv"), row.names = FALSE)
  draw_forest(fa, make_sets(brainstem), file.path(out8, "Forest_aSyn_brainstem.pdf"))
}

make_figure8 <- function() {
  if (has[["fig8_microglia"]]) fig8_expression()
  if (has[["fig8_indd"]]) {
    indd <- load_indd()
    if (has[["fig8_microglia"]]) fig8_ad_pathology(indd)   # also needs SEAAD_dataset.xlsx (supplement)
    fig8_regional_amyloid(indd)
    fig8_forest(indd)
  }
  message("Figure 8 done -> ", out8, " and ", outS)
}

if (RUN_FIG7 && has[["fig7"]]) make_figure7()
if (RUN_FIG8) make_figure8()
