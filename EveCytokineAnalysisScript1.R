# =============================================================================
#  Eve Technologies cytokine panel: GPNMB WT / HET / KO microglia
# =============================================================================
#
#  What this script does
#  ---------------------
#  Step 1   Mean concentration of every analyte, to spot analytes that sit
#           unusually high alongside the media components (M-CSF, IL-34).
#  Step 1b  Table of analytes with unreadable (NA) wells.
#  Step 2   Technical-replicate variability; analytes whose a/b wells differ
#           by > 50% of their mean are removed from Steps 3-4.
#  Step 3   KO vs WT differential analysis with limma (moderated t-test).
#  Step 4   Significant protein hits plotted next to KO vs WT bulk RNA-seq.
#
#  Inputs
#  ------
#  eve_file  Eve export (.xlsx), Obs Conc. columns in pg/mL
#  rna_file  DESeq2 results, columns: X (Ensembl), log2FoldChange, padj,
#            baseMean, gene_name. NOTE: this contrast is WT vs KO (GPNMB is
#            positive), so Step 4 flips its sign to match the protein data.
#
#  Outputs (all written to out_dir)
#  --------------------------------
#  Step1_analyte_means.pdf
#  Step1b_missing_wells.csv
#  Step2_technical_variability.pdf / .csv
#  Step3_KO_vs_WT.pdf / Step3_KO_vs_WT_limma.csv
#  Step4_sig_hits_vs_RNA.pdf / .csv
#  Step4_all_analytes_vs_RNA.csv
#
#  Sample layout
#  -------------
#  3 genotypes x 2 biological replicates (1, 2) x 2 technical wells (a, b)
#  = 12 wells. Samples 13-14 are pilots and are excluded from all analyses.
# =============================================================================


# =============================================================================
# Setup: packages, file paths, output folder
# =============================================================================
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

eve_file <- "/Volumes/PC60/EVECytokines_edited.xlsx"
rna_file <- "/Volumes/PC60/bulkRNASeq/RNaSeqResultsTables/MicrogliaKOData/BackgroundDifferences.csv"
out_dir  <- "/Volumes/PC60/EveCytokineAnalyses"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# Base pdf() can't encode Greek letters (α, β, γ) and silently drops them.
# cairo_pdf embeds Unicode fonts, so every figure is saved through this helper.
# Base pdf() drops Greek letters (α, β, γ); macOS's Quartz PDF device keeps them
save_pdf <- function(file, plot, width, height) {
  ggsave(file.path(out_dir, file), plot, width = width, height = height,
         device = function(filename, width, height, ...)
           grDevices::quartz(file = filename, width = width, height = height, type = "pdf"))
}
# Analytes flagged as possibly above the upper limit of quantitation
over_range <- c("MCP-1", "IL-8", "IL-34", "MCP-4")

# Analysis thresholds (change here, not further down)
tech_cut <- 30     # Step 2: max mean |a - b| as % of the analyte's mean
sig_cut  <- 0.05   # Steps 3-4: BH-adjusted p (q) threshold for protein hits

# Eve sample number -> sample name
key <- data.frame(
  sample_num = 1:14,
  sample = c("WT PBS 1a",  "WT PBS 1b",  "WT PBS 2a",  "WT PBS 2b",
             "HET PBS 1a", "HET PBS 1b", "HET PBS 2a", "HET PBS 2b",
             "KO PBS 1a",  "KO PBS 1b",  "KO PBS 2a",  "KO PBS 2b",
             "Pilot 1 (WT)", "Pilot 2 (KO)")
)


# =============================================================================
# Load the Eve export
# =============================================================================
# The sheet has stacked header rows, so it is read with no column names and
# every cell as text (keeps "---" / "OOR <" intact until we convert on purpose).
#
# Layout after dropping empty columns:
#   row 1 = units, row 2 = analyte names, row 3 = FI / Obs Conc. / Exp Conc.
#   then standards, then one row per sample (first column = sample number).
# Each analyte occupies 3 columns: FI | Obs Conc. | Exp Conc.
# The analyte name sits above the FI column, one column LEFT of Obs Conc.

raw <- read_excel(eve_file, col_names = FALSE, col_types = "text")
raw <- raw[, colSums(!is.na(raw)) > 0]           # Eve leaves column A blank

analyte_row <- unlist(raw[2, ])                   # EGF, Eotaxin, ...
metric_row  <- unlist(raw[3, ])                   # FI, Obs Conc., Exp Conc., ...
samples     <- raw[grepl("^[0-9]+$", raw[[1]]), ] # sample rows only (drops standards)

# Keep the Obs Conc. column of each analyte; non-numeric cells become NA
obs_cols <- which(metric_row == "Obs Conc.")
conc <- as.data.frame(lapply(samples[, obs_cols], as.numeric))
names(conc) <- analyte_row[obs_cols - 1]

# Attach sample names
conc <- data.frame(
  sample = key$sample[match(as.integer(samples[[1]]), key$sample_num)],
  conc,
  check.names = FALSE                             # keep names like "IL-1β" as-is
)

# Long format used by every step: one row per well x analyte, pilots removed
wells <- conc %>%
  filter(!grepl("Pilot", sample)) %>%
  pivot_longer(-sample, names_to = "analyte", values_to = "conc") %>%
  mutate(genotype = sub(" .*", "", sample),              # WT / HET / KO
         bio_rep  = sub(".* (\\d)[ab]$", "\\1", sample)) # biological replicate 1 / 2
# =============================================================================
# Step 1: Mean level of each analyte (10 highest and 10 lowest)
# =============================================================================
# One point per analyte = mean of the 12 wells; grey line = min to max.
# Orange = media components (M-CSF, IL-34)
# Blue   = possibly above the upper limit of detection (set in `over_range`)
# Analytes with no readable wells are excluded here (see Step 1b).

avg <- wells %>%
  group_by(analyte) %>%
  summarise(mean = mean(conc, na.rm = TRUE),
            min  = suppressWarnings(min(conc, na.rm = TRUE)),
            max  = suppressWarnings(max(conc, na.rm = TRUE)),
            n    = sum(!is.na(conc)),
            .groups = "drop") %>%
  filter(n > 0) %>%
  mutate(category = case_when(analyte %in% c("M-CSF", "IL-34") ~ "Media component",
                              analyte %in% over_range          ~ "Above detection range",
                              TRUE                             ~ "Other"))

top_bottom <- bind_rows(
  avg %>% slice_max(mean, n = 20) %>% mutate(group = "Highest 20"),
  avg %>% slice_min(mean, n = 20) %>% mutate(group = "Lowest 20")
) %>%
  mutate(group   = factor(group, levels = c("Highest 20", "Lowest 20")),
         analyte = reorder(analyte, mean))

p1 <- ggplot(top_bottom, aes(x = mean, y = analyte)) +
  geom_segment(aes(x = min, xend = max, yend = analyte), colour = "grey75", linewidth = 0.4) +
  geom_point(aes(colour = category), size = 2.2) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_colour_manual(values = c("Media component"       = "#eb6834",
                                 "Above detection range" = "#2a78d6",
                                 "Other"                 = "grey30"),
                      name = NULL) +
  labs(x = "Concentration, pg/mL (log scale)\npoint = mean of 12 wells, line = min to max",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"),
        legend.position = "top")

p1
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step1_analyte_means_top_bottom10.pdf", p1,
         width = 6, height = 6)
write.csv(top_bottom, file.path(out_dir, "Step1_top_bottom10.csv"), row.names = FALSE)




# =============================================================================
# Step 1b: Analytes with unreadable wells
# =============================================================================
# Lists every analyte with at least one NA well, which wells, and what Eve
# actually wrote in those cells. "OOR <", "OOR >" and "---" are expected.
# NOTE: Eve marks extrapolated values with "*" (e.g. "*1.23"). as.numeric()
# turns those into NA as well; if any appear in raw_text, they are real
# (less certain) measurements and could be recovered by stripping the "*".

raw_obs <- data.frame(sample = conc$sample, samples[, obs_cols], check.names = FALSE)
names(raw_obs)[-1] <- names(conc)[-1]

na_table <- wells %>%
  left_join(raw_obs %>% pivot_longer(-sample, names_to = "analyte", values_to = "raw_text"),
            by = c("sample", "analyte")) %>%
  group_by(analyte) %>%
  summarise(n_readable    = sum(!is.na(conc)),
            n_missing     = sum(is.na(conc)),
            missing_wells = if (all(is.na(conc))) "all" else paste(sample[is.na(conc)], collapse = ", "),
            raw_text      = paste(unique(raw_text[is.na(conc)]), collapse = "; "),
            .groups = "drop") %>%
  filter(n_missing > 0) %>%
  mutate(status = ifelse(n_readable == 0, "Dropped (no readable wells)", "Kept, partly missing")) %>%
  arrange(desc(n_missing))

na_table
write.csv(na_table, file.path(out_dir, "Step1b_missing_wells.csv"), row.names = FALSE)


# =============================================================================
# Step 2: Technical-replicate variability (intra-assay CV) and filter
# =============================================================================
# For each analyte:
#   1. CV for each of the 6 duplicate-well pairs (WT 1, WT 2, HET 1, HET 2, KO 1, KO 2):
#        CV = SD(a, b) / mean(a, b) x 100
#   2. average the pair CVs                               -> mean_cv (%)
# Analytes with mean_cv > tech_cut (set in Setup) are excluded from Steps 3-4.
# Each pair is scaled by its own mean, so genotype differences in level don't
# inflate the CV. The filter uses ALL analytes; the figure shows the 20 most
# and 20 least variable.

tech_cut   # CV cutoff (%) from Setup

tech_var <- wells %>%
  group_by(analyte, genotype, bio_rep) %>%
  summarise(cv = 100 * sd(conc) / mean(conc), .groups = "drop") %>%   # NA if a well is missing
  filter(is.finite(cv)) %>%
  group_by(analyte) %>%
  summarise(mean_cv = mean(cv), n_pairs = n(), .groups = "drop") %>%  # n_pairs = pairs used (max 6)
  mutate(status = ifelse(mean_cv > tech_cut,
                         paste0("Removed (CV > ", tech_cut, "%)"), "Kept"))

keep_analytes <- tech_var$analyte[tech_var$mean_cv <= tech_cut]
removed       <- tech_var$analyte[tech_var$mean_cv >  tech_cut]
removed                                                           # analytes excluded here
write.csv(tech_var, file.path(out_dir, "Step2_technical_variability_CV.csv"), row.names = FALSE)

# ---- Figure: 20 most and 20 least variable analytes -------------------------
var_top_bottom <- bind_rows(
  tech_var %>% slice_max(mean_cv, n = 20) %>% mutate(group = "Most variable (top 20)"),
  tech_var %>% slice_min(mean_cv, n = 20) %>% mutate(group = "Least variable (bottom 20)")
) %>%
  mutate(group   = factor(group, levels = c("Most variable (top 20)", "Least variable (bottom 20)")),
         analyte = reorder(analyte, mean_cv))

p2 <- ggplot(var_top_bottom, aes(x = mean_cv, y = analyte, colour = status)) +
  geom_vline(xintercept = tech_cut, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(size = 2) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = setNames(c("grey30", "#eb6834"),
                                        c("Kept", paste0("Removed (CV > ", tech_cut, "%)"))),
                      name = NULL) +
  labs(x = "Mean intra-assay CV of duplicate wells (%)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"),
        legend.position = "top")

p2
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step2_technical_variability_CV_top_bottom20.pdf", p2,
         width = 6.5, height = 8)



# =============================================================================
# Step 2: Technical-replicate variability (intra-assay CV) and filter
# =============================================================================
# For each analyte:
#   1. CV for each of the 6 duplicate-well pairs (WT 1, WT 2, HET 1, HET 2, KO 1, KO 2):
#        CV = SD(a, b) / mean(a, b) x 100
#   2. average the pair CVs                               -> mean_cv (%)
# Analytes with mean_cv > tech_cut (set in Setup) are excluded from Steps 3-4.
# Each pair is scaled by its own mean, so genotype differences in level don't
# inflate the CV. The filter uses ALL analytes; the figure shows the 20 most
# and 20 least variable.

tech_cut   # CV cutoff (%) from Setup

tech_var <- wells %>%
  group_by(analyte, genotype, bio_rep) %>%
  summarise(cv = 100 * sd(conc) / mean(conc), .groups = "drop") %>%   # NA if a well is missing
  filter(is.finite(cv)) %>%
  group_by(analyte) %>%
  summarise(mean_cv = mean(cv), n_pairs = n(), .groups = "drop") %>%  # n_pairs = pairs used (max 6)
  mutate(status = ifelse(mean_cv > tech_cut,
                         paste0("Removed (CV > ", tech_cut, "%)"), "Kept"))

keep_analytes <- tech_var$analyte[tech_var$mean_cv <= tech_cut]
removed       <- tech_var$analyte[tech_var$mean_cv >  tech_cut]
removed                                                           # analytes excluded here
write.csv(tech_var, file.path(out_dir, "Step2_technical_variability_CV.csv"), row.names = FALSE)

# ---- Figure: 20 most and 20 least variable analytes -------------------------
var_top_bottom <- bind_rows(
  tech_var %>% slice_max(mean_cv, n = 20) %>% mutate(group = "Most variable (top 20)"),
  tech_var %>% slice_min(mean_cv, n = 20) %>% mutate(group = "Least variable (bottom 20)")
) %>%
  mutate(group   = factor(group, levels = c("Most variable (top 20)", "Least variable (bottom 20)")),
         analyte = reorder(analyte, mean_cv))

p2 <- ggplot(var_top_bottom, aes(x = mean_cv, y = analyte, colour = status)) +
  geom_vline(xintercept = tech_cut, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(size = 2) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = setNames(c("grey30", "#eb6834"),
                                        c("Kept", paste0("Removed (CV > ", tech_cut, "%)"))),
                      name = NULL) +
  labs(x = "Mean intra-assay CV of duplicate wells (%)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"),
        legend.position = "top")

p2
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step2_technical_variability_CV_top_bottom20.pdf", p2,
         width = 6.5, height = 8)



# =============================================================================
# Step 3: KO vs WT with limma (well-level, n = 4 wells per genotype)
# =============================================================================
# Each of the 12 wells is treated as an independent data point
# (2 cultures x 2 technical wells per genotype).
# limma fits a linear model per analyte, then shrinks each variance toward the
# value typical for the panel (empirical Bayes); trend = TRUE lets that typical
# variance depend on abundance. Tests are moderated t-tests with BH correction.
# Refs: Smyth 2004, Stat Appl Genet Mol Biol, doi:10.2202/1544-6115.1027
#       Ritchie et al. 2015, Nucleic Acids Res, doi:10.1093/nar/gkv007
#
# logFC = log2(KO) - log2(WT): positive = higher in KO.
#
# Outputs: Step3_KO_vs_WT_limma.csv          (all tested analytes)
#          Step3_KO_vs_WT.pdf                (all tested analytes)
#          Step3_KO_vs_WT_significant.pdf    (q < sig_cut only)

# ---- 3a. Choose the analytes to test ----------------------------------------
# Need a readable well in every genotype, otherwise KO vs WT can't be estimated
# (e.g. 6CKine: readable only in WT 1). Media components are excluded.
has_all_geno <- wells %>%
  group_by(analyte, genotype) %>%
  summarise(any_data = any(!is.na(conc)), .groups = "drop") %>%
  group_by(analyte) %>%
  summarise(ok = all(any_data), .groups = "drop") %>%
  filter(ok) %>%
  pull(analyte)

setdiff(keep_analytes, has_all_geno)                  # removed here (expect 6CKine)
test_analytes <- setdiff(intersect(keep_analytes, has_all_geno), c("M-CSF", "IL-34"))
length(test_analytes)                                 # expect 78

# ---- 3b. Well-level log2 values ---------------------------------------------
well_dat <- wells %>%
  filter(analyte %in% test_analytes) %>%
  mutate(lconc = log2(conc + 1))                      # +1 pg/mL avoids log(0)

# Analyte x well matrix (12 columns: "WT PBS 1a", "WT PBS 1b", ...)
mat <- well_dat %>%
  select(analyte, sample, lconc) %>%
  pivot_wider(names_from = sample, values_from = lconc) %>%
  tibble::column_to_rownames("analyte") %>%
  as.matrix()

geno <- factor(sub(" .*", "", colnames(mat)), levels = c("WT", "HET", "KO"))
table(geno)                                           # expect 4 wells per genotype

# ---- 3c. Fit the model ------------------------------------------------------
design <- model.matrix(~ 0 + geno)
colnames(design) <- levels(geno)

contr <- makeContrasts(KO_vs_WT = KO - WT, HET_vs_WT = HET - WT, levels = design)
fit   <- eBayes(contrasts.fit(lmFit(mat, design), contr), trend = TRUE)

fit$df.prior                                          # degrees of freedom borrowed from the panel

# ---- 3d. Results table ------------------------------------------------------
res <- topTable(fit, coef = "KO_vs_WT", number = Inf, sort.by = "none", confint = TRUE) %>%
  tibble::rownames_to_column("analyte") %>%
  filter(!is.na(logFC)) %>%
  mutate(sig = ifelse(adj.P.Val < sig_cut, paste("q <", sig_cut), paste("q ≥", sig_cut)))

write.csv(res, file.path(out_dir, "Step3_KO_vs_WT_limma.csv"), row.names = FALSE)

# ---- 3e. Fold-change figures ------------------------------------------------
# Point = log2FC, grey bar = 95% CI, colour = p-value,
# filled = q < sig_cut, open = not significant
plot_fc <- function(df) {
  df <- df %>% mutate(analyte = reorder(analyte, logFC))
  ggplot(df, aes(x = logFC, y = analyte)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_linerange(aes(xmin = CI.L, xmax = CI.R), colour = "grey75", linewidth = 0.4) +
    geom_point(aes(colour = -log10(P.Value), shape = sig), size = 2) +
    scale_colour_gradient(low = "#86b6ef", high = "#104281", name = "-log10(p)") +
    scale_shape_manual(values = setNames(c(16, 1), c(paste("q <", sig_cut), paste("q ≥", sig_cut))),
                       name = NULL) +
    labs(x = "log2 fold change, KO vs WT (95% CI; positive = higher in KO)", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
}

# All tested analytes
p3 <- plot_fc(res)
p3
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step3_KO_vs_WT.pdf", p3,
         width = 6, height = 0.13 * nrow(res) + 1.5)

# Significant analytes only
res_sig <- res %>% filter(adj.P.Val < sig_cut)
nrow(res_sig)                                         # number of significant analytes

p3_sig <- plot_fc(res_sig)
p3_sig
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step3_KO_vs_WT_significant.pdf", p3_sig,
         width = 6, height = 0.2 * nrow(res_sig) + 1.5)


# =============================================================================
# Step 3b: Top 6 up and top 6 down in KO, with WT / HET / KO levels
# =============================================================================
# Ranked by KO vs WT log2FC among analytes passing q < sig_cut. Analytes in
# `force_down` are always shown in the "Down in KO" row (marked with *); the
# remaining slots are filled by the strongest significant decreases.
# HET vs WT log2FC and q come from the same limma fit as Step 3.
library(patchwork)

n_show     <- 6
force_down <- "TNFα"
grep("TNF", res$analyte, value = TRUE)                # check the exact spelling

het <- topTable(fit, coef = "HET_vs_WT", number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("analyte") %>%
  select(analyte, HET_log2FC = logFC, HET_q = adj.P.Val)

sig_res <- res %>% filter(adj.P.Val < sig_cut)

up_hits <- sig_res %>%
  filter(logFC > 0) %>%
  slice_max(logFC, n = n_show)

down_hits <- bind_rows(
  res %>% filter(analyte %in% force_down),                            # added by hand
  sig_res %>% filter(logFC < 0, !analyte %in% force_down) %>%
    slice_min(logFC, n = n_show - length(force_down))
) %>%
  arrange(logFC)                                                      # strongest decrease first

top_hits <- bind_rows(
  up_hits   %>% mutate(direction = "Up in KO"),
  down_hits %>% mutate(direction = "Down in KO")
) %>%
  mutate(added_manually = analyte %in% force_down,
         label          = ifelse(added_manually, paste0(analyte, "*"), analyte))

# Mean concentration per genotype (mean of the 4 wells)
geno_levels <- well_dat %>%
  filter(analyte %in% top_hits$analyte) %>%
  group_by(analyte, genotype) %>%
  summarise(mean_pg_mL = mean(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = mean_pg_mL)

top_table <- top_hits %>%
  select(direction, analyte, KO_log2FC = logFC, KO_q = adj.P.Val, added_manually) %>%
  left_join(het, by = "analyte") %>%
  left_join(geno_levels, by = "analyte") %>%
  select(direction, analyte, WT, HET, KO, KO_log2FC, KO_q, HET_log2FC, HET_q, added_manually) %>%
  mutate(across(c(WT, HET, KO), ~ round(.x, 1)),
         over_range = analyte %in% over_range)        # flagged in Setup

top_table
write.csv(top_table, file.path(out_dir, "Step3b_top_up_down.csv"), row.names = FALSE)

# ---- Figure: one row of 6 panels per direction, stacked into a 2 x 6 grid ----
# Bars = genotype mean, points = the 4 individual wells.
plot_row <- function(dir) {
  hits_dir <- top_hits %>% filter(direction == dir)
  df <- well_dat %>%
    filter(analyte %in% hits_dir$analyte, !is.na(conc)) %>%
    left_join(hits_dir %>% select(analyte, label), by = "analyte") %>%
    mutate(genotype = factor(genotype, levels = c("WT", "HET", "KO")),
           label    = factor(label, levels = hits_dir$label))   # strongest effect on the left
  
  ggplot(df, aes(x = genotype, y = conc, fill = genotype)) +
    stat_summary(fun = mean, geom = "col", width = 0.65) +
    geom_point(size = 1.6, colour = "grey15",
               position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    facet_wrap(~ label, scales = "free_y", nrow = 1) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","),
                       expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = c(WT = "grey70", HET = "#86b6ef", KO = "#2a78d6")) +
    labs(title = dir, x = NULL, y = "pg/mL") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(1, "lines"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
}

p3b <- plot_row("Up in KO") / plot_row("Down in KO")
p3b
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step3b_top_up_down.pdf", p3b, width = 12, height = 6)
# =============================================================================
# Step 4: Protein vs bulk RNA-seq, all analytes with an RNA match
# =============================================================================
# Every analyte tested in Step 3 whose gene is found in the RNA results is
# shown, significant or not. Both measures are on the KO vs WT scale.
#
# Outputs: Step4_all_analytes_vs_RNA.csv   (all tested analytes, matched or not)
#          Step4_RNA_matched.csv           (matched analytes only, as plotted)
#          Step4_RNA_matched.pdf

# ---- 4a. Load RNA results and put them on the KO vs WT scale ----------------
# GPNMB is higher in WT, so on a KO-vs-WT scale its log2FC must be negative.
# If the file reports it as positive, the file is WT vs KO and is flipped.
RNAResults <- read.csv(rna_file)

gpnmb_lfc <- RNAResults$log2FoldChange[RNAResults$gene_name == "GPNMB"][1]
rna_sign  <- ifelse(gpnmb_lfc > 0, -1, 1)
gpnmb_lfc; rna_sign                       # expect 2.93 and -1

rna <- RNAResults %>%
  transmute(gene       = toupper(trimws(gene_name)),
            RNA_log2FC = rna_sign * log2FoldChange,      # now KO vs WT, same as protein
            RNA_padj   = padj,                           # NA = filtered by DESeq2 (low counts)
            baseMean) %>%
  filter(!is.na(gene), gene != "") %>%
  group_by(gene) %>%
  slice_max(baseMean, n = 1, with_ties = FALSE) %>%     # one row per gene (alt-locus duplicates)
  ungroup() %>%
  select(-baseMean)

# ---- 4b. Map Eve protein names to gene symbols ------------------------------
# Lookup for names that don't convert directly; everything else is converted
# automatically (Greek letters -> A/B/G/D, hyphens removed: IL-1β -> IL1B).
alias <- c(
  "Eotaxin" = "CCL11", "Eotaxin-2" = "CCL24", "Eotaxin-3" = "CCL26",
  "FLT-3L" = "FLT3LG", "Fractalkine" = "CX3CL1",
  "G-CSF" = "CSF3", "GM-CSF" = "CSF2", "M-CSF" = "CSF1",
  "GROα" = "CXCL1", "IL-1RA" = "IL1RN", "IL-8" = "CXCL8",
  "IL-12p40" = "IL12B", "IL-12p70" = "IL12A", "IL-17E/IL-25" = "IL25",
  "IL-23" = "IL23A", "IL-28A" = "IFNL2", "IL-29" = "IFNL1", "IL-35" = "EBI3",
  "IP-10" = "CXCL10", "MCP-1" = "CCL2", "MCP-2" = "CCL8", "MCP-3" = "CCL7",
  "MCP-4" = "CCL13", "MDC" = "CCL22", "MIG" = "CXCL9",
  "MIP-1α" = "CCL3", "MIP-1β" = "CCL4", "MIP-1δ" = "CCL15",
  "MIP-3α" = "CCL20", "MIP-3β" = "CCL19", "RANTES" = "CCL5",
  "PDGF-AA" = "PDGFA", "PDGF-AB/BB" = "PDGFB", "TGFα" = "TGFA",
  "TNFα" = "TNF", "TNFβ" = "LTA", "sCD40L" = "CD40LG", "TARC" = "CCL17",
  "VEGF-A" = "VEGFA", "6CKine" = "CCL21", "BCA-1" = "CXCL13",
  "CTACK" = "CCL27", "ENA-78" = "CXCL5", "I-309" = "CCL1", "I-TAC" = "CXCL11",
  "SCF" = "KITLG", "SDF-1" = "CXCL12", "TPO" = "THPO", "TRAIL" = "TNFSF10"
)
auto_symbol <- function(x) toupper(gsub("[- ]", "", chartr("αβγδ", "ABGD", x)))
to_gene     <- function(x) ifelse(x %in% names(alias), unname(alias[x]), auto_symbol(x))

# ---- 4c. Reference table: every tested analyte with its RNA result ----------
all_vs_rna <- res %>%
  mutate(analyte = as.character(analyte),
         gene    = to_gene(analyte)) %>%
  select(analyte, gene, prot_log2FC = logFC, CI.L, CI.R, P.Value, adj.P.Val) %>%
  left_join(rna, by = "gene") %>%
  arrange(adj.P.Val)

all_vs_rna %>% filter(is.na(RNA_log2FC)) %>% pull(analyte)   # no RNA match: add to `alias` if needed
write.csv(all_vs_rna, file.path(out_dir, "Step4_all_analytes_vs_RNA.csv"), row.names = FALSE)

# ---- 4d. Analytes with an RNA match, significant in protein OR RNA ----------
# Kept if protein q < sig_cut OR RNA padj < 0.05; dropped if neither is significant.
matched <- all_vs_rna %>%
  filter(!is.na(RNA_log2FC)) %>%
  mutate(prot_sig = adj.P.Val < sig_cut,
         rna_sig  = !is.na(RNA_padj) & RNA_padj < 0.05) %>%
  filter(prot_sig | rna_sig) %>%
  mutate(label     = ifelse(auto_symbol(analyte) == gene, analyte,
                            paste0(analyte, " (", gene, ")")),       # e.g. "MCP-1 (CCL2)"
         direction = ifelse(sign(prot_log2FC) == sign(RNA_log2FC),
                            "Same direction", "Opposite direction"),
         label     = factor(label, levels = label[order(prot_log2FC)]))

nrow(matched)                                    # number of analytes plotted
table(matched$direction)
count(matched, prot_sig, rna_sig)                # protein only / RNA only / both

# Protein-RNA agreement among these analytes
cor.test(matched$prot_log2FC, matched$RNA_log2FC, method = "spearman", exact = FALSE)

write.csv(matched, file.path(out_dir, "Step4_RNA_matched.csv"), row.names = FALSE)

# ---- 4e. Bar graph: protein and RNA log2FC side by side ---------------------
# Blue = protein (with 95% CI), orange = RNA.
# Solid = significant (protein q < sig_cut; RNA padj < 0.05), faded = not.
# Panels split analytes by whether protein and RNA change in the same direction.
matched_long <- bind_rows(
  matched %>% transmute(label, direction, assay = "Protein (Eve)",
                        log2FC = prot_log2FC, lo = CI.L, hi = CI.R,
                        sig = adj.P.Val < sig_cut),
  matched %>% transmute(label, direction, assay = "RNA (bulk RNA-seq)",
                        log2FC = RNA_log2FC, lo = NA_real_, hi = NA_real_,
                        sig = !is.na(RNA_padj) & RNA_padj < 0.05)
)

dodge <- position_dodge(width = 0.8)

p4 <- ggplot(matched_long, aes(x = log2FC, y = label, fill = assay, group = assay)) +
  geom_col(aes(alpha = sig), position = dodge, width = 0.7) +
  geom_errorbar(aes(xmin = lo, xmax = hi), position = dodge,
                width = 0.25, linewidth = 0.3, na.rm = TRUE) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Protein (Eve)" = "#2a78d6", "RNA (bulk RNA-seq)" = "#eb6834")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3),
                     labels = c(`TRUE` = "Significant", `FALSE` = "Not significant")) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey40"))) +
  labs(x = "log2 fold change, KO vs WT (positive = higher in KO)",
       y = NULL, fill = NULL, alpha = NULL,
       caption = paste0("Protein: limma, solid = q < ", sig_cut, ", error bars = 95% CI. ",
                        "RNA: solid = padj < 0.05; sign oriented so GPNMB is negative (higher in WT).")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"),
        legend.position = "top")

p4
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step4_RNA_matched.pdf", p4,
         width = 7, height = 0.25 * nrow(matched) + 2)




#finally let's check how many of our 80 analytes are in our olink data
'/Users/pennerc/Desktop/Desktop - BGS-STU23-138/2_Cognitive Decline Project/data sets/olink data 2/20220729_olink2_QC.csv'


# =============================================================================
# Step 5: PCA of WT and KO wells
# =============================================================================
# Input: log2 well-level values from Step 3 (`mat`), WT and KO wells only.
# Each analyte is centred and scaled to unit variance, so every analyte counts
# equally rather than the most variable ones dominating.

pca_mat <- mat[, grepl("^(WT|KO) ", colnames(mat))]
pca_mat <- pca_mat[complete.cases(pca_mat), ]          # PCA can't use missing values
pca_mat <- pca_mat[apply(pca_mat, 1, sd) > 0, ]        # constant analytes can't be scaled
nrow(pca_mat)                                          # analytes used in the PCA

pca     <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(pca$x[, 1:2], sample = rownames(pca$x)) %>%
  mutate(genotype = sub(" .*", "", sample),
         culture  = factor(paste(genotype, sub(".* (\\d)[ab]$", "\\1", sample)),
                           levels = c("WT 1", "WT 2", "KO 1", "KO 2")),
         well     = sub(".*\\d([ab])$", "\\1", sample))

rep_cols <- c("WT 1" = "#b2182b", "WT 2" = "#f08080",    # reds
              "KO 1" = "#00796b", "KO 2" = "#5ec4b6")    # teals

p5 <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = culture, shape = well)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 4) +
  scale_colour_manual(values = rep_cols, name = "Culture") +
  scale_shape_manual(values = c(a = 16, b = 17), name = "Well") +
  labs(x = paste0("PC1 (", var_exp[1], "% of variance)"),
       y = paste0("PC2 (", var_exp[2], "% of variance)")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p5
save_pdf("/Volumes/PC60/EveCytokineAnalyses/Step5_PCA_WT_vs_KO.pdf", p5, width = 6, height = 4.5)

# Analytes driving PC1 (largest absolute loadings)
head(sort(abs(pca$rotation[, "PC1"]), decreasing = TRUE), 10)














# =============================================================================
# Overlap with the Olink dataset
# =============================================================================
olink_file <- "/Users/pennerc/Desktop/Desktop - BGS-STU23-138/2_Cognitive Decline Project/data sets/olink data 2/20220729_olink2_QC.csv"
olink_cols <- names(read.csv(olink_file, nrows = 0, check.names = FALSE))   # header only
length(olink_cols)

# The 80 analysed analytes (tested in limma, media components excluded)
analytes80 <- setdiff(as.character(res$analyte), c("M-CSF", "IL-34"))
length(analytes80)                                   # should be 80

# Normalise names: drop punctuation/spaces, upper case; Greek -> letter or spelled out
norm       <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
greek_abbr <- function(x) chartr("αβγδω", "ABGDW", x)
greek_word <- function(x) {
  for (g in names(gw <- c(α = "alpha", β = "beta", γ = "gamma", δ = "delta", ω = "omega")))
    x <- gsub(g, gw[[g]], x, fixed = TRUE)
  x
}

olink_key <- norm(greek_word(olink_cols))

olink_overlap <- tibble(analyte = analytes80) %>%
  mutate(gene     = to_gene(analyte),                       # from Step 4b
         key_name = norm(greek_abbr(analyte)),
         key_word = norm(greek_word(analyte)),
         key_gene = norm(gene),
         olink_column = case_when(
           key_name %in% olink_key ~ olink_cols[match(key_name, olink_key)],
           key_word %in% olink_key ~ olink_cols[match(key_word, olink_key)],
           key_gene %in% olink_key ~ olink_cols[match(key_gene, olink_key)],
           TRUE ~ NA_character_),
         in_olink = !is.na(olink_column)) %>%
  select(analyte, gene, in_olink, olink_column)

sum(olink_overlap$in_olink)                                   # number of the 80 found in Olink
olink_overlap %>% filter(in_olink)                            # matched, with the Olink column name
olink_overlap %>% filter(!in_olink) %>% pull(analyte)         # not matched: skim for naming misses

write.csv(olink_overlap, file.path(out_dir, "Olink_overlap.csv"), row.names = FALSE)




