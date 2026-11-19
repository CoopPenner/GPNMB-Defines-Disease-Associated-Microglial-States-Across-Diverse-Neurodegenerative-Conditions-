# this script reproduces figure 6 
# manuscript version: Novemember 5th 

library(dplyr)
library(readxl)
library(MASS)
library(data.table)
library(ggplot2)
library(patchwork)
library(ggnewscale)
select <- dplyr::select

# this script works based on our box folder, this should be updated before submission
box <- c(file.path(Sys.getenv("USERPROFILE"), "Box"), "~/Library/CloudStorage/Box-Box", "~/Box")
box <- box[dir.exists(box)][1]
clin <- file.path(box, "Cooper + Marc", "GPNMB microglia paper", "Data", "Clinical")
out <- file.path(clin, "outputs", "Figure6")
dir.create(out, showWarnings = FALSE)
set.seed(1)

#read SEA-AD metadata, one row per donor
meta <- read_excel(file.path(clin, "inputs", "SEA-AD", "Cohort_Metadata.xlsx")) %>%
  transmute(donor = `Donor ID`,
            age = as.numeric(sub("90+", "90", `Age at Death`, fixed = TRUE)),
            sex = factor(Sex),
            pmi = as.numeric(PMI),
            apoe4 = factor(ifelse(grepl("4", `APOE Genotype`), "E4", "noE4"), c("noE4", "E4")),
            dementia = factor(`Cognitive Status`, c("No dementia", "Dementia")),
            Braak = factor(sub("Braak ", "", Braak), c("0", "I", "II", "III", "IV", "V", "VI"), ordered = TRUE),
            Thal = factor(sub("Thal ", "", Thal), c("0", "1", "2", "3", "4", "5"), ordered = TRUE),
            CERAD = factor(`CERAD score`, c("Absent", "Sparse", "Moderate", "Frequent"), ordered = TRUE),
            interval = as.numeric(`Interval from last CASI in months`))
# read cognitive slops (memory, language) and keep one per donor
cog <- read_excel(file.path(clin, "inputs", "SEA-AD", "Cognitive_Scores.xlsx")) %>%
  group_by(donor = `Donor ID`) %>%
  summarise(Memory = na.omit(slope_zmem0)[1], Language = na.omit(slope_zlan0)[1], .groups = "drop")
# read per-nucleus table, average each metric per donorxregion, keep samples w/ only 50 nuclei
sea <- fread(file.path(clin, "outputs", "SEAAD_cells.csv.gz"), select = c("donor", "region", "GPNMB", "TREM2", "brown", "yellow")) %>%
  as.data.frame() %>% rename(GNDAM = brown, Homeo = yellow)
sea_reg <- sea %>% group_by(donor, region) %>%
  summarise(n = n(), GPNMB = mean(GPNMB), TREM2 = mean(TREM2), GNDAM = mean(GNDAM), Homeo = mean(Homeo), .groups = "drop")
keep <- sea_reg %>% filter(region %in% c("MTG", "DFC", "MEC"), n >= 50) %>% select(donor, region)
# pool per donor, z-score gpnmb and trem2 across the 84 donors, and join metadata and slopes
sea_don <- sea %>% semi_join(keep, by = c("donor", "region")) %>% group_by(donor) %>%
  summarise(nuclei = n(), GPNMB = mean(GPNMB), TREM2 = mean(TREM2), .groups = "drop") %>%
  mutate(GPNMB = as.numeric(scale(GPNMB)), TREM2 = as.numeric(scale(TREM2))) %>%
  left_join(meta, by = "donor") %>% left_join(cog, by = "donor")
write.csv(sea_don, file.path(out, "SEAAD_donors_MTG_DFC_MEC.csv"), row.names = FALSE)

## PANELS C AND D (SEA_AD)
# ordinal logistic regression (polr) with age, sex, PMI and APOE4, for Braak/Thal/CERAD 
# for memor and language slopes, same but using a linear model
res <- NULL
for (g in c("GPNMB", "TREM2")) {
  for (y in c("Braak", "Thal", "CERAD")) {
    d <- sea_don
    d$stage <- droplevels(d[[y]])
    fit <- polr(as.formula(paste("stage ~", g, "+ age + sex + pmi + apoe4")), data = d, Hess = TRUE)
    s <- coef(summary(fit))[g, ]
    res <- rbind(res, data.frame(cohort = "SEA-AD", area = "Donor", gene = g, kind = "pathology", outcome = y,
                                 est = exp(s[[1]]), p = 2 * pnorm(-abs(s[[3]])), n = nrow(d)))
  }
  for (y in c("Memory", "Language")) {
    fit <- lm(as.formula(paste(y, "~", g, "+ age + sex + pmi + apoe4 + interval")), data = sea_don)
    s <- summary(fit)$coefficients[g, ]
    res <- rbind(res, data.frame(cohort = "SEA-AD", area = "Donor", gene = g, kind = "cognition", outcome = y,
                                 est = s[[1]], p = s[[4]], n = nobs(fit)))
  }
}

# PANEL B (SEA-AD)
# define region groups and then average z-scores and join metadata
sea_z <- sea_reg %>% filter(n >= 50) %>% group_by(region) %>%
  mutate(GPNMB = as.numeric(scale(GPNMB)), TREM2 = as.numeric(scale(TREM2)), GNDAM = as.numeric(scale(GNDAM)), Homeo = as.numeric(scale(Homeo))) %>%
  ungroup()
sea_areas <- list(MTL = c("LEC", "MEC", "HIP"), Cortex = c("ITG", "FI", "MTG", "STG", "AnG", "DFC", "V1C"))
for (a in names(sea_areas)) {
  d <- sea_z %>% filter(region %in% sea_areas[[a]]) %>% group_by(donor) %>%
    summarise(GPNMB = mean(GPNMB), TREM2 = mean(TREM2), GNDAM = mean(GNDAM), Homeo = mean(Homeo), .groups = "drop") %>%
    mutate(GPNMB = as.numeric(scale(GPNMB)), TREM2 = as.numeric(scale(TREM2)), GNDAM = as.numeric(scale(GNDAM)), Homeo = as.numeric(scale(Homeo))) %>%
    inner_join(meta, by = "donor") %>% left_join(cog, by = "donor")
  for (g in c("GNDAM", "GPNMB", "TREM2", "Homeo")) {
    for (y in c("Braak", "Thal", "CERAD")) {
      d$stage <- droplevels(d[[y]])
      fit <- polr(as.formula(paste("stage ~", g, "+ age + sex + pmi + apoe4")), data = d, Hess = TRUE)
      s <- coef(summary(fit))[g, ]
      res <- rbind(res, data.frame(cohort = "SEA-AD", area = a, gene = g, kind = "pathology", outcome = y,
                                   est = exp(s[[1]]), p = 2 * pnorm(-abs(s[[3]])), n = nrow(d)))
    }
    for (y in c("Memory", "Language")) {
      fit <- lm(as.formula(paste(y, "~", g, "+ age + sex + pmi + apoe4 + interval")), data = d)
      s <- summary(fit)$coefficients[g, ]
      res <- rbind(res, data.frame(cohort = "SEA-AD", area = a, gene = g, kind = "cognition", outcome = y,
                                   est = s[[1]], p = s[[4]], n = nobs(fit)))
    }
  }
}

# PANELS B and E (AMP-PD)
# read AMP-PD metadata, keep only inidividuals with a LB Braak stage and keep only "PD" diagnosis
# keep all regions, samples >50 nuclei
# models are same as above
amp_meta <- read.csv(file.path(clin, "inputs", "AMP-PD", "clinical_metadata.csv")) %>%
  filter(case_control_other_latest == "Case" | (case_control_other_latest == "Control" & !is.na(path_braak_lb))) %>%
  transmute(donor = participant_id,
            dx = factor(ifelse(case_control_other_latest == "Case", "PD", "Control"), c("Control", "PD")),
            age = age_at_baseline, sex = factor(sex), LB = path_braak_lb, NFT = path_braak_nft)

amp <- fread(file.path(clin, "outputs", "AMPPD_cells.csv.gz"), select = c("donor", "region", "GPNMB", "TREM2", "brown", "yellow")) %>%
  as.data.frame() %>% rename(GNDAM = brown, Homeo = yellow) %>%
  group_by(donor, region) %>%
  summarise(n = n(), GPNMB = mean(GPNMB), TREM2 = mean(TREM2), GNDAM = mean(GNDAM), Homeo = mean(Homeo), .groups = "drop") %>%
  filter(region %in% c("DMNX", "GPI", "PFC", "PMC", "PVC"), n >= 50) %>%
  inner_join(amp_meta, by = "donor") %>%
  group_by(region) %>%
  mutate(GPNMB = as.numeric(scale(GPNMB)), TREM2 = as.numeric(scale(TREM2)), GNDAM = as.numeric(scale(GNDAM)), Homeo = as.numeric(scale(Homeo))) %>%
  ungroup()
write.csv(amp, file.path(out, "AMPPD_donor_regions.csv"), row.names = FALSE)

amp_areas <- list(DMNX = "DMNX", PVC = "PVC", Brainstem = "DMNX", `Basal ganglia` = "GPI", Cortex = c("PFC", "PMC", "PVC"))
for (a in names(amp_areas)) {
  d <- amp %>% filter(region %in% amp_areas[[a]]) %>% group_by(donor) %>%
    summarise(GPNMB = mean(GPNMB), TREM2 = mean(TREM2), GNDAM = mean(GNDAM), Homeo = mean(Homeo),
              age = first(age), sex = first(sex), LB = first(LB), NFT = first(NFT), .groups = "drop")
  if (length(amp_areas[[a]]) > 1) {
    d <- d %>% mutate(GPNMB = as.numeric(scale(GPNMB)), TREM2 = as.numeric(scale(TREM2)), GNDAM = as.numeric(scale(GNDAM)), Homeo = as.numeric(scale(Homeo)))
  }
  genes <- c("GNDAM", "GPNMB", "TREM2", "Homeo")
  if (a %in% c("DMNX", "PVC")) genes <- c("GPNMB", "TREM2")
  for (g in genes) {
    for (y in c("LB", "NFT")) {
      dd <- d[!is.na(d[[y]]), ]
      dd$stage <- droplevels(factor(dd[[y]], ordered = TRUE))
      fit <- polr(as.formula(paste("stage ~", g, "+ age + sex")), data = dd, Hess = TRUE)
      s <- coef(summary(fit))[g, ]
      res <- rbind(res, data.frame(cohort = "AMP-PD", area = a, gene = g, kind = "pathology", outcome = y,
                                   est = exp(s[[1]]), p = 2 * pnorm(-abs(s[[3]])), n = nrow(dd)))
    }
  }
}
# Benjamini-Hochberg correction within each cohort x area x metric
# save as Figure6_stats.csv
res <- res %>% group_by(cohort, area, gene, kind) %>% mutate(fdr = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(stars = ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", "ns"))),
         cell = ifelse(fdr < 0.05, stars, ifelse(fdr <= 0.1, sprintf("%.2f", fdr), "ns")),
         label = paste0(ifelse(kind == "pathology", sprintf("OR %.2f", est), sprintf("β per SD %.3f", est)),
                        ifelse(fdr >= 0.05 & fdr < 0.1, sprintf(" (FDR = %.3f)", fdr), paste0(" (", stars, ")"))))
write.csv(res, file.path(out, "Figure6_stats.csv"), row.names = FALSE)

# establish theme and colors for the graphs
blue <- "#2a78d6"
grey <- "#f0efec"
red <- "#c93c3c"
dem_cols <- c("No dementia" = "#5B8DC9", "Dementia" = "#D9534F")
pd_cols <- c("Control" = "#5B8DC9", "PD" = "#D9534F")
th <- theme_classic(base_size = 7.5) +
  theme(axis.text = element_text(colour = "black", size = 6.5), axis.title = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3))
hm_th <- theme_minimal(base_size = 6.5) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black", size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text = element_text(face = "bold", size = 6, lineheight = 0.9),
        plot.title = element_text(face = "bold", size = 6.5, hjust = 0.5), legend.title = element_text(size = 5.5),
        legend.text = element_text(size = 5), legend.key.height = unit(2, "mm"), legend.key.width = unit(6, "mm"), legend.position = "bottom")
or_fill <- scale_fill_gradient2(low = blue, mid = grey, high = red, midpoint = 0, limits = c(-2, 2), breaks = c(-2, 0, 2),
                                labels = c("0.25", "1", "4"), oob = scales::squish, name = "Pathology score (OR per SD)",
                                guide = guide_colourbar(direction = "horizontal", title.position = "top", order = 1))

ad <- res %>% filter(cohort == "SEA-AD", area %in% c("MTL", "Cortex")) %>%
  mutate(block = factor(paste(area, ifelse(kind == "pathology", "Pathology", "Cognition"), sep = "\n"),
                        c("MTL\nPathology", "MTL\nCognition", "Cortex\nPathology", "Cortex\nCognition")),
         outcome = factor(outcome, c("Braak", "Thal", "CERAD", "Memory", "Language")),
         gene = factor(gene, c("Homeo", "TREM2", "GPNMB", "GNDAM")),
         bold = ifelse(fdr < 0.05, "bold", "plain"))
pd <- res %>% filter(cohort == "AMP-PD", area %in% c("Brainstem", "Basal ganglia", "Cortex")) %>%
  mutate(block = factor(area, c("Brainstem", "Basal ganglia", "Cortex"), c("Brainstem", "Basal\nganglia", "Cortex")),
         outcome = factor(outcome, c("LB", "NFT"), c("Braak LB", "Braak NFT")),
         gene = factor(gene, c("Homeo", "TREM2", "GPNMB", "GNDAM")),
         bold = ifelse(fdr < 0.05, "bold", "plain"))

# Build panel B
hm_ad <- ggplot(ad, aes(outcome, gene)) +
  geom_tile(data = filter(ad, kind == "pathology"), aes(fill = log2(est)), colour = "white", linewidth = 0.6) + or_fill +
  new_scale_fill() +
  geom_tile(data = filter(ad, kind == "cognition"), aes(fill = -est), colour = "white", linewidth = 0.6) +
  scale_fill_gradient2(low = blue, mid = grey, high = red, midpoint = 0, limits = c(-0.025, 0.025), breaks = c(-0.025, 0, 0.025),
                       labels = c("slower", "0", "faster"), oob = scales::squish, name = "Cognitive decline (beta per SD)",
                       guide = guide_colourbar(direction = "horizontal", title.position = "top", order = 2)) +
  geom_text(aes(label = cell, fontface = bold), size = 2) +
  facet_grid(. ~ block, scales = "free_x", space = "free_x") + scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Alzheimer's disease (SEA-AD)") + hm_th
hm_pd <- ggplot(pd, aes(outcome, gene)) +
  geom_tile(aes(fill = log2(est)), colour = "white", linewidth = 0.6) + or_fill +
  geom_text(aes(label = cell, fontface = bold), size = 2) +
  facet_grid(. ~ block, scales = "free_x", space = "free_x") + scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Parkinson's disease (AMP-PD)") + hm_th + theme(axis.text.y = element_blank())
ggsave(file.path(out, "Fig6B_heatmap.pdf"), hm_ad + hm_pd + plot_layout(widths = c(10, 6), guides = "collect") & theme(legend.position = "bottom"),
       width = 170, height = 62, units = "mm", device = cairo_pdf)

# Build panels C and D
titles <- c(Braak = "Tau (Braak stage)", Thal = "Amyloid (Thal phase)", CERAD = "Plaques (CERAD)")
cols <- list(Braak = c("#EAF1F8", "#08306B"), Thal = c("#FCE3D6", "#99000D"), CERAD = c("#EDF8E9", "#006D2C"))
pc <- list()
pdd <- list()
for (g in c("GPNMB", "TREM2")) {
  for (y in c("Braak", "Thal", "CERAD")) {
    lab <- res$label[res$area == "Donor" & res$gene == g & res$outcome == y]
    d <- sea_don %>% mutate(stage = droplevels(factor(.data[[y]], ordered = FALSE)), val = .data[[g]])
    p <- ggplot(d, aes(stage, val)) +
      geom_boxplot(aes(fill = stage), outlier.shape = NA, width = 0.7, linewidth = 0.3, show.legend = FALSE) +
      geom_jitter(aes(colour = dementia), width = 0.18, height = 0, size = 0.7) +
      scale_fill_manual(values = colorRampPalette(cols[[y]])(nlevels(d$stage))) +
      scale_colour_manual(values = dem_cols, name = NULL) + guides(colour = guide_legend(override.aes = list(size = 2))) +
      coord_cartesian(ylim = c(-2.6, 3.8)) +
      labs(x = NULL, y = if (y == "Braak") paste("Microglial", g, "(z)") else NULL,
           title = if (g == "GPNMB") paste0(titles[[y]], "\n", lab) else lab) + th
    if (y == "CERAD") p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    pc[[paste(g, y)]] <- p
  }
  for (y in c("Memory", "Language")) {
    lab <- res$label[res$area == "Donor" & res$gene == g & res$outcome == y]
    d <- sea_don %>% filter(!is.na(.data[[y]])) %>% mutate(val = .data[[g]], slope = .data[[y]])
    pdd[[paste(g, y)]] <- ggplot(d, aes(val, slope)) +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey80", linewidth = 0.45) +
      geom_point(aes(colour = dementia), size = 0.8) +
      scale_colour_manual(values = dem_cols, name = NULL) + guides(colour = guide_legend(override.aes = list(size = 2))) +
      labs(x = paste("Microglial", g, "(z)"), y = if (y == "Memory") "Annual slope (z units)" else NULL,
           title = if (g == "GPNMB") paste0(y, "\n", lab) else lab) + th
  }
}
ggsave(file.path(out, "Fig6C_SEAAD_pathology.pdf"), wrap_plots(pc, nrow = 2, widths = c(7, 6, 4.2), guides = "collect") & theme(legend.position = "bottom"),
       width = 120, height = 118, units = "mm", device = cairo_pdf)
ggsave(file.path(out, "Fig6D_SEAAD_cognition.pdf"), wrap_plots(pdd, nrow = 2, guides = "collect") & theme(legend.position = "bottom"),
       width = 95, height = 118, units = "mm", device = cairo_pdf)

# build panel E
where <- c(DMNX = "Brainstem (DMNX)", PVC = "Primary visual cortex (PVC)")
pe <- list()
for (g in c("GPNMB", "TREM2")) {
  for (r in c("DMNX", "PVC")) {
    for (y in c("LB", "NFT")) {
      lab <- res$label[res$cohort == "AMP-PD" & res$area == r & res$gene == g & res$outcome == y]
      d <- amp %>% filter(region == r, !is.na(.data[[y]])) %>% mutate(val = pmin(.data[[g]], 3.5))
      if (y == "LB") {
        d$stage <- factor(d$LB, 0:6)
        fills <- colorRampPalette(c("#FFF7BC", "#8C2D04"))(7)
        top <- "Lewy bodies (Braak LB)"
      } else {
        d$stage <- factor(d$NFT, 0:3, c("0", "I", "II", "III"))
        fills <- colorRampPalette(c("#EAF1F8", "#08306B"))(4)
        top <- "Tau (Braak NFT)"
      }
      pe[[paste(g, r, y)]] <- ggplot(d, aes(stage, val)) +
        geom_boxplot(aes(fill = stage), outlier.shape = NA, width = 0.7, linewidth = 0.3, show.legend = FALSE) +
        geom_jitter(aes(colour = dx), width = 0.18, height = 0, size = 0.7) +
        scale_fill_manual(values = fills, drop = FALSE) + scale_x_discrete(drop = FALSE) +
        scale_colour_manual(values = pd_cols, name = NULL, labels = c(Control = "Control", PD = "Parkinson's")) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        coord_cartesian(ylim = c(-3.3, 4.4)) +
        labs(x = NULL, y = if (r == "DMNX" && y == "LB") paste("Microglial", g, "(z)") else NULL,
             title = if (g == "GPNMB") paste0(where[[r]], "\n", top, "\n", lab) else lab) + th
    }
  }
}
ggsave(file.path(out, "Fig6E_AMPPD_DMNX_vs_PVC.pdf"), wrap_plots(pe, nrow = 2, widths = c(7, 4.5, 7, 4.5), guides = "collect") & theme(legend.position = "bottom"),
       width = 210, height = 125, units = "mm", device = cairo_pdf)
