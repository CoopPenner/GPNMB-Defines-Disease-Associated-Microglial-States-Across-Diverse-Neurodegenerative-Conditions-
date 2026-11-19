# the goasl of this is to build SEA-AD clinical metadata tables + UCELL microglia scores 

library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(data.table)


# linking the script ot the box... 
box <- c(file.path(Sys.getenv("USERPROFILE"), "Box"), "~/Library/CloudStorage/Box-Box", "~/Box")
box <- box[dir.exists(box)][1]
root <- file.path(box, "Cooper + Marc", "GPNMB microglia paper", "Data", "Clinical")
input <- file.path(root, "inputs/SEA-AD")
shared <- file.path(root, "inputs/shared")
output <- file.path(root, "outputs")

modules <- read.csv(file.path(shared, "Consensus_modules_all_members_long.csv")) # load cooper's modules <3
modules <- modules[modules$module != "blue", ] # remove the neuronal ambient module
# how many genes to use to define a module?? I am putting here 20 - because that is what cooper did for Figure 1. This can be later changed.
top20 <- modules %>% group_by(module) %>% arrange(desc(kME), .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup() 

signatures <- split(top20$gene_name, top20$module)
signatures$brown_noGPNMB <- setdiff(signatures$brown, "GPNMB") # remove GPNMB from GNDAM module, save as new module

## Load SEA-AD:
obj <- readRDS(file.path(input, "SEA-AD_Microglia_v2.rds"))
if (!"data" %in% Layers(obj)) obj <- NormalizeData(obj)
# add UCELL module scores
scores <- ScoreSignatures_UCell(GetAssayData(obj, layer = "counts"), features = signatures, maxRank = 1500,
                                chunk.size = 2000, BPPARAM = BiocParallel::SerialParam(), name = "")
# pull GPNMB & TREM2 expression
expression <- FetchData(obj, c("GPNMB", "TREM2"), layer = "data")
umap <- Embeddings(obj, "umap")
# this builds a table of a per cell basis info, then I rmeove it to save disk space
cells <- data.frame(cell = colnames(obj),
                    donor = obj$sample_id,
                    region = obj$brain_region,
                    cluster = as.character(obj$seurat_clusters),
                    UMAP_1 = umap[, 1],
                    UMAP_2 = umap[, 2],
                    GPNMB = expression$GPNMB,
                    TREM2 = expression$TREM2,
                    scores[colnames(obj), names(signatures)])
rm(obj)
gc()

# I want to keep apart GNDAM_noGPNMB so that it does not compete w GUNDAM when assigning module identity 
module_names <- setdiff(names(signatures), "brown_noGPNMB")
# assign every nucleus a score for each module, z-score, then grab the top z-score
# & only assign if winning module is 1SD above average; otherwise we leave unassigned (gray)
assign_module <- function(score_table, brown_column) {
  score_matrix <- as.matrix(score_table[, c(setdiff(module_names, "brown"), brown_column)])
  colnames(score_matrix)[ncol(score_matrix)] <- "brown"
  z <- scale(score_matrix)
  best <- max.col(z, ties.method = "first")
  best_z <- z[cbind(seq_len(nrow(z)), best)]
  ifelse(best_z >= 1, colnames(z)[best], "unassigned")
}

# this can be refined... i would like to group modules into broad biological programs....
# this is something that I can make this analysis stand out... consider whether other modules may contribute to pro-inf...
broad_state <- function(module) {
  case_when(module %in% c("midnightblue", "magenta") ~ "Pro-inflammatory",
            module %in% c("yellow", "black") ~ "Homeostatic",
            module == "brown" ~ "GNDAM",
            TRUE ~ "Other")
}

main_regions <- c("MTG", "DFC", "MEC")
in_main <- cells$region %in% main_regions

cells$module_3regions <- NA
cells$module_3regions_withGPNMB <- NA
cells$module_3regions[in_main] <- assign_module(cells[in_main, ], "brown_noGPNMB") 
cells$module_3regions_withGPNMB[in_main] <- assign_module(cells[in_main, ], "brown") # to compare GNDAM vs GNDAM minus GPNMB
cells$module_allregions <- assign_module(cells, "brown_noGPNMB")
cells$state_3regions <- ifelse(in_main, broad_state(cells$module_3regions), NA)
cells$state_allregions <- broad_state(cells$module_allregions)
# write the massive cells csv
fwrite(cells, file.path(output, "SEAAD_cells.csv.gz"))

# add metadata from file
metadata <- read_excel(file.path(input, "Cohort_Metadata.xlsx"))
lewy_body_stage <- c("Not Identified (olfactory bulb not assessed)" = "None",
                     "Not Identified (olfactory bulb assessed)" = "None",
                     "Olfactory bulb only" = "Focal",
                     "Brainstem-predominant" = "Focal",
                     "Amygdala-predominant" = "Focal",
                     "Limbic (Transitional)" = "Limbic",
                     "Neocortical (Diffuse)" = "Neocortical")
clinical <- metadata %>%
  transmute(donor = `Donor ID`,
            age = as.numeric(sub("90+", "90", `Age at Death`, fixed = TRUE)),
            sex = Sex,
            pmi = as.numeric(PMI),
            rin = as.numeric(RIN),
            apoe_genotype = `APOE Genotype`,
            apoe_e4_dose = lengths(regmatches(`APOE Genotype`, gregexpr("4", `APOE Genotype`))),
            apoe_e4 = ifelse(apoe_e4_dose > 0, "E4", "noE4"),
            dementia = `Cognitive Status`,
            braak = sub("Braak ", "", Braak),
            braak_n = match(braak, c("0", "I", "II", "III", "IV", "V", "VI")) - 1,
            thal = as.numeric(sub("Thal ", "", Thal)),
            cerad = `CERAD score`,
            cerad_n = match(cerad, c("Absent", "Sparse", "Moderate", "Frequent")) - 1,
            adnc = `Overall AD neuropathological Change`,
            adnc_n = match(adnc, c("Not AD", "Low", "Intermediate", "High")) - 1,
            lewy_body = `Highest Lewy Body Disease`,
            lewy_body_stage = lewy_body_stage[`Highest Lewy Body Disease`],
            late = LATE,
            caa = `Overall CAA Score`,
            casi_interval_months = as.numeric(`Interval from last CASI in months`))

cognition <- read_excel(file.path(input, "Cognitive_Scores.xlsx")) %>%
  group_by(donor = `Donor ID`) %>%
  summarise(slope_memory = first(na.omit(slope_zmem0)),
            slope_executive = first(na.omit(slope_zexf0)),
            slope_language = first(na.omit(slope_zlan0)),
            slope_visuospatial = first(na.omit(slope_zvsp0)),
            slope_casi = first(na.omit(slope_zcasi_irt0)),
            .groups = "drop")

genotype <- read.csv(file.path(input, "rs199347_WGS_genotypes.csv")) %>%
  transmute(donor,
            rs199347 = c("0/0" = "AA", "0/1" = "AG", "1/1" = "GG")[GT],
            rs199347_A_dosage = c("0/0" = 2, "0/1" = 1, "1/1" = 0)[GT]) # AA homozygous risk allele

staging <- read.csv(file.path(input, "SEAAD_region_staging.csv")) # merges clinical metadata, rs199347 genotype, and cognitive slopes into one row pr donor

donor_info <- clinical %>% left_join(genotype, by = "donor") %>% left_join(cognition, by = "donor")

# this gets the avg values 
score_means <- function(data) {
  data %>% summarise(n_nuclei = n(),
                     GPNMB = mean(GPNMB),
                     GPNMB_pct_positive = 100 * mean(GPNMB > 0),
                     TREM2 = mean(TREM2),
                     GNDAM = mean(brown),
                     GNDAM_noGPNMB = mean(brown_noGPNMB),
                     HOMEO = mean(yellow),
                     HOMEO2 = mean(black),
                     across(all_of(module_names), mean, .names = "score_{.col}"),
                     .groups = "drop")
}

# whereas this gets the prportion of icroglia in X state
state_percent <- function(data, column, prefix) {
  data %>% count(donor, .data[[column]]) %>%
    group_by(donor) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = all_of(column), values_from = pct, names_prefix = prefix, values_fill = 0)
}

# Donor_region = one row per donor per region, 10 regions = this is to get amyloid vs tau specific effects (see ppt)
donor_region <- cells %>% group_by(donor, region) %>% score_means()
region_states <- cells %>% count(donor, region, state_allregions) %>%
  group_by(donor, region) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% select(-n) %>%
  pivot_wider(names_from = state_allregions, values_from = pct, names_prefix = "pct_state_", values_fill = 0)
region_modules <- cells %>% count(donor, region, module_allregions) %>%
  group_by(donor, region) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% select(-n) %>%
  pivot_wider(names_from = module_allregions, values_from = pct, names_prefix = "pct_module_", values_fill = 0)
donor_region <- donor_region %>%
  left_join(region_states, by = c("donor", "region")) %>%
  left_join(region_modules, by = c("donor", "region")) %>%
  left_join(staging %>% select(region, tau_onset, amy_onset, region_group = group), by = "region") %>%
  left_join(donor_info, by = "donor") %>%
  mutate(tau_local = as.integer(braak_n >= tau_onset),
         amyloid_local = as.integer(thal >= amy_onset),
         local_status = case_when(tau_local == 0 & amyloid_local == 0 ~ "Neither",
                                  tau_local == 0 & amyloid_local == 1 ~ "Amyloid only",
                                  tau_local == 1 & amyloid_local == 0 ~ "Tau only",
                                  TRUE ~ "Both"),
         tau_dose = pmax(0, braak_n - tau_onset + 1),
         amyloid_dose = pmax(0, thal - amy_onset + 1))

# Donor sheet = one row per donor, pooling the 3 regions together
main_cells <- cells %>% filter(region %in% main_regions)
donor <- main_cells %>% group_by(donor) %>% score_means() %>%
  mutate(across(c(GPNMB, TREM2, GNDAM, GNDAM_noGPNMB, HOMEO, HOMEO2), ~ as.numeric(scale(.x)), .names = "{.col}_z")) %>%
  left_join(state_percent(main_cells, "state_3regions", "pct_state_"), by = "donor") %>%
  left_join(state_percent(main_cells, "module_3regions", "pct_module_"), by = "donor") %>%
  left_join(state_percent(main_cells, "module_3regions_withGPNMB", "pct_module_withGPNMB_") %>%
              select(donor, pct_module_withGPNMB_brown), by = "donor") %>%
  left_join(main_cells %>% group_by(donor) %>% summarise(pct_cluster5 = 100 * mean(cluster == "5"), .groups = "drop"), by = "donor") %>%
  left_join(donor_info, by = "donor")

signature_table <- data.frame(signature = rep(names(signatures), lengths(signatures)),
                              gene = unlist(signatures),
                              rank = unlist(lapply(signatures, seq_along)))

write_xlsx(list(donor = donor,
                donor_region = donor_region,
                signatures = signature_table),
           file.path(output, "SEAAD_dataset.xlsx"))
