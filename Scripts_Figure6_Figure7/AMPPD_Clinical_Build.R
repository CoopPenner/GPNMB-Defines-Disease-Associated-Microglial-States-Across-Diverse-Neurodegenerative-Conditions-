# same but for the AMP_PD dataset

library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(writexl)
library(data.table)

# linking the script to the box...
box <- c(file.path(Sys.getenv("USERPROFILE"), "Box"), "~/Library/CloudStorage/Box-Box", "~/Box")
box <- box[dir.exists(box)][1]
root <- file.path(box, "Cooper + Marc", "GPNMB microglia paper", "Data", "Clinical")
input <- file.path(root, "inputs/AMP-PD")
shared <- file.path(root, "inputs/shared")
output <- file.path(root, "outputs")

modules <- read.csv(file.path(shared, "Consensus_modules_all_members_long.csv")) # load cooper's modules <3
modules <- modules[modules$module != "blue", ] # remove the neuronal ambient module
# how many genes to use to define a module?? I am putting here 20 - because that is what cooper did for Figure 1. This can be later changed.
top20 <- modules %>% group_by(module) %>% arrange(desc(kME), .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup()

signatures <- split(top20$gene_name, top20$module)
signatures$brown_noGPNMB <- setdiff(signatures$brown, "GPNMB") # remove GPNMB from GNDAM module, save as new module

## Load AMP-PD:
obj <- readRDS(file.path(input, "AMPPD_Microglia_v1.rds"))
if (!"data" %in% Layers(obj)) obj <- NormalizeData(obj)
signatures <- lapply(signatures, function(genes) intersect(genes, rownames(obj))) # keep only the module genes that are in the AMP-PD object
# add UCELL module scores
scores <- ScoreSignatures_UCell(GetAssayData(obj, layer = "counts"), features = signatures, maxRank = 1500,
                                chunk.size = 2000, BPPARAM = BiocParallel::SerialParam(), name = "")
# pull GPNMB & TREM2 expression
expression <- FetchData(obj, c("GPNMB", "TREM2"), layer = "data")
umap <- Embeddings(obj, "umap")
# this builds a table of a per cell basis info, then I remove it to save disk space
cells <- data.frame(cell = colnames(obj),
                    donor = obj$participant_id,
                    region = obj$brain_region,
                    cluster = as.character(obj$seurat_clusters),
                    UMAP_1 = umap[, 1],
                    UMAP_2 = umap[, 2],
                    GPNMB = expression$GPNMB,
                    TREM2 = expression$TREM2,
                    scores[colnames(obj), names(signatures)])
rm(obj)
gc()

# I want to keep apart GNDAM_noGPNMB so that it does not compete w GNDAM when assigning module identity
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

# in AMP-PD we use all regions together (DMNX, GPI, PFC, PMC, PVC), so only one assignment
cells$module <- assign_module(cells, "brown_noGPNMB")
cells$module_withGPNMB <- assign_module(cells, "brown") # to compare GNDAM vs GNDAM minus GPNMB
cells$state <- broad_state(cells$module)
# write the massive cells csv
fwrite(cells, file.path(output, "AMPPD_cells.csv.gz"))

# add metadata from file
clinical <- read.csv(file.path(input, "clinical_metadata.csv")) %>%
  transmute(donor = participant_id,
            group = case_control_other_latest,
            diagnosis = diagnosis_latest,
            PD = as.integer(case_control_other_latest == "Case"),
            age = age_at_baseline,
            sex,
            braak_lb = ifelse(is.na(path_braak_lb) & case_control_other_latest == "Control", 0, path_braak_lb), # controls with no LB score = stage 0
            braak_nft = path_braak_nft,
            dlb_likelihood = sub(" Likelihood", "", path_dlb_prob),
            mmse = mmse_total_score,
            updrs3 = updrs3_motor_examination_score,
            hoehn_yahr = upd2hy_hoehn_and_yahr_stage,
            age_at_diagnosis)

genotype <- read.csv(file.path(input, "AMPPD_rs199347_calls_merged.csv")) %>%
  transmute(donor = sub("-BLM0.*", "", sample),
            rs199347,
            rs199347_A_dosage = c(AA = 2, AG = 1, GG = 0)[rs199347], # AA homozygous risk allele
            rs199347_call_method = method)

donor_info <- clinical %>% left_join(genotype, by = "donor") # merges clinical metadata and rs199347 genotype into one row per donor

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

# whereas this gets the proportion of microglia in X state
state_percent <- function(data, column, prefix) {
  data %>% count(donor, .data[[column]]) %>%
    group_by(donor) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = all_of(column), values_from = pct, names_prefix = prefix, values_fill = 0)
}

# Donor_region = one row per donor per region, 5 regions = this is to get region-specific effects 
donor_region <- cells %>% group_by(donor, region) %>% score_means()
region_states <- cells %>% count(donor, region, state) %>%
  group_by(donor, region) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% select(-n) %>%
  pivot_wider(names_from = state, values_from = pct, names_prefix = "pct_state_", values_fill = 0)
region_modules <- cells %>% count(donor, region, module) %>%
  group_by(donor, region) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% select(-n) %>%
  pivot_wider(names_from = module, values_from = pct, names_prefix = "pct_module_", values_fill = 0)
donor_region <- donor_region %>%
  left_join(region_states, by = c("donor", "region")) %>%
  left_join(region_modules, by = c("donor", "region")) %>%
  left_join(donor_info, by = "donor")

# Donor sheet = one row per donor, pooling all regions together
donor <- cells %>% group_by(donor) %>% score_means() %>%
  mutate(across(c(GPNMB, TREM2, GNDAM, GNDAM_noGPNMB, HOMEO, HOMEO2), ~ as.numeric(scale(.x)), .names = "{.col}_z")) %>%
  left_join(state_percent(cells, "state", "pct_state_"), by = "donor") %>%
  left_join(state_percent(cells, "module", "pct_module_"), by = "donor") %>%
  left_join(state_percent(cells, "module_withGPNMB", "pct_module_withGPNMB_") %>%
              select(donor, pct_module_withGPNMB_brown), by = "donor") %>%
  left_join(donor_info, by = "donor")

## Load all AMP-PD cells (not only microglia); this is the small per cell table with GPNMB counts, built from the big h5ad
all_cells <- readRDS(file.path(input, "AMPPD_global_celltable.rds")) %>%
  filter(ct_tier != "dropped") %>%
  mutate(cell_class = case_when(cell_type %in% c("EN", "IN", "DMNX_Neu", "GPI_Neu") ~ "Neurons", # combine all neuron labels
                                cell_type == "Astro" ~ "Astrocytes",
                                cell_type == "OPC" ~ "OPCs",
                                cell_type == "Microglia" ~ "Microglia",
                                cell_type == "Oligo" ~ "Oligodendrocytes",
                                TRUE ~ "Other"))

# pseudobulk = sum all GPNMB reads of each cell type per donor per region, then normalize (CPM)
# pct_of_cells = how much of each region is that cell type (cell type composition)
celltype_pseudobulk <- all_cells %>%
  group_by(donor = participant_id, region = brain_region, cell_class) %>%
  summarise(n_cells = n(),
            GPNMB_umi = sum(GPNMB),
            library_size = sum(total_counts),
            GPNMB_pct_positive = 100 * mean(GPNMB > 0),
            .groups = "drop") %>%
  mutate(GPNMB_log2cpm = log2(1e6 * GPNMB_umi / library_size + 1)) %>%
  group_by(donor, region) %>% mutate(pct_of_cells = 100 * n_cells / sum(n_cells)) %>% ungroup() %>%
  left_join(donor_info, by = "donor")

signature_table <- data.frame(signature = rep(names(signatures), lengths(signatures)),
                              gene = unlist(signatures),
                              rank = unlist(lapply(signatures, seq_along)))

write_xlsx(list(donor = donor,
                donor_region = donor_region,
                celltype_pseudobulk = celltype_pseudobulk,
                signatures = signature_table),
           file.path(output, "AMPPD_dataset.xlsx"))
