options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(purrr)
library(stringr)
library(stats)
library(ggrepel)

data <- read.csv("./Cleaned/20200106_Puglielli_ATaseKO_tc_tidy_peptide_stoich_output.csv")

### Preparing the data for ANOVA analysis
ms <- data %>% select(genotype, replicate, PG.ProteinGroups, PG.Genes, PEP.StrippedSequence,
                      EG.ModifiedPeptide, k_site, stoich_corr)
ms$condition_num <- 0
ms$condition_num[which(ms$genotype == "ATaseWT")] <- 1
ms$condition_num[which(ms$genotype == "ATase1KO")] <- 2
ms$condition_num[which(ms$genotype == "ATase2KO")] <- 3

ms$condition_num <- as.numeric(ms$condition_num)

ms$name <- paste(ms$PG.ProteinGroups, ms$k_site, ms$EG.ModifiedPeptide, sep = "_")

mss <- split(ms, paste(ms$PG.ProteinGroups, ms$k_site, ms$EG.ModifiedPeptide, sep = "_"))

################################

tmpfn <- function(dat) {
  dat <- dat %>% filter(condition_num == 1 | condition_num == 2)
  if(length(unique(dat$condition_num)) < 2)
    return(NULL)
  fit <- aov(stoich_corr ~ factor(condition_num), dat)
  Tukey <- TukeyHSD(fit)
}
test_ATase1KO <- mss %>%
  map(tmpfn) %>%
  compact()
p_adj_ATase1KO <- data.frame(sapply(unlist(test_ATase1KO, recursive = FALSE), '[',,"p adj"))
names(p_adj_ATase1KO)[1] <- "p_adj_ATase1KO"
p_adj_ATase1KO$names <- rownames(p_adj_ATase1KO)
p_adj_ATase1KO <- separate(p_adj_ATase1KO, col = names, sep = ".fac", into = c("names", "factor"))
p_adj_ATase1KO$factor <- NULL

###########################################################################

tmpfn <- function(dat) {
  dat <- dat %>% filter(condition_num == 1 | condition_num == 3)
  if(length(unique(dat$condition_num)) < 2)
    return(NULL)
  fit <- aov(stoich_corr ~ factor(condition_num), dat)
  Tukey <- TukeyHSD(fit)
}
test_ATase2KO <- mss %>%
  map(tmpfn) %>%
  compact()
p_adj_ATase2KO <- data.frame(sapply(unlist(test_ATase2KO, recursive = FALSE), '[',,"p adj"))
names(p_adj_ATase2KO)[1] <- "p_adj_ATase2KO"
p_adj_ATase2KO$names <- rownames(p_adj_ATase2KO)
p_adj_ATase2KO <- separate(p_adj_ATase2KO, col = names, sep = ".fac", into = c("names", "factor"))
p_adj_ATase2KO$factor <- NULL

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

p_adj_ATase1KO[is.nan(p_adj_ATase1KO)] <- NA
p_adj_ATase2KO[is.nan(p_adj_ATase2KO)] <- NA

#########################################################
data <- data %>% 
  select(genotype, replicate, everything())  %>% 
  arrange(PG.ProteinGroups, k_site, k_count, genotype, replicate)

summary <- data %>% 
  group_by(genotype, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           k_site, k_count, PEP.StrippedSequence, EG.ModifiedPeptide) %>% 
  summarise(stoich_mean = mean(stoich_corr, na.rm = TRUE), 
            stoich_sd = sd(stoich_corr, na.rm = TRUE),
            stoich_n = n()) %>% 
  arrange(PG.ProteinGroups, k_site, genotype) %>% 
  ungroup()

summary$names <- paste(summary$PG.ProteinGroups, summary$k_site, summary$EG.ModifiedPeptide, sep = "_")

###########################################################################

merge_ATaseKO <- merge(p_adj_ATase1KO, p_adj_ATase2KO, all = TRUE)

data_stats <- merge(summary, merge_ATaseKO, by = "names", all = TRUE)

data_stats <- data_stats %>% subset(stoich_mean <= 1)

rm(p_adj_ATase1KO, p_adj_ATase2KO, merge_ATaseKO, ms, mss, summary, 
   test_ATase1KO, test_ATase2KO);gc()


####
wide_mean <- data_stats %>%
  select(PG.ProteinGroups, k_site, PG.Genes, PG.ProteinDescriptions,
         EG.ModifiedPeptide, PEP.StrippedSequence, k_count, genotype,
         p_adj_ATase1KO, p_adj_ATase2KO, stoich_mean) %>% 
  spread(genotype, stoich_mean)
names(wide_mean)[10] <- "ATase1KO_mean"
names(wide_mean)[11] <- "ATase2KO_mean"
names(wide_mean)[12] <- "ATaseWT_mean"

wide_sd <- data_stats %>%
  select(PG.ProteinGroups, k_site, PG.Genes, PG.ProteinDescriptions,
         EG.ModifiedPeptide, PEP.StrippedSequence, k_count, genotype, 
         p_adj_ATase1KO, p_adj_ATase2KO, stoich_sd) %>% 
  spread(genotype, stoich_sd)
names(wide_sd)[10] <- "ATase1KO_sd"
names(wide_sd)[11] <- "ATase2KO_sd"
names(wide_sd)[12] <- "ATaseWT_sd"

wide_n <- data_stats %>%
  select(PG.ProteinGroups, k_site, PG.Genes, PG.ProteinDescriptions, 
         EG.ModifiedPeptide, PEP.StrippedSequence, k_count, genotype, 
         p_adj_ATase1KO, p_adj_ATase2KO, stoich_n) %>% 
  spread(genotype, stoich_n)
names(wide_n)[10] <- "ATase1KO_n"
names(wide_n)[11] <- "ATase2KO_n"
names(wide_n)[12] <- "ATaseWT_n"

wide <- merge(wide_mean, wide_sd, all = TRUE)
wide <- merge(wide, wide_n, all = TRUE)

rm(wide_mean, wide_sd, wide_n); gc()

wide$ATase1KO_WT <- wide$ATase1KO_mean - wide$ATaseWT_mean
wide$ATase2KO_WT <- wide$ATase2KO_mean - wide$ATaseWT_mean

wide$ATase1KO_WT <- round(wide$ATase1KO_WT, 4)
wide$ATase2KO_WT <- round(wide$ATase2KO_WT, 4)
wide$p_adj_ATase1KO <- round(wide$p_adj_ATase1KO, 4)
wide$p_adj_ATase2KO <- round(wide$p_adj_ATase2KO, 4)
wide$ATaseWT_mean <- round(wide$ATaseWT_mean, 4)
wide$ATase1KO_mean <- round(wide$ATase1KO_mean, 4)
wide$ATase2KO_mean <- round(wide$ATase2KO_mean, 4)


###################
H_present <- read.csv("H_pep_present_addition.csv")
H_present$n_Hsum_include <- NULL
H_present$n_Lsum_isnot_NA <- NULL

H_present_wide <- H_present %>%
  select(PG.ProteinGroups, k_site, PG.Genes, PG.ProteinDescriptions,
         EG.ModifiedPeptide, PEP.StrippedSequence, k_count, genotype,
         stoich_assesment) %>% 
  spread(genotype, stoich_assesment)
names(H_present_wide)[8] <- "ATase1KO_st_check"
names(H_present_wide)[9] <- "ATase2KO_st_check"
names(H_present_wide)[10] <- "ATaseWT_st_check"


wide <- merge(wide, H_present_wide, all.x = TRUE)

wide$ATase1KO_st_check[which(!is.na(wide$ATase1KO_mean))] <- "calc_stoich"
wide$ATase1KO_st_check[which((is.na(wide$ATase1KO_st_check) & is.na(wide$ATase1KO_mean)))] <- "missing_peptide"
wide$ATase1KO_st_check[which(((wide$ATase1KO_st_check == "calc_stoich") & is.na(wide$ATase1KO_mean)))] <- "missing_peptide"

wide$ATase2KO_st_check[which((!is.na(wide$ATase2KO_mean)))] <- "calc_stoich"
wide$ATase2KO_st_check[which((is.na(wide$ATase2KO_st_check) & is.na(wide$ATase2KO_mean)))] <- "missing_peptide"
wide$ATase2KO_st_check[which(((wide$ATase2KO_st_check == "calc_stoich") & is.na(wide$ATase2KO_mean)))] <- "missing_peptide"

wide$ATaseWT_st_check[which((!is.na(wide$ATaseWT_mean)))] <- "calc_stoich"
wide$ATaseWT_st_check[which((is.na(wide$ATaseWT_st_check) & is.na(wide$ATaseWT_mean)))] <- "missing_peptide"
wide$ATaseWT_st_check[which(((wide$ATaseWT_st_check == "calc_stoich") & is.na(wide$ATaseWT_mean)))] <- "missing_peptide"

length(which(wide$ATaseWT_n >= 3 & is.na(wide$ATase1KO_mean)))
length(which(wide$ATaseWT_n >= 3 & wide$ATase1KO_st_check == "H_pep_present"))

length(which(wide$ATaseWT_n >= 3 & is.na(wide$ATase2KO_mean)))
length(which(wide$ATaseWT_n >= 3 & wide$ATase2KO_st_check == "H_pep_present"))

length(which(wide$ATase1KO_n >= 3 & is.na(wide$ATaseWT_mean)))
length(which(wide$ATase1KO_n >= 3 & wide$ATaseWT_st_check == "H_pep_present"))

length(which(wide$ATase2KO_n >= 3 & is.na(wide$ATaseWT_mean)))
length(which(wide$ATase2KO_n >= 3 & wide$ATaseWT_st_check == "H_pep_present"))

length(which((wide$ATase1KO_n >= 3 | wide$ATase2KO_n >= 3) & is.na(wide$ATaseWT_mean)))
length(which((wide$ATase1KO_n >= 3 | wide$ATase2KO_n >= 3) & wide$ATaseWT_st_check == "H_pep_present"))


# ##########
# Loading the ER protein list -  need to convert from human to mouse...
# human_ER <- read.csv(file = "subcell_location_Endoplasmic.csv")
# human_ER_small <- human_ER %>% select(Gene, Ensembl, Gene.description, Uniprot, UniProt.evidence) %>%
#   subset(Uniprot != "")

human_ER_uniprot <- read.csv(file = "ER_updated.csv")
human_ER_uniprot$PG.Genes_human <- as.character(unlist(lapply(human_ER_uniprot$Gene.names, function(x){
  unlist(strsplit(x, split = " "))[1]
})))
human_ER_uniprot$PG.Genes <- str_to_title(human_ER_uniprot$PG.Genes_human, locale = "en")
human_ER_uniprot <- human_ER_uniprot %>%
  select(PG.Genes) %>%
  distinct() %>%
  mutate(ER_proteinatlas = "TRUE")

uniprot_localization <- read.csv("./uniprot-proteome_mous_subcell_location.csv")
uniprot_localization$ER <- FALSE
uniprot_localization$ER[grep("ER", uniprot_localization$Subcellular.location..CC.)] <- TRUE
uniprot_localization$ER[grep("Endoplasmic", uniprot_localization$Subcellular.location..CC.)] <- TRUE
uniprot_localization$ER[grep("endoplasmic", uniprot_localization$Subcellular.location..CC.)] <- TRUE
uniprot_localization$ER[grep("ER", uniprot_localization$Gene.ontology..cellular.component.)] <- TRUE
uniprot_localization$ER[grep("Endoplasmic", uniprot_localization$Gene.ontology..cellular.component.)] <- TRUE
uniprot_localization$ER[grep("endoplasmic", uniprot_localization$Gene.ontology..cellular.component.)] <- TRUE

length(which(uniprot_localization$ER == TRUE))

uniprot_loc_small <- uniprot_localization %>% select(Entry, ER, Subcellular.location..CC., Gene.ontology..cellular.component.) %>% 
  rename("PG.ProteinGroups" = "Entry",
         "ER_uniprot" = "ER",
         "Subcell_CC" = "Subcellular.location..CC.",
         "GOTerm_CC" = "Gene.ontology..cellular.component.")

merge_ER <- merge(wide, human_ER_uniprot, all.x = TRUE)
merge_ER <- merge(merge_ER, uniprot_loc_small, all.x = TRUE)

merge_ER$ER_proteinatlas[which(is.na(merge_ER$ER_proteinatlas))] <- FALSE
merge_ER$ER_uniprot[which(is.na(merge_ER$ER_uniprot))] <- FALSE

merge_ER$ER_all <- FALSE
merge_ER$ER_all[grep("TRUE", merge_ER$ER_proteinatlas)] <- TRUE
merge_ER$ER_all[grep("TRUE", merge_ER$ER_uniprot)] <- TRUE

merge_ER$Golgi <- FALSE
merge_ER$Golgi[grep("Golgi", merge_ER$Subcell_CC)] <- TRUE
merge_ER$Golgi[grep("golgi", merge_ER$Subcell_CC)] <- TRUE
merge_ER$Golgi[grep("Golgi", merge_ER$GOTerm_CC)] <- TRUE
merge_ER$Golgi[grep("golgi", merge_ER$GOTerm_CC)] <- TRUE

merge_ER$Secreted <- FALSE
merge_ER$Secreted[grep("Secreted", merge_ER$Subcell_CC)] <- TRUE
#merge_ER$Secreted[grep("Secreted", merge_ER$GOTerm_CC)] <- TRUE

merge_ER$ribosome <- FALSE
merge_ER$ribosome[grep("ribosomal protein", merge_ER$PG.ProteinDescriptions)] <- TRUE

length(which(merge_ER$ER_all != merge_ER$Golgi))

length(which(merge_ER$ER_all == TRUE))
length(which(merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE))
length(which(merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE | merge_ER$Secreted == TRUE))

length(which(merge_ER$ER_all == TRUE & 
               (merge_ER$p_adj_ATase1KO <= 0.05 & merge_ER$ATase1KO_WT <= -0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE) & 
               (merge_ER$p_adj_ATase1KO <= 0.05 & merge_ER$ATase1KO_WT <= -0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE | merge_ER$Secreted == TRUE) & 
               (merge_ER$p_adj_ATase1KO <= 0.05 & merge_ER$ATase1KO_WT <= -0.05)))

length(which(merge_ER$ER_all == TRUE & 
               (merge_ER$p_adj_ATase1KO <= 0.05 & merge_ER$ATase1KO_WT >= 0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE) & 
               (merge_ER$p_adj_ATase1KO <= 0.05 & merge_ER$ATase1KO_WT >= 0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE | merge_ER$Secreted == TRUE) & 
               (merge_ER$p_adj_ATase1KO <= 0.05 & merge_ER$ATase1KO_WT >= 0.05)))

length(which(merge_ER$ER_all == TRUE & 
               (merge_ER$p_adj_ATase2KO <= 0.05 & merge_ER$ATase2KO_WT <= -0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE) & 
               (merge_ER$p_adj_ATase2KO <= 0.05 & merge_ER$ATase2KO_WT <= -0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE | merge_ER$Secreted == TRUE) & 
               (merge_ER$p_adj_ATase2KO <= 0.05 & merge_ER$ATase2KO_WT <= -0.05)))

length(which(merge_ER$ER_all == TRUE & 
               (merge_ER$p_adj_ATase2KO <= 0.05 & merge_ER$ATase2KO_WT >= 0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE) & 
               (merge_ER$p_adj_ATase2KO <= 0.05 & merge_ER$ATase2KO_WT >= 0.05)))
length(which((merge_ER$Golgi == TRUE | merge_ER$ER_all == TRUE | merge_ER$Secreted == TRUE) & 
               (merge_ER$p_adj_ATase2KO <= 0.05 & merge_ER$ATase2KO_WT >= 0.05)))

#write.table(merge_ER, "20200117_merge_ER_all_sites.csv", sep = ",", row.names = FALSE)

secretory_pathway <- merge_ER %>% subset((ER_all == TRUE | Golgi == TRUE | Secreted == TRUE) &
                                           ((p_adj_ATase1KO <= 0.1 & ATase1KO_WT <= -0.05) | 
                                              (p_adj_ATase2KO <= 0.1 & ATase2KO_WT <= -0.05)))

length(which(merge_ER$p_adj_ATase1KO <= 0.05))
length(which(merge_ER$p_adj_ATase2KO <= 0.05))

ggplot(subset(merge_ER, ER_all == TRUE | Golgi == TRUE | Secreted == TRUE), 
       aes(ATase1KO_WT, -log10(p_adj_ATase1KO))) +
  geom_point(aes(color = (p_adj_ATase1KO <= 0.05) & ATase1KO_WT <= -0.05), show.legend = FALSE) +
  #facet_wrap(~ (ER_all == TRUE | Golgi == TRUE)) +
  geom_text_repel(data = subset(merge_ER, (p_adj_ATase1KO <= 0.05) & ATase1KO_WT <= -0.05 &
                                  ((ER_all == TRUE | Golgi == TRUE | Secreted == TRUE) & ribosome == FALSE)),
                  aes(ATase1KO_WT, -log10(p_adj_ATase1KO), label = PG.Genes)) +
  labs(y = "-log10(p value)") +
  theme_bw(base_size = 13)
#ggsave("ATase1KO_volcanoplot.pdf", height = 6, width = 6)

ggplot(subset(merge_ER, (ER_all == TRUE | Golgi == TRUE | Secreted == TRUE) & p_adj_ATase2KO != 0), 
       aes(ATase2KO_WT, -log10(p_adj_ATase2KO))) +
  geom_point(aes(color = (p_adj_ATase2KO <= 0.05) & ATase2KO_WT <= -0.05), show.legend = FALSE) +
  #facet_wrap(~ (ER_all == TRUE | Golgi == TRUE)) +
  geom_text_repel(data = subset(merge_ER, (p_adj_ATase2KO <= 0.05) & ATase2KO_WT <= -0.05 &
                                  ((ER_all == TRUE | Golgi == TRUE | Secreted == TRUE) & p_adj_ATase2KO != 0 & ribosome == FALSE)),
                  aes(ATase2KO_WT, -log10(p_adj_ATase2KO), label = PG.Genes)) +
  labs(y = "-log10(p value)") +
  theme_bw(base_size = 13)
#ggsave("ATase2KO_volcanoplot.pdf", height = 6, width = 6)

ggplot(subset(merge_ER, !is.na(ATase2KO_WT) & (ER_all == TRUE)), 
       aes(ATase2KO_WT >= 0)) +
  geom_bar(aes(fill = ATase2KO_WT >= 0)) +
  facet_wrap(~p_adj_ATase2KO <= 0.05, nrow = 2, scales = "free_y")

ggplot(subset(merge_ER, p_adj_ATase1KO <= 0.05), aes(ATase1KO_WT)) +
  geom_histogram(aes(fill = (ER_all == TRUE | Golgi == TRUE | Secreted == TRUE))) +
  labs(fill = "Secretory Pathway")

ggplot(subset(merge_ER, p_adj_ATase2KO <= 0.05), aes(ATase2KO_WT)) +
  geom_histogram(aes(fill = (ER_all == TRUE | Golgi == TRUE | Secreted == TRUE))) +
  labs(fill = "Secretory Pathway")



pathway_david <- merge_ER %>% 
  subset((p_adj_ATase1KO <= 0.05 & ATase1KO_WT <= -0.05) |
           (ATase2KO_WT <= -0.05 & p_adj_ATase2KO <= 0.05))

pathway_david2 <- merge_ER %>% 
  subset(((p_adj_ATase1KO <= 0.1 & ATase1KO_WT <= -0.05) |
            (ATase2KO_WT <= -0.05 & p_adj_ATase2KO <= 0.1)) & 
           (ER_all == TRUE | Golgi == TRUE))

#write.table(pathway_david, "pathway_david_sig_decreased.csv", sep = ",", row.names = FALSE)

