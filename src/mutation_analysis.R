# File: mutation_analysis.R
# Author: Tomas Bencomo
# Description: Calculate mutation rates for each gene used to
# rank genes in Figure 1A.
# The exon_info.csv file was generated by querying bioMart for
# all human ensembl transcripts with a cds_length > 0.



library(readr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(tximport)
library(annotables)
library(stringr)

# Load differentiation RNASeq data that measures expression levels
# of differentiating skin at days 0, 3, and 6. Used to filter out
# genes with low expression levels
get_rnaseq_data <- function() {
  files <- c('d0_403879.genes.results', 'd3_403880.genes.results', 'd6_403881.genes.results')
  files <- file.path("data/rnaseq", files)
  samples <- c('day0', 'day3', 'day6')
  names(files) <- samples
  txi.rsem <- tximport(files, type = 'rsem', txIn=FALSE, txOut = FALSE)
  txi.rsem$length[txi.rsem$length == 0] <- 1
  condition <- factor(c('d0', 'd3', 'd6'))
  sample_info <- data.frame(sample = samples,
                            condition = condition)
  dds <- DESeqDataSetFromTximport(
    txi.rsem, 
    colData=sample_info, 
    design = ~ condition
  )
  dds <- estimateSizeFactors(dds)
  mat <- data.frame(counts(dds, normalized=TRUE))
  mat$ensembl_id <- rownames(mat)
  mat$ensembl_id <- str_extract(mat$ensembl_id, "[A-Z]+[0-9]+")
  mat <- mat %>%
    left_join(grch38, by = c("ensembl_id" = "ensgene")) %>%
    drop_na(ensembl_id) %>%
    select(symbol, ensembl_id, day0, day3, day6) %>%
    dplyr::rename(Hugo_Symbol = symbol)
}

# Includes CDS length for relevant transcripts and total exon length
# (length of all combined exons for the gene). This info is used
# to create missense mutations/bp column. CDS length is ultimately 
# used to normalize number of missense mutations by gene
get_gene_data <- function(mutations) {
  exons_cds <- read_csv("data/exon_info.csv") #includes CDS length
  exons_total <- read_csv("data/exon_lengths.csv", 
                          col_names = c("Hugo_Symbol", "TotalExonLength"))
  gene_info <- mutations %>%
    select(Hugo_Symbol, Transcript_ID) %>%
    distinct() %>%
    inner_join(exons_cds, by = c("Transcript_ID" = "ensembl_transcript_id")) %>%
    inner_join(exons_total)
}

mutation_types <- c("Missense_Mutation", "Splice_Site", 
                    "Nonsense_Mutation", "Silent")
patients <- read_csv("data/patient_metadata.csv", col_types = c(patient = "f"))
mutations <- read_delim("data/mutations.maf.gz", 
                        "\t", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = c(
                          t_ref_count = "c", t_alt_count = "c", t_depth = "c", 
                          n_depth = "c", n_alt_count = "c", n_ref_count = "c"
                        ))

gene_info <- get_gene_data(mutations)
rnaseq <- get_rnaseq_data()

# Compute Missense to Silent Ratio and missense to bases ratios
snp_rates <- mutations %>%
  filter(Variant_Classification %in% mutation_types) %>%
  mutate(Variant_Classification = factor(Variant_Classification, levels=mutation_types)) %>%
  dplyr::count(Hugo_Symbol, Variant_Classification, .drop = F) %>%
  spread(Variant_Classification, n) %>%
  inner_join(gene_info) %>%
  mutate(
    missense.vs.silent = Missense_Mutation / Silent,
    mutations_norm_cds = Missense_Mutation / cds_length,
    mutations_norm_total = Missense_Mutation / TotalExonLength
  )

# Compute proportion of patients with important mutation for each gene
sample_rates <- mutations %>%
  mutate(patient = factor(patient)) %>%
  filter(
    Variant_Classification %in% c("Missense_Mutation", "Splice_Site", 
                                  "Nonsense_Mutation")) %>%
  mutate(
    Variant_Classification = factor(Variant_Classification, 
                                    levels = c("Missense_Mutation", "Splice_Site", 
                                               "Nonsense_Mutation")
                                    )) %>%
  dplyr::count(Hugo_Symbol, patient, .drop = FALSE) %>%
  mutate(mutated = case_when(
    n > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  group_by(Hugo_Symbol) %>%
  summarise(n = sum(mutated),
            proportion = n / 100)

# Join gene information, proportion of SCCs mutated, and missense ratios
gene_rates <- sample_rates %>%
  inner_join(snp_rates) %>%
  dplyr::rename(samples_mutated = n) %>%
  select(Hugo_Symbol, samples_mutated, proportion, missense.vs.silent, 
         mutations_norm_cds, mutations_norm_total, Missense_Mutation,
         Silent, Nonsense_Mutation, Splice_Site, Transcript_ID,
         cds_length, TotalExonLength) %>%
  filter(Silent > 0) %>%
  left_join(rnaseq) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)


# Use thresholds to prioritize genes with top biological relevance
expr_thresh <- 20
scc_thresh <- 15
miss.vs.silent_thresh <- 5

res <- gene_rates %>%
  filter(day0 > expr_thresh, day3 > expr_thresh, day6 > expr_thresh,
         samples_mutated > scc_thresh, missense.vs.silent >= 5
         ) %>%
  arrange(desc(mutations_norm_cds))
  
## Add HRAS and NOTCH1 for reference
final_res <- rbind(head(res, 7), filter(gene_rates, Hugo_Symbol %in% c("HRAS", "NOTCH1")))
print(final_res)


