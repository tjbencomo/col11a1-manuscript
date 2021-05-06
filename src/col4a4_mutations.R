## Create heatmap similar to Figure 1C for COL4A4
## for reviewer request.

library(readr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)


get_region_element <- function(position) {
  if (position < 39) {
    return("Pre-7S Domain")
  } else if (position >= 39 & position <= 64) {
    return("7S Domain")
  } else if (position >= 65 & position <= 1459) {
    return("Triple-helical region")
  } else if (position >= 1460 & position <= 1464) {
    return("Post Triple Helix")
  } else {
    return("Collagen IV NC1")
  }
}

get_region <- Vectorize(get_region_element)


data_dir <- "data"
infp <- file.path(data_dir, "mutations.maf.gz")

muts <- read_tsv(infp)

coding <- c("Missense_Mutation", "Nonsense_Mutation")
col4a4 <- muts %>%
  filter(
    Hugo_Symbol == "COL4A4", 
    Variant_Classification %in% coding
  )

aachanges <- col4a4$HGVSp_Short
col4a4 <- col4a4 %>%
  mutate(
    refAA = str_extract(aachanges, "[A-Z]+"),
    position = as.numeric(str_extract(aachanges, "[0-9]+")),
    mutAA = str_extract(stringi::stri_reverse(aachanges), "[A-Z]+|\\*")
  ) %>%
  mutate(region = get_region(position))

region.lengths <- c(38, 64 - 39 + 1, 1459 - 65 + 1, 5, 1690 - 1465 + 1)
region.names <- c("Pre-7S Domain", "7S Domain", "Triple-helical region", 
                  "Post Triple Helix", "Collagen IV NC1")
region.info <- data.frame(
  region = region.names,
  size = region.lengths
)

df <- col4a4 %>%
  count(region, name = "n_mutations") %>%
  right_join(region.info) %>%
  mutate(n_mutations = tidyr::replace_na(n_mutations, 0)) %>%
  mutate(rate = n_mutations / size)

melted <- df %>%
  select(region, rate) %>%
  melt(id.var = "region")
  
p <- melted %>%
  inner_join(region.info) %>%
  mutate(region = factor(region, levels = region.names)) %>%
  ggplot(aes(variable, region)) + 
  geom_tile(aes(fill = value), width = .4) +
  coord_flip() +
  scale_fill_distiller(palette = "Reds", direction = 1, name = "Mutations/AA") +
  theme_classic() +
  labs(x = "", y = "Region")
print(p)
