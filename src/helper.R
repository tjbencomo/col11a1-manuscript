library(dplyr)
library(stringr)

## Turn HGVSp_Short column into separate columns with the
## reference and mutant amino acids and AA position
parse_HGVSp <- function(mutations) {
  stopifnot("HGVSp_Short" %in% colnames(mutations))
  s <- str_sub(mutations$HGVSp_Short, 3)
  aa <- str_split(s, "[0-9]+", simplify = T)
  mutations <- mutations %>%
    mutate(
      refAA = aa[, 1],
      position = as.numeric(str_extract(HGVSp_Short, "[0-9]+")),
      mutAA = aa[, 2]
    ) %>%
    mutate(
      mutAA = case_when(
        str_detect(mutAA, "_splice") ~ "splice",
        str_detect(HGVSp_Short, "=") ~ refAA,
        TRUE ~ mutAA
      ))
  return(mutations)
}

## Anotate whether a coding mutation affects a Proline or
## Glycine residue
annotate_pg_mutations <- function(mutations) {
  stopifnot(c("refAA", "position", "mutAA") %in% colnames(mutations))
  mutations <- mutations %>%
    mutate(pg_mutation = case_when(
      refAA %in% c("P", "G") & !(mutAA %in%c("P", "G")) ~ 1,
      TRUE ~ 0
    ))
  return(mutations)
}

## Determine which region of the COL11A1 protein a position
## falls in
get_col11a1_regions_element <- function(position) {
  if (is.na(position)) {
    return(NA)
  }
  if (position < 230) {
    return("Pre-Nonhelical region")
  } else if (position >= 230 && position <= 419) {
    return("Nonhelical region")
  } else if (position >= 420 && position <= 508) {
    return("Triple-helical region (interrupted)")
    # return("Triple-helical region")
  } else if (position >= 509 && position <= 511) {
    return("Short nonhelical segment")
  } else if (position >= 512 && position <= 528) {
    return("Telopeptide")
  } else if (position >= 529 && position <= 1542) {
    return("Triple-helical region")
  } else if (position >= 1543 && position <= 1563) {
    return("Nonhelical region (C-terminal)")
  } else {
    return("Fibrillar collagen NC1")
  }
}
get_col11a1_regions <- Vectorize(get_col11a1_regions_element)

## Determines whether an AA position for a given collagen
## falls in the triple helix region as specified by UniProt
isInHelix <- function(gene, position) {
  if(is.na(position)) {
    return(NA)
  }
  if (gene == "COL2A1" & position>= 201 & position <= 1214) {
    return(TRUE)
  } else if (gene == "COL11A2" & position >= 487 & position <= 1500) {
    return(TRUE)
  } else if (gene == "COL5A1" & position >= 559 & position <= 1570) {
    return(TRUE)
  } else if (gene == "COL11A1" & ((position >= 420 & position <= 508) | (position >= 529 & position <= 1542))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
