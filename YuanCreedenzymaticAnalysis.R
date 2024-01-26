# Creedenzymatic Analysis

library(tidyverse)
library(creedenzymatic)

process_creedenzymatic <-
  function(krsa_path,uka_path, peptide_path) {
    krsa_data <- read_csv(krsa_path, show_col_types = FALSE) |>
      select(Kinase, Score = AvgZ) |>
      read_krsa(trns = "abs", sort = "desc")

    uka_data <- read_tsv(uka_path, show_col_types = FALSE) |>
     select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
    read_uka(trns = "abs", sort = "desc")

    peptide_data <-
      read_csv(peptide_path, show_col_types = FALSE) |>
      select(Peptide, Score = totalMeanLFC)

    kea3_data <-
      read_kea(
        peptide_data,
        sort = "asc",
        trns = "abs",
        method = "MeanRank",
        lib = "kinase-substrate"
      )

    ptmsea_data <-
      read_ptmsea(peptide_data)

    combined <- combine_tools(
      KRSA_df = krsa_data,
      #UKA_df = uka_data,
      KEA3_df = kea3_data,
      PTM_SEA_df = ptmsea_data
    )

    combined
  }

krsa_files <- c(
  "results/HPC-krsa_table_full_Exer_HPC_CTL_HPC_STK.csv",
  "results/HPC-krsa_table_full_HPC_Exer_HPC_CTL_PTK.csv",
  "results/STR-krsa_table_full_Exer_STR_CTL_STR_STK.csv",
  "results/STR-krsa_table_full_STR_Exer_STR_CTL_PTK.csv"
)

uka_files <- c(
"kinome_data/UKA-STK/HPC/Summaryresults 20230404-1530.txt",
"kinome_data/UKA-PTK/HPC/Summaryresults 20240126-1025.txt",
"kinome_data/UKA-STK/STR/Summaryresults 20230404-1531.txt",
"kinome_data/UKA-PTK/STR/Summaryresults 20240126-1027.txt"
)

peptide_files <- c(
  "results/HPC-dpp_Exer_HPC_CTL_HPC-STK.csv",
  "results/HPC-dpp_HPC_Exer_HPC_CTL-PTK.csv",
  "results/STR-dpp_Exer_STR_CTL_STR-STK.csv",
  "results/STR-dpp_STR_Exer_STR_CTL-PTK.csv"
)

result <-
  list(
    krsa_path = krsa_files,
    uka_path = uka_files,
    peptide_path = peptide_files
  ) |>
  pmap(process_creedenzymatic) |>
  set_names(c("HPCSTK", "HPCSTK", "STRSTK","STRPTK")) |>
  imap_dfr(~ write_csv(.x, str_glue("results/{.y}_EBYcreedenzymatic.csv")), .id = "Comparison")

