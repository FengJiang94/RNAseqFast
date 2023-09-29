## code to prepare `data` dataset goes here

built_in_length_mm10 <- read.csv("data-raw/built_in_length_mm10.csv", row.names = 1, header = T)
built_in_length_hg38 <- read.csv("data-raw/built_in_length_hg38.csv", row.names = 1, header = T)

built_in_lnc_ref_hg38 <- readRDS("data-raw/built_in_lnc_ref_hg38.RData")
built_in_lnc_ref_m39 <- readRDS("data-raw/built_in_lnc_ref_m39.RData")

# Data cleaning code here...
# (Do NOT put data analysis code here!)

# This should be the last line.
# If overwrite = T so everytime I run the script the
# updated objects are saved, but the default is overwrite = F
usethis::use_data(built_in_length_mm10, overwrite = T)
usethis::use_data(built_in_length_hg38, overwrite = T)
usethis::use_data(built_in_lnc_ref_hg38, overwrite = T)
usethis::use_data(built_in_lnc_ref_m39, overwrite = T)
