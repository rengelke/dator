## code to prepare `DATASET` dataset goes here

# # set up data-raw folder containing DATSET.R
# usethis::use_data_raw(name = "DATASET")

library(magrittr)

# Wiki Pathways -----------------------------------------------------------

wiki_human <- clusterProfiler::read.gmt(
    rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", format = "gmt") )%>%
    tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")

wiki_mouse <- clusterProfiler::read.gmt(
    rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", format = "gmt") )%>%
    tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")

list.files() %>% .[grepl("wikipathways", .)] %>% lapply(file.remove)



# MSigDb ------------------------------------------------------------------

msigdb_H_human <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
msigdb_C2_human <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
msigdb_C5_human <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5")

msigdb_H_mouse <- msigdbr::msigdbr(species = "Mus musculus", category = "H")
msigdb_C2_mouse <- msigdbr::msigdbr(species = "Mus musculus", category = "C2")
msigdb_C5_mouse <- msigdbr::msigdbr(species = "Mus musculus", category = "C5")


# Color names -------------------------------------------------------------

color_names <- read.csv("~/../Dropbox/color_names.csv", stringsAsFactors = FALSE)




# Data sets ---------------------------------------------------------------

data_sirt <- readRDS("~/R_wd/dator/data-raw/data_sirt.rds")


# SAVE to R/sysdata.rda ---------------------------------------------------



usethis::use_data(wiki_human,
                  wiki_mouse,
                  msigdb_H_human,
                  msigdb_C2_human,
                  msigdb_C5_human,
                  msigdb_H_mouse,
                  msigdb_C2_mouse,
                  msigdb_C5_mouse,
                  color_names,
                  data_sirt,
                  internal = TRUE, overwrite = TRUE)
