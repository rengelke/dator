## code to prepare `DATASET` dataset goes here

# # set up data-raw folder containing DATSET.R
# usethis::use_data_raw(name = "DATASET")

library(magrittr)

# Wiki Pathways -----------------------------------------------------------

wiki_mouse <- clusterProfiler::read.gmt(
    rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", format = "gmt") )%>%
    tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")

wiki_human <- clusterProfiler::read.gmt(
    rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", format = "gmt") )%>%
    tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")

list.files() %>% .[grepl("wikipathways", .)] %>% lapply(file.remove)



# Color names -------------------------------------------------------------

color_names <- read.csv("~/../Dropbox/color_names.csv", stringsAsFactors = FALSE)




# Data sets ---------------------------------------------------------------

data_sirt <- readRDS("~/R_wd/dator/data-raw/data_sirt.rds")


# SAVE to R/sysdata.rda ---------------------------------------------------



usethis::use_data(wiki_human,
                  wiki_mouse,
                  color_names,
                  data_sirt,
                  internal = TRUE, overwrite = TRUE)
