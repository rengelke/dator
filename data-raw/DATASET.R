## code to prepare `DATASET` dataset goes here

# # set up data-raw folder containing DATSET.R
# usethis::use_data_raw(name = "DATASET")

library(magrittr)


# UniProt -----------------------------------------------------------------

uniprot_human <- UniProt.ws::UniProt.ws(taxId = 9606)
uniprot_mouse <- UniProt.ws::UniProt.ws(taxId = 10090)


# Wiki Pathways -----------------------------------------------------------

wiki_human <- clusterProfiler::read.gmt(
    rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", format = "gmt") )%>%
    tidyr::separate(ont, c("name", "version", "wpid", "org"), "%") %>%
    dplyr::rename(id = wpid)

cross_ref <- UniProt.ws::select(uniprot_human,
                   keys = wiki_human$gene,
                   keytype = "ENTREZ_GENE",
                   columns = c("UNIPROTKB"))

wiki_human %<>% dplyr::right_join(., cross_ref, by = c("gene" = "ENTREZ_GENE")) %>%
    dplyr::rename(uniprot = UNIPROTKB) %>%
    dplyr::select(id, name, uniprot, gene)


wiki_mouse <- clusterProfiler::read.gmt(
    rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", format = "gmt") )%>%
    tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")




# MSigDb ------------------------------------------------------------------

msigdb_H_human <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
msigdb_C2_human <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
msigdb_C5_human <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5")

msigdb_H_mouse <- msigdbr::msigdbr(species = "Mus musculus", category = "H")
msigdb_C2_mouse <- msigdbr::msigdbr(species = "Mus musculus", category = "C2")
msigdb_C5_mouse <- msigdbr::msigdbr(species = "Mus musculus", category = "C5")


# Reactome ----------------------------------------------------------------

download.file("https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt",
              "reactome.txt")
reactome <- readr::read_delim("reactome.txt", "\t", col_names = FALSE,
                              escape_double = FALSE, trim_ws = TRUE) %>%
    `colnames<-`(c("uniprot", "id", "link", "name", "type", "organism")) %>%
    dplyr::mutate(gene = NA)

reactome_human <- reactome %>% dplyr::filter(organism == "Homo sapiens") %>%
    dplyr::select(id, name, uniprot, gene)
reactome_mouse <- reactome %>% dplyr::filter(organism == "Mus musculus") %>%
    dplyr::select(id, name, uniprot, gene)


# Color names -------------------------------------------------------------

#color_names <- read.csv("~/../Dropbox/color_names.csv", stringsAsFactors = FALSE)




# Data sets ---------------------------------------------------------------

data_sirt <- readRDS("~/R_wd/dator/data-raw/data_sirt.rds")
data_dp <- readRDS("~/R_wd/dator/data-raw/data_dp.rds")


# Clean -------------------------------------------------------------------

list.files() %>% .[grepl("wikipathways|reactome", .)] %>%
    lapply(file.remove)


# SAVE to R/sysdata.rda ---------------------------------------------------



usethis::use_data(wiki_human,
                  wiki_mouse,
                  msigdb_H_human,
                  msigdb_C2_human,
                  msigdb_C5_human,
                  msigdb_H_mouse,
                  msigdb_C2_mouse,
                  msigdb_C5_mouse,
                  reactome_human,
                  reactome_mouse,
                  uniprot_human,
                  uniprot_mouse,
                  #color_names,
                  data_sirt,
                  data_dp,
                  internal = TRUE,
                  overwrite = TRUE)

