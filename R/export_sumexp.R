

#' Export SummarizedExperiment to Excel file
#'
#' @param object SummarizedExperiment object
#' @param file file name
#'
#' @return
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @examples export_sumexp(object, file = "./results/study_results.xlsx")
export_sumexp <- function (object, file) {

    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))

    n_contrast <- length(object@metadata$contrastdefs)

    uniprot_ac <- fdata(object) %>%
        dplyr::select(feature_id, feature_uniprot)

    if (n_contrast == 1) {
        limma_tmp <- object %>%
            autonomics.import::extract_limma_dt()
        limma_tmp %<>% dplyr::mutate(contrast = names(object@metadata$contrastdefs))
        limma_res <- limma_tmp %>%
            dplyr::left_join(., uniprot_ac, by = c("fid" = "feature_id")) %>%
            dplyr::select(contrast, fid, feature_uniprot, fname,
                          effect, p, fdr, bonf) %>%
            `colnames<-`(c("Contrast", "ID", "Uniprot AC", "Gene", "Effect [log ratio]",
                           "P", "FDR", "Bonf")) %>%
            dplyr::mutate(" " = " ")

    } else {
        seq_len(n_contrast) %>% lapply(function (x) {
            limma_tmp <- object %>%
                autonomics.import::extract_limma_dt() %>%
                dplyr::filter(contrast == names(object@metadata$contrastdefs[x]))
            limma_tmp %>% dplyr::left_join(., uniprot_ac,
                                           by = c("fid" = "feature_id")) %>%
                dplyr::select(contrast, fid, feature_uniprot, fname,
                              effect, p, fdr, bonf) %>%
                `colnames<-`(c("Contrast", "ID", "Uniprot AC", "Gene", "Effect [log ratio]",
                               "P", "FDR", "Bonf")) %>%
                dplyr::mutate(" " = " ")
        }) %>%
            do.call(cbind, .) %>%
            `<<-`(limma_res, .)
    }


    raw_data <- autonomics.import::exprs(object) %>%
        as.data.frame() %>%
        cbind(Gene = autonomics.import::fdata(object)$feature_name,
              Uniprot = autonomics.import::fdata(object)$feature_uniprot,
              .) %>%
        tibble::rownames_to_column(var = "ID")

    xlsx::write.xlsx(limma_res, file = file, sheetName = "statistics", row.names = FALSE)
    xlsx::write.xlsx(raw_data, file = file, sheetName = "data", append = TRUE, row.names = FALSE)
}



