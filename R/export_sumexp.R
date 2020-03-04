

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

    limma_res <- object %>%
       autonomics.import::extract_limma_dt() %>%
       .[order(.$p), ] %>%
       `colnames<-`(c("ID", "Gene", "contrast", "effect [log ratio]",
                     "p-value", "FDR", "Bonf"))

    raw_data <- autonomics.import::exprs(object) %>%
       as.data.frame() %>%
       cbind(Gene = autonomics.import::fdata(object)$feature_name,
             Uniprot = autonomics.import::fdata(object)$feature_uniprot,
             .) %>%
       tibble::rownames_to_column(var = "ID")

   xlsx::write.xlsx(limma_res, file = file, sheetName = "statistics", row.names = FALSE)
   xlsx::write.xlsx(raw_data, file = file, sheetName = "data", append = TRUE, row.names = FALSE)
}
