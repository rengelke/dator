#' Set cluster definitions
#'
#' @param object SummarizedExperiment object
#' @param cluster_type type of cluster to generate: "pval", "pval_bins", "bins"
#' @param contrastdefs contrast definition(s) for which cluster should be generated
#' @param breaks numeric vector indicating p-value cutoffs for "pval_bins"
#' @param sig_level significance cutoff for p-value based cluster "pval"
#' @param k number of bins for bin cluster "bins"
#'
#' @return updated SummarizedExperiment object
#' @export
#'
#' @import magrittr
#' @import dplyr
#' @import autonomics
#'
#' @examples add_cluster(object)
add_cluster <- function (object,
                         cluster_type = c("pval", "pval_bins", "bins"),
                         contrastdefs = NULL,
                         breaks = c(0.05, 0.10, 0.20),
                         sig_level = 0.05,
                         k = 5) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    if (assertthat::has_name(object@metadata, c("contrastdefs", "limma")) == FALSE) {
        stop("Object metadata does not contain limma results")
    }

    assertive.types::assert_is_any_of(cluster_type, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(contrastdefs, classes = c("NULL", "character"))
    sig_level %>% assertive.types::assert_is_a_number()
    k %>% assertive.types::assert_is_a_number() %>%
        assertive.numbers::assert_all_are_whole_numbers() %>%
        assertive.numbers::assert_all_are_greater_than(0)

    cluster_type <- match.arg(cluster_type, c("pval", "pval_bins", "bins"), several.ok = TRUE)

    if (is.null(contrastdefs)) {
        contrastdefs <- names(object@metadata$contrastdefs)
    }


    clusterdefs_list <- contrastdefs %>% lapply(function (x) {

        clusterdefs <- object@metadata$limma[, x, ] %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "feature_id")

        if ("pval" %in% cluster_type) {
            clusterdefs %<>%
                dplyr::mutate(clust_pval =  ifelse(p <= sig_level & effect <= 0, 1,
                                                   ifelse(p <= sig_level & effect >= 0, 2,
                                                          3)))

        }

        if ("pval_bins" %in% cluster_type) {

            clusterdefs %<>%
                dplyr::mutate(clust_pbins = ifelse(p <= 0.05 & effect <= 0, 1,
                                                   ifelse(p <= 0.10 & p > 0.05 & effect <= 0, 2,
                                                          ifelse(p <= 0.20 & p > 0.10 & effect <= 0, 3,
                                                                 ifelse(p <= 1 & p > 0.20, 4,
                                                                        ifelse(p <= 0.20 & p > 0.10 & effect >= 0, 5,
                                                                               ifelse(p <= 0.10 & p > 0.05 & effect >= 0, 6,
                                                                                      ifelse(p <= 0.05 & effect >= 0, 7,
                                                                                             NA))))))))
        }

        if ("bins" %in% cluster_type) {

            clusterdefs %<>%
                dplyr::mutate(effect_tmp = ifelse(!is.na(p), effect, NA)) %>%
                dplyr::mutate(clust_bins = cut(effect_tmp,
                                               quantile(effect_tmp, (0:k)/k, na.rm = TRUE),
                                               include.lowest = TRUE, labels = FALSE)) %>%
                dplyr::select(-effect_tmp)
        }

        clusterdefs %<>%
            tibble::column_to_rownames(var = "feature_id") %>%
            as.matrix()
    })

    names(clusterdefs_list) <- contrastdefs

    if (!all(rownames(clusterdefs_list[[1]]) == fdata(object)$feature_id)) {
        stop("Not able to align feature names of clusterdefs.")
    }

    S4Vectors::metadata(object)$clusterdefs <- array(dim = c(nrow(clusterdefs_list[[1]]),
                                                             length(contrastdefs),
                                                             ncol(clusterdefs_list[[1]])),
                                                  dimnames = list(feature = rownames(clusterdefs_list[[1]]),
                                                                  contrast = contrastdefs,
                                                                  quantity = colnames(clusterdefs_list[[1]])))

    for (i in seq_along(clusterdefs_list)) {
        S4Vectors::metadata(object)$clusterdefs[, names(clusterdefs_list)[i], ] <- clusterdefs_list[[i]]
    }

    object

}

