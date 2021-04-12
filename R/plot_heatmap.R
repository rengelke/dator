

#' Plot heatmap
#'
#' @param object SummarizedExperiment object
#' @param contrasts  contrasts, as defined in clustedefs of SummarizedExperiment object
#' @param p_sig_level p-value cutoff for features to be used for plotting
#' @param fdr_sig_level fdr cutoff for features to be used for plotting
#' @param name heatmap legend name
#' @param scale logical, scale expession values if TRUE
#' @param center logical, center expression values if TRUE
#' @param row_km number of row clusters 'slices' to be generated
#' @param breaks heatmap scale midpoint and limits
#' @param colors heatmap colors; 3 values for low, midpoint and high
#' @param ... further parameters passed to ComplexHeatmap::Heatmap function
#'
#' @return Heatmap containing expression values
#'
#' @importFrom magrittr %<>% %>%
#' @import autonomics
#'
#' @export
#'
#' @examples plot_heatmap(object, contrasts = "Sirt_vs_Mock", p_cutoff = 0.001, scale = TRUE, cluster_columns = FALSE, row_km = 3, name = c("log ratio"), show_row_names = FALSE, border = TRUE, breaks = c(-1.5, 0, 1.5), colors = c('#046D68', "#E8E8E8", '#B13407'))
plot_heatmap <- function (object,
                          contrasts = NULL,
                          p_sig_level = 0.05,
                          fdr_sig_level = 1,
                          name = " ",
                          scale = "none",
                          row_km = NULL,
                          breaks = c(-2.5, 0, 2.5),
                          colors = c("#306866", "#F8F2E2", "#E54219"),
                          ...) {

    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    if (assertthat::has_name(object@metadata, c("contrastdefs", "limma")) == FALSE) {
        stop("Object metadata does not contain limma results")
    }
    assertive.types::assert_is_any_of(contrasts, classes = c("NULL", "character"))


    if (is.null(contrasts)) {
        contrasts <- names(object@metadata$contrastdefs)
    }

    scale <- match.arg(scale, c("none", "row", "column"))

    message(paste0("Creating heatmap with following parameters ...",
                   "\ncontrast: ", paste(contrasts, collapse = " / "),
                   "\np-value cutoff: ", p_sig_level,
                   "\nfdr cutoff: ", fdr_sig_level,
                   "\nscale: ", scale))


    idx1 <- object@metadata$limma[, contrasts, "p"] %>%
        apply(1, function (x) sum(x <= p_sig_level, na.rm = TRUE) >= 1)

    idx2 <-  object@metadata$limma[, contrasts, "fdr"] %>%
        apply(1, function (x) sum(x <= fdr_sig_level, na.rm = TRUE) >= 1)


    idx <- idx1 & idx2 & !is.na(idx1) & !is.na(idx2)

    mat <- exprs(object) %>%
        `rownames<-`(fdata(object)$feature_name %>% make.names(unique = TRUE)) %>%
        as.matrix()
    mat <- mat[idx, ]

    if (scale != "none") {
        mat <- switch(scale,
                      none = mat,
                      row = dator::scale_rows(mat),
                      column = t(dator::scale_rows(t(mat))))
    }

    mat %<>% replace(is.na(.), 0)


    col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors, space = "LUV")

    mat %>% ComplexHeatmap::Heatmap(name = name, row_km = row_km, col = col_fun, ...) %>%
        ComplexHeatmap::draw() -> hmap


    if (!is.null(row_km)) {
        custom_clust <- ComplexHeatmap::row_order(hmap) %>%
            plyr::ldply(., cbind, .id = "clust_custom") %>%
            .[order(.[, 2]), ] %>%
            dplyr::mutate(feature_id = rownames(mat)) %>%
            dplyr::select(feature_id, clust_custom)

        hmap@layout$custom_clust <- custom_clust
        message("Cluster definition stored in heatmap object (`hmap@layout$custom_clust`)")
    }

    hmap

}





