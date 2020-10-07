



#' Plot features from over-representation analysis results
#'
#' @param object SummarizedExperiment object
#' @param plot_type plot of type: "violinplot", "heatplot"
#' @param cluster_type over-representation analysis results from cluster type: "pval", "pbins", "bins"
#' @param contrasts contrasts, as defined in clustedefs of SummarizedExperiment object
#' @param db over-representation analysis results from annotation database: "wiki", "kegg", "go_bp", "go_mf", "go_cc", "msigH", "msigC1", "msigC2"
#' @param exclude exclude cluster (not in use at the moment)
#' @param sig_level feature p-value cutoff for heatplot
#' @param p_cutoff p-value cutoff for over-representation analysis results
#' @param q_cutoff p-value cutoff for over-representation analysis results
#'
#' @return Plot features of over-representation analysis results
#' @export
#'
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#' @import ggplot2
#' @import autonomics
#'
#' @examples plot_ora_features(object, plot_type = "heatplot", db = "wiki", cluster_type = "pval")
plot_ora_features <- function(object,
                              plot_type = NULL,
                              cluster_type = NULL,
                              contrasts = NULL,
                              db = NULL,
                              exclude = NULL,
                              select_term_id = NULL,
                              sig_level = 0.05,
                              p_cutoff = 0.05,
                              q_cutoff = 0.5) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    assertive.types::assert_is_any_of(plot_type, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(cluster_type, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(contrasts, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(db, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(exclude, classes = c("NULL", "numeric"))
    assertive.types::assert_is_any_of(select_term_id, classes = c("NULL", "character"))
    assertive.types::assert_is_a_number(sig_level)
    assertive.types::assert_is_a_number(p_cutoff)
    assertive.types::assert_is_a_number(q_cutoff)

    plot_type <- match.arg(plot_type, c("violinplot", "heatplot"), several.ok = FALSE)

    cluster_type <- match.arg(cluster_type, c("pval", "pbins", "bins", "custom"), several.ok = FALSE)
    db <- match.arg(db, c("wiki", "kegg", "go_bp", "go_mf", "go_cc",
                          "msigH", "msigC2", "msigC5", "custom"), several.ok = FALSE)

    if (is.null(contrasts)) {
        contrasts <- names(object@metadata$contrastdefs)[1]
    }

    message(paste0("Creating ", plot_type,  " with following parameters ...",
                   "\ncontrast: ", contrasts,
                   "\ncluster type: ", cluster_type,
                   "\ndatabase: ", db))


    if (any(grepl(names(object@metadata), pattern = "ora_")) == FALSE) {
        stop("Object does not contain over-representation analysis results")
    }
    if (any(grepl(names(object@metadata), pattern = paste0("ora_", db))) == FALSE) {
        stop(paste0("Object does not contain ", db, " over-representation analysis results"))
    }
    cluster_col <- c(paste0("clust_", cluster_type), paste0("clust_", cluster_type, "_label"))
    if (any(grepl(colnames(object@metadata[[paste0("ora_", db)]][[contrasts]][[cluster_type]]),
                  pattern = cluster_col[1])) == FALSE) {
        stop("Can't find cluster definition column in over-representation analysis results")
    }


    n_clust <- object@metadata$clusterdefs[[contrasts]][, cluster_col[1]] %>%
        table()

    label1 <- paste0("C", n_clust %>% names(), "\n(n = ", n_clust, ")") %>%
        `names<-`(n_clust %>% names())
    label2 <- object@metadata$clusterdefs[[contrasts]][, cluster_col[2]] %>%
        levels(.) %>%
        paste0(., "\n(n = ", n_clust, ")") %>%
        `names<-`(n_clust %>% names())

    palette <- c("#5ca4a9", "#ed6a5a", "#8fbc62", "#AA7484",
                 "#f49958", "#a994d1", "#aaaaaa")[1:length(n_clust)]



    if (plot_type == "violinplot") {

        ## supposed to work only for cluster_type == "pval"
        if (cluster_type != "pval") {
            message("Heatplot only available for 'pval' cluster type.
                    \nConverting cluster type to 'pval'")
            cluster_type <- "pval"
        }

        ora_tmp <- object@metadata[[paste0("ora_", db)]][[contrasts]][[cluster_type]] %>%
            dplyr::filter(pvalue <= p_cutoff & qvalue <= q_cutoff)

        genes2term <- ora_tmp %>%
            dplyr::select(ID, Description, geneID) %>%
            apply(1, function (x) {
                x["geneID"] %>% stringr::str_split(pattern = "/", simplify = TRUE) %>%
                    t() %>% rbind.data.frame(stringsAsFactors = FALSE) %>%
                    dplyr::rename(geneID = V1) %>%
                    dplyr::mutate(ID = as.character(x["ID"]),
                                  Description = as.character(x["Description"]))
            }) %>%
            do.call(rbind, .)

        if (!is.null(select_term_id)) {
            genes2term %<>% dplyr::filter(ID %in% select_term_id)
        }

        gene_oi <- genes2term %$% geneID %>% unlist() %>% paste(collapse = "/") %>%
            stringr::str_split(pattern = "/") %>% unlist()

        gene_all <- fdata(object)$ENTREZID %>%
            stringr::str_split(pattern = ";", simplify = TRUE) %>% .[, 1]

        gene_oi_idx <- gene_all %in% gene_oi

        genes2effect <- object@metadata$limma[, contrasts, ] %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "feature_id") %>%
            dplyr::mutate(geneID = gene_all) %>%
            dplyr::filter(geneID %in% gene_oi) %>%
            dplyr::mutate(feature_name = fdata(object)$feature_name[gene_oi_idx]) %>%
            dplyr::select(geneID, feature_id, feature_name, p, effect)

        tmp <- dplyr::left_join(genes2term, genes2effect, by = "geneID") %>%
            dplyr::filter(p <= sig_level) %>%
            dplyr::mutate(logp = -1*log(p, 10)) %>%
            dplyr::select(Description, feature_name, effect, logp) %>%
            dplyr::arrange(effect, feature_name) %>%
            dplyr::arrange(Description) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.)))) %>%
            dplyr::mutate(feature_name = feature_name %>%
                              factor(., levels = unique(.)))

        p <- tmp %>%
            ggplot(aes(x = Description, y = effect, color = logp)) +
            geom_violin(width = 1.3) +
            coord_flip() +
            geom_jitter(alpha = 0.2, width = .02) +
            scale_color_gradientn(colours = c('#e8cea6', '#e0bc8e', '#d7ab7c', '#cd996d', '#c38760', '#b97655', '#af644a', '#a45240', '#994037')) +
            #ggbeeswarm::geom_quasirandom(alpha = 0.2, width = 0.2) +
            theme_bw()

    }


    if (plot_type == "heatplot") {

        ## supposed to work only for cluster_type == "pval"
        if (cluster_type != "pval") {
            message("Heatplot only available for 'pval' cluster type.
                    \nConverting cluster type to 'pval'")
            cluster_type <- "pval"
        }


        ora_tmp <- object@metadata[[paste0("ora_", db)]][[contrasts]][[cluster_type]] %>%
            dplyr::filter(pvalue <= p_cutoff & qvalue <= q_cutoff)

        genes2term <- ora_tmp %>%
            dplyr::select(ID, Description, geneID) %>%
            apply(1, function (x) {
                x["geneID"] %>% stringr::str_split(pattern = "/", simplify = TRUE) %>%
                    t() %>% rbind.data.frame(stringsAsFactors = FALSE) %>%
                    dplyr::rename(geneID = V1) %>%
                    dplyr::mutate(ID = as.character(x["ID"]),
                                  Description = as.character(x["Description"]))
            }) %>%
            do.call(rbind, .)

        if (!is.null(select_term_id)) {
            genes2term %<>% dplyr::filter(ID %in% select_term_id)
        }

        gene_oi <- genes2term %$% geneID %>% unlist() %>% paste(collapse = "/") %>%
            stringr::str_split(pattern = "/") %>% unlist()

        gene_all <- fdata(object)$ENTREZID %>%
            stringr::str_split(pattern = ";", simplify = TRUE) %>% .[, 1]

        gene_oi_idx <- gene_all %in% gene_oi

        genes2effect <- object@metadata$limma[, contrasts, ] %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "feature_id") %>%
            dplyr::mutate(geneID = gene_all) %>%
            dplyr::filter(geneID %in% gene_oi) %>%
            dplyr::mutate(feature_name = fdata(object)$feature_name[gene_oi_idx]) %>%
            dplyr::select(geneID, feature_id, feature_name, p, effect)


        tmp <- dplyr::left_join(genes2term, genes2effect, by = "geneID") %>%
            dplyr::filter(p <= sig_level) %>%
            dplyr::mutate(logp = -1*log(p, 10)) %>%
            dplyr::select(Description, feature_name, effect, logp) %>%
            dplyr::arrange(effect, feature_name) %>%
            dplyr::arrange(Description) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.)))) %>%
            dplyr::mutate(feature_name = feature_name %>%
                              factor(., levels = unique(.)))

        p <- tmp %>% ggplot(aes(x = feature_name, y = Description)) +
            geom_tile(aes(fill= effect)) +
            scale_colour_gradient2(low = "#1c807d", mid = "white",
                                   high = "#B7412D", midpoint = 0, space = "Lab",
                                   na.value = "white", guide = "colourbar", aesthetics = "fill") +
            #theme_bw() +
            theme(panel.grid.major = element_blank(),
                axis.text.x = element_text(angle = 60, hjust = 1)) +
            labs(x = NULL, y = NULL, legend = "effect") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

p

}



