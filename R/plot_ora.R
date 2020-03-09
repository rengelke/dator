
# library(magrittr)
# library(autonomics)
# library(ggplot2)
#
#
# object <- dator:::data_sirt
#
#
# object %<>% dator::add_ora(db = "wiki")
# object %<>% dator::add_ora(db = "kegg")
#
# object %<>% add_ora(.,
#                   db = "kegg",
#                   cluster_type = "pbins",
#                   contrasts = NULL,
#                   sig_level = 0.05,
#                   p_cutoff = 0.05,
#                   q_cutoff = 0.2,
#                   simplify_go = TRUE,
#                   breaks = c(0.01, 0.05, 0.1, 0.2, 0.5),
#                   # breaks = c(0.05),
#                   k = 3,
#                   custom_db = NULL)


#' Plot over-representation analysis results
#'
#' @param object SummarizedExperiment object
#' @param plot_type plot of type: "barplot", "dotplot1", "dotplot2", "heatmap", "heatplot"
#' @param cluster_type over-representation analysis results from cluster type: "pval", "pbins", "bins"
#' @param contrasts contrasts, as defined in clustedefs of SummarizedExperiment object
#' @param db over-representation analysis results from annotation database: "wiki", "kegg", "go_bp", "go_mf", "go_cc", "msigH", "msigC1", "msigC2"
#' @param exclude exclude cluster (not in use at the moment)
#' @param sig_level p-value cutoff for heatplot
#' @param p_cutoff p-value cutoff for over-representation analysis results
#' @param q_cutoff p-value cutoff for over-representation analysis results
#'
#' @return Plot of over-representation analysis results
#' @export
#'
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#' @import ggplot2
#' @import autonomics
#'
#' @examples plot_ora(object, plot_type = "barplot", db = "wiki", cluster_type = "pval")
plot_ora <- function(object,
                     plot_type = NULL,
                     cluster_type = NULL,
                     contrasts = NULL,
                     db = NULL,
                     exclude = NULL,
                     sig_level = 0.05,
                     p_cutoff = 0.05,
                     q_cutoff = 0.5) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    assertive.types::assert_is_any_of(plot_type, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(cluster_type, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(contrasts, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(db, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(exclude, classes = c("NULL", "numeric"))
    assertive.types::assert_is_a_number(sig_level)
    assertive.types::assert_is_a_number(p_cutoff)
    assertive.types::assert_is_a_number(q_cutoff)

    plot_type <- match.arg(plot_type, c("barplot", "dotplot1", "dotplot2",
                                        "heatmap", "heatplot"), several.ok = FALSE)

    cluster_type <- match.arg(cluster_type, c("pval", "pbins", "bins"), several.ok = FALSE)
    db <- match.arg(db, c("wiki", "kegg", "go_bp", "go_mf", "go_cc",
                          "msigH", "msigC2", "msigC5", "custom"), several.ok = FALSE)

    if (is.null(contrasts)) {
        contrasts <- names(object@metadata$contrastdefs)[1]
    }

    message(paste0("Plotting with following parameters ...",
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



    if (plot_type == "barplot") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrasts]][[cluster_type]] %>%
            dplyr::arrange(pvalue) %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff & qvalue <= q_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup()
        if (nrow(tmp) == 0) {stop("No significant terms to plot.")}
        tmp %<>%
            dplyr::arrange(pvalue) %>%
            dplyr::arrange_at(dplyr::vars(dplyr::starts_with("clust_"))) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.))))

        p <- tmp %>%
            ggplot(aes(x = Description, y = logp,
                       fill = as.factor(!!rlang::sym(cluster_col[1])))) +
            geom_bar(stat = "identity") +
            facet_wrap(cluster_col[1],
                       nrow = 1,
                       labeller = as_labeller(label2)) +
            coord_flip() +
            labs(x = "Term", y = expression(-log[10]~p~value)) +
            geom_hline(yintercept = 1.3, linetype = "dashed",
                       color = "grey", size = 0.75) +
            scale_fill_manual(values = c("#5CA4A9", "#ED6A5A")) +
            #scale_fill_manual(values = c("#00B2CA", "#D1495B")) +
            #scale_fill_manual(values = c("#00B2CA", "#D96C06")) +
            theme_bw(base_size = 14) +
            theme(legend.position = "none")

        }



    if (plot_type == "dotplot1") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrasts]][[cluster_type]] %>%
            dplyr::arrange(pvalue) %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff & qvalue <= q_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup()
        if (nrow(tmp) == 0) {stop("No significant terms to plot.")}
        tmp %<>%
            dplyr::arrange(pvalue) %>%
            dplyr::arrange_at(dplyr::vars(dplyr::starts_with("clust_"))) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.)))) %>%
            dplyr::mutate(gene_ratio = GeneRatio %>%
                              sapply(., function (x) eval(parse(text=x)))) %>%
            dplyr::mutate_(cluster = cluster_col) %>%
            dplyr::mutate(cluster = factor(cluster)) %>%
            dplyr::mutate(logp = ifelse(logp >= 10, 10, logp))

        p <- tmp %>%
            ggplot(aes(x = cluster, y = Description, size = gene_ratio, color = logp)) +
            geom_point() +
            scale_x_discrete(labels = label2) +
            scale_color_continuous(low = "#85c1c5",  high = "#ED6A5A",
                                   space = "Lab", limit = c(0, 10),
                                   breaks = c(1.3, 3, 5, 8),
                                   labels = c("0.05", "0.001", "1e-05", "1e-08")) +
            labs(x = "", y = "Term", color = "p-value", size = "gene ratio") +
            theme_bw()

    }


    if (plot_type == "dotplot2") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrasts]][[cluster_type]] %>%
            dplyr::arrange(pvalue) %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff & qvalue <= q_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup()
        if (nrow(tmp) == 0) {stop("No significant terms to plot.")}
        tmp %<>%
            dplyr::mutate(gene_ratio = GeneRatio %>%
                              sapply(., function (x) eval(parse(text=x)))) %>%
            dplyr::mutate(logp = ifelse(logp >= 10, 10, logp)) %>%
            dplyr::arrange(dplyr::desc(gene_ratio)) %>%
            dplyr::arrange_at(dplyr::vars(dplyr::starts_with("clust_"))) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.))))

        p <- tmp %>%
            ggplot(aes(x = gene_ratio, y = Description, size = Count, color = logp)) +
            facet_wrap(colnames(tmp)[colnames(tmp) %>% grepl("clust_", .)],
                       nrow = 1,
                       labeller = as_labeller(label2)) +
            geom_point() +
            scale_color_gradient(low = "grey",  high = "red2",
                                   space = "Lab", limit = c(0, 10),
                                   breaks = c(1.3, 3, 5, 8),
                                   labels = c("0.05", "0.001", "1e-05", "1e-08")) +
            labs(x = "gene ratio", y = "Term", color = "p-value", size = "count") +
            theme_bw()

        }


        if (plot_type == "heatmap") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrasts]][[cluster_type]] %>%
            dplyr::arrange(pvalue) %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff & qvalue <= q_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup()
        if (nrow(tmp) == 0) {stop("No significant terms to plot.")}
        tmp %<>%
            dplyr::select("Description", "logp", cluster_col[1]) %>%
            dplyr::mutate(logp = ifelse(logp >= 10, 10, logp)) %>%
            reshape2::dcast(as.formula(paste("Description", cluster_col[1], sep="~")),
                            value.var = "logp") %>%
            replace(is.na(.), 0) %>%
            tibble::column_to_rownames(var = "Description")


        p <- tmp %>% pheatmap::pheatmap(cluster_cols = FALSE,
                                        labels_col = label2,
                                        angle_col = 0,
                                        color = colorRampPalette(
                                            RColorBrewer::brewer.pal(n = 7, name = "PuBu"))(100),
                                        silent = TRUE)
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

        genes2term <- ora_tmp %>% dplyr::select(Description, geneID) %>%
            apply(1, function (x) {
                x["geneID"] %>% stringr::str_split(pattern = "/", simplify = TRUE) %>%
                    t() %>% rbind.data.frame(stringsAsFactors = FALSE) %>%
                    dplyr::rename(geneID = V1) %>%
                    dplyr::mutate(Description = as.character(x["Description"]))
                }) %>%
            do.call(rbind, .)

        gene_oi <- ora_tmp %$% geneID %>% unlist() %>% paste(collapse = "/") %>%
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
            dplyr::select(Description, feature_name, effect) %>%
            dplyr::arrange(effect, feature_name) %>%
            dplyr::arrange(Description) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.)))) %>%
            dplyr::mutate(feature_name = feature_name %>%
                              factor(., levels = unique(.)))

        p <- tmp %>% ggplot(aes(x = feature_name, y = Description)) +
            geom_tile(aes(fill= effect)) +
            scale_colour_gradient2(low = scales::muted("red"), mid = "white",
                                   high = scales::muted("blue"), midpoint = 0, space = "Lab",
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



