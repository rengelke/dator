
library(magrittr)
library(autonomics)
library(ggplot2)


object <- dator:::data_sirt
object %<>% dator::add_ora(db = "wiki")
object %<>% dator::add_ora(db = "kegg")




plot_ora <- function(object,
                     plot_type = NULL,
                     cluster_type = NULL,
                     contrastdefs = NULL,
                     db = NULL,
                     exclude = NULL,
                     sig_level = 0.05,
                     p_cutoff = 0.05,
                     q_cutoff = 0.05) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    assertive.types::assert_is_any_of(plot_type, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(cluster_type, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(contrastdefs, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(db, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(exclude, classes = c("NULL", "numeric"))
    assertive.types::assert_is_a_number(sig_level)
    assertive.types::assert_is_a_number(p_cutoff)
    assertive.types::assert_is_a_number(q_cutoff)

    plot_type <- match.arg(plot_type, c("barplot", "dotplot1", "dotplot2",
                                        "heatmap", "heatplot"), several.ok = FALSE)
    cluster_type <- match.arg(cluster_type, c("pval", "pval_bins", "bins"), several.ok = FALSE)
    db <- match.arg(db, c("wiki", "kegg", "go",
                          "msigH", "msigC2", "msigC5"), several.ok = FALSE)

    if (is.null(contrastdefs)) {
        contrastdefs <- names(object@metadata$contrastdefs)[1]
    }

    message(paste0("Plotting with following parameters ...",
                   "\ncontrast: ", contrastdefs,
                   "\ncluster type: ", cluster_type,
                   "\ndatabase: ", db))


    if (any(grepl(names(object@metadata), pattern = "ora_")) == FALSE) {
        stop("Object metadata does not contain any over-representation analysis results")
    }
    if (any(grepl(names(object@metadata), pattern = paste0("ora_", db))) == FALSE) {
        stop(paste0("Object metadata does not contain ", db, " over-representation analysis results"))
    }
    cluster_col <- paste0("clust_", cluster_type)
    if (any(grepl(colnames(object@metadata[[paste0("ora_", db)]][[contrastdefs]][[cluster_type]]),
                  pattern = cluster_col)) == FALSE) {
        stop("Can't find cluster definition column in over-representation analysis results")
    }


    n_clust <- object@metadata$clusterdefs[, contrastdefs, paste0("clust_", cluster_type)] %>%
        table()

    labels <- paste0("C", n_clust %>% names(), "\n(n = ", n_clust, ")") %>%
        `names<-`(n_clust %>% names())



    if (plot_type == "barplot") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrastdefs]][[cluster_type]] %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(pvalue) %>%
            dplyr::arrange_at(vars(dplyr::starts_with("clust_"))) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.))))

        p <- tmp %>%
            ggplot(aes(x = Description, y = logp)) +
            geom_bar(stat = "identity") +
            facet_wrap(colnames(tmp)[colnames(tmp) %>% grepl("clust_", .)],
                       nrow = 1,
                       labeller = as_labeller(labels)) +
            coord_flip() +
            labs(x = "Term", y = expression(-log[10]~P~value))
    }



    if (plot_type == "dotplot1") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrastdefs]][[cluster_type]] %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(pvalue) %>%
            dplyr::arrange_at(vars(dplyr::starts_with("clust_"))) %>%
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
            scale_x_discrete(labels = labels) +
            scale_color_continuous(low = "mediumblue",  high = "red2",
                                   space = "Lab", limit = c(0, 10),
                                   breaks = c(1.3, 3, 5, 8),
                                   labels = c("0.05", "0.001", "1e-05", "1e-08")) +
            labs(x = "", y = "Term", color = "p-value", size = "gene ratio")
    }


    if (plot_type == "dotplot2") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrastdefs]][[cluster_type]] %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(gene_ratio = GeneRatio %>%
                              sapply(., function (x) eval(parse(text=x)))) %>%
            dplyr::mutate(logp = ifelse(logp >= 10, 10, logp)) %>%
            dplyr::arrange(dplyr::desc(gene_ratio)) %>%
            dplyr::arrange_at(vars(dplyr::starts_with("clust_"))) %>%
            dplyr::mutate(Description = Description %>%
                              factor(., levels = rev(unique(.))))

        p <- tmp %>%
            ggplot(aes(x = gene_ratio, y = Description, size = Count, color = logp)) +
            facet_wrap(colnames(tmp)[colnames(tmp) %>% grepl("clust_", .)],
                       nrow = 1,
                       labeller = as_labeller(labels)) +
            geom_point() +
            scale_color_continuous(low = "mediumblue",  high = "red2",
                                   space = "Lab", limit = c(0, 10),
                                   breaks = c(1.3, 3, 5, 8),
                                   labels = c("0.05", "0.001", "1e-05", "1e-08")) +
            labs(x = "gene ratio", y = "Term", color = "p-value", size = "count")
    }


        if (plot_type == "heatmap") {

        tmp <- object@metadata[[paste0("ora_", db)]][[contrastdefs]][[cluster_type]] %>%
            dplyr::group_by(Description) %>%
            dplyr::mutate(seen_one = dplyr::cumany(pvalue <= p_cutoff)) %>%
            dplyr::filter(seen_one) %>%
            dplyr::ungroup() %>%
            dplyr::select("Description", "logp", cluster_col) %>%
            dplyr::mutate(logp = ifelse(logp >= 10, 10, logp)) %>%
            reshape2::dcast(as.formula(paste("Description", cluster_col, sep="~")),
                            value.var = "logp") %>%
            replace(is.na(.), 0) %>%
            tibble::column_to_rownames(var = "Description")


        p <- tmp %>% pheatmap::pheatmap(cluster_cols = FALSE,
                                        labels_col = labels,
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


        ora_tmp <- object@metadata[[paste0("ora_", db)]][[contrastdefs]][[cluster_type]] %>%
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

        genes2effect <- object@metadata$limma[, contrastdefs, ] %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "feature_id") %>%
            dplyr::mutate(geneID = gene_id) %>%
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
            theme_bw() +
            theme(#panel.grid.major = element_blank(),
                axis.text.x = element_text(angle = 60, hjust = 1)) +
            labs(x = NULL, y = NULL, legend = "effect")
    }


p


}



