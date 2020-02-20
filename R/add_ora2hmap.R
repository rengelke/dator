

library(magrittr)
library(autonomics)
library(ggplot2)


object <- dator:::data_sirt
object %<>% dator::add_ora(db = "wiki")
object %<>% dator::add_ora(db = "kegg")


robust_dist <- function(x, y) {
    qx = quantile(x, c(0.1, 0.9))
    qy = quantile(y, c(0.1, 0.9))
    l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
    x = x[l]
    y = y[l]
    sqrt(sum((x - y)^2))
}


splitFacet <- function(x){
    facet_vars <- names(x$facet$params$facets)         # 1
    x$facet    <- ggplot2::ggplot()$facet              # 2
    datasets   <- split(x$data, x$data[facet_vars])    # 3
    new_plots  <- lapply(datasets,function(new_data) { # 4
        x$data <- new_data
        x})
}


expr(object)
fdata(object) %>% View()


mat <- autonomics::exprs(object)[object@metadata$limma[, "Sirt_vs_Mock", "p"] <= 0.02 &
                              !is.na(object@metadata$limma[, "Sirt_vs_Mock", "p"]), ] %>%
    replace(is.na(.), 0)

mat %>% pheatmap::pheatmap(scale = "row")


# mat %>% ComplexHeatmap::Heatmap(.,
#                                 clustering_distance_rows = robust_dist)





mat %>% ComplexHeatmap::Heatmap(row_km = 4) %>%
    ComplexHeatmap::draw()

ht = ComplexHeatmap::Heatmap(mat, row_km = 4)
ht = ComplexHeatmap::draw(ht)  # draw first to fix k-means clustering
ComplexHeatmap::row_order(ht)

ComplexHeatmap::row_dend(ht)



ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_simple(cbind(1:10, 10:1), pch = 1:2))


add_dclust_ora <- function(object,
                           contrast = NULL,
                           db = c("wiki", "kegg", "go", "msigH", "msigC2", "msigC5"),
                           k = 5,
                           sig_level = 0.05,
                           p_cutoff = 0.05,
                           q_cutoff = 0.05,
                           seed = 2203,
                           ...) {






    mat <- autonomics::exprs(object)[object@metadata$limma[, contrast, "p"] <= sig_level &
                                         !is.na(object@metadata$limma[, contrast, "p"]), ] %>%
        replace(is.na(.), 0)

    set.seed(seed)
    hmap <- mat %>% ComplexHeatmap::Heatmap(row_km = k) %>%
        ComplexHeatmap::draw()

    hmap_order <- ComplexHeatmap::row_order(hmap) %>%
        lapply(data.frame) %>%
        dplyr::bind_rows(.id = "clust_hmap") %>%
        dplyr::arrange_at(ncol(.)) %>%
        dplyr::mutate(feature_id = rownames(mat)) %>%
        dplyr::select(clust_hmap, feature_id)


    clusterdefs <- object@metadata$limma[, contrast, ] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "feature_id")

    clusterdefs %<>%
        dplyr::left_join(., hmap_order, by = "feature_id")


    cluster_col <- "clust_hmap"

    tmp_fdata <- fdata(object)
    gene_universe <- fdata(object)$ENTREZID[clusterdefs$p %>% is.na() %>% `!`] %>%
        dator::tidy_keys()

    tmp_fdata$clust_id <- as.numeric(clusterdefs[, cluster_col])

    no_of_groups <- tmp_fdata$clust_id %>%
        unique() %>% .[!is.na(.)] %>%
        sort()




    if (db == "wiki") {

        ora_result <- no_of_groups %>% lapply(function (z) {

            gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                ENTREZID %>%
                dator::tidy_keys()

            clusterProfiler::enricher(gene = gene_oi,
                                      universe = gene_universe,
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1,
                                      TERM2GENE = db_wiki %>%
                                          dplyr::select(wpid, gene),
                                      TERM2NAME = db_wiki %>%
                                          dplyr::select(wpid, name))@result %>%
                dplyr::mutate(logp = -1*log(pvalue, 2),
                              !!cluster_col := z)

        })
        ora_result %<>% do.call(rbind, .)

    }


    if (db == "kegg") {

        ora_result <- no_of_groups %>% lapply(function (z) {

            gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                ENTREZID %>%
                dator::tidy_keys()

            clusterProfiler::enrichKEGG(gene = gene_oi,
                                        universe = gene_universe,
                                        organism = organism_shrt,
                                        pvalueCutoff = 1,
                                        qvalueCutoff = 1,
                                        use_internal_data = FALSE)@result %>%
                dplyr::mutate(logp = -1*log(pvalue, 2),
                              !!cluster_col := z)

        })
        ora_result %<>% do.call(rbind, .)

    }


    if (db == "go") {

        ora_result <- no_of_groups %>% lapply(function (z) {

            gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                ENTREZID %>%
                dator::tidy_keys()

            clusterProfiler::enrichGO(gene = gene_oi,
                                      universe = gene_universe,
                                      OrgDb = org_db,
                                      ont = "ALL",
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1,
                                      readable = TRUE)@result %>%
                dplyr::mutate(logp = -1*log(pvalue, 2),
                              !!cluster_col := z)

        })
        ora_result %<>% do.call(rbind, .)

    }

    if (db == "msigH") {

        ora_result <- no_of_groups %>% lapply(function (z) {

            gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                ENTREZID %>%
                dator::tidy_keys()

            clusterProfiler::enricher(gene = gene_oi,
                                      universe = gene_universe,
                                      TERM2GENE = db_msigH %>%
                                          dplyr::select(gs_name, entrez_gene),
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)@result %>%
                dplyr::mutate(logp = -1*log(pvalue, 2),
                              !!cluster_col := z)

        })
        ora_result %<>% do.call(rbind, .)

    }

    if (db == "msigC2") {

        ora_result <- no_of_groups %>% lapply(function (z) {

            gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                ENTREZID %>%
                dator::tidy_keys()

            clusterProfiler::enricher(gene = gene_oi,
                                      universe = gene_universe,
                                      TERM2GENE = db_msigC2 %>%
                                          dplyr::select(gs_name, entrez_gene),
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)@result %>%
                dplyr::mutate(logp = -1*log(pvalue, 2),
                              !!cluster_col := z)

        })
        ora_result %<>% do.call(rbind, .)

    }


    if (db == "msigC5") {

        ora_result <- no_of_groups %>% lapply(function (z) {

            gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                ENTREZID %>%
                dator::tidy_keys()

            clusterProfiler::enricher(gene = gene_oi,
                                      universe = gene_universe,
                                      TERM2GENE = db_msigC5 %>%
                                          dplyr::select(gs_name, entrez_gene),
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)@result %>%
                dplyr::mutate(logp = -1*log(pvalue, 2),
                              !!cluster_col := z)

        })
        ora_result %<>% do.call(rbind, .)

    }

    ora_result


    tmp <- ora_result %>%
        dplyr::filter(pvalue <= p_cutoff & qvalue <= q_cutoff) %>%
        #dplyr::filter(!!sym(cluster_col) == gx) %>%
        dplyr::arrange(dplyr::desc(logp)) %>%
        dplyr::arrange(clust_hmap) %>%
        dplyr::mutate(Description = Description %>%
                          factor(., levels = rev(unique(.))))

    #initial groups
    n_groups <- ComplexHeatmap::row_order(hmap) %>% names() %>% as.numeric()
    #groups after ora
    n_groups_now <- tmp$clust_hmap %>% unique()


    plot_list <- n_groups %>%
        lapply(function (xp) {

            if (xp %in% n_groups_now) {
                tmp %>%
                    ggplot(aes(Description, logp)) +
                    geom_bar(stat = "identity") +
                    ggforce::facet_wrap_paginate(clust_hmap ~., ncol = 1, nrow = 1,
                                                 page = grep(xp, n_groups_now))
            } else {
                tmp %>%
                    ggplot(aes(Description, logp)) +
                    ggforce::facet_wrap_paginate(clust_hmap ~.,
                                                 ncol = 1, nrow = 1, page = 1)
            }
            })

    plot_list

    g_facet <- tmp %>%
        ggplot(aes(Description, logp)) +
        geom_bar(stat = "identity") +
        facet_wrap(clust_hmap ~., ncol = 1)
    g_facet



    ComplexHeatmap::Heatmap(mat, row_split = 4)

    g <- g_facet

    Heatmap(mat, right_annotation = HeatmapAnnotation(ggplot = anno_empty(height = unit(6, "cm"))))

    decorate_annotation("ggplot", {
        vp = current.viewport()$name
        print(g, vp = vp)
    })
    cowplot::plot_grid(hmap, g_facet, align = "h")




    gb_hmap <- grid::grid.grabExpr(draw(hmap))
    gb_ggplot <- grid::grid.grabExpr({print(g_facet)})


    grid.newpage()
    pushViewport(viewport(x = 0, y = 1, width = 0.5, height = 1, just = c("left", "top")))
    grid.draw(gb_hmap)
    popViewport()

    pushViewport(viewport(x = 0.5, y = 1, width = 0.5, height = 1, just = c("left", "top")))
    grid.draw(gb_ggplot)
    popViewport()



}













# ### ---------------------------------------------------------------------




add_dclust_ora <- function(object,
                    dclust_object,
                    contrast = NULL,
                    db = c("wiki", "kegg", "go", "msigH", "msigC2", "msigC5"),
                    k = 5,
                    verbose = TRUE,
                    ...) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    assertive.types::assert_is_any_of(dclust_object, classes = c("pheatmap"))
    assertive.types::assert_is_any_of(db, classes = c("NULL", "character"))
    db <- match.arg(db, c("wiki", "kegg", "go",
                          "msigH", "msigC2", "msigC5"), several.ok = FALSE)

    if (is.null(contrast)) {
        contrast <- names(object@metadata$contrastdefs)[1]
        message(paste0("Using '", contrast, "' contrast for analysis"))
    }

    organism <- autonomics.annotate::infer_organism(fdata(object)$feature_id[1:3],
                                                    "uniprot")
    if (organism == "Homo sapiens") {
        organism_shrt <- "hsa"
        org_db <- "org.Hs.eg.db"
        if (db == "wiki") db_wiki <- dator:::wiki_human
        if (db == "msigH") db_msigH <- dator:::msigdb_H_human
        if (db == "msigC2") db_msigC1 <- dator:::msigdb_C2_human
        if (db == "msigC5") db_msigC2 <- dator:::msigdb_C5_human
    }
    if (organism == "Mus musculus") {
        organism_shrt <- "mmu"
        org_db <- "org.Mm.eg.db"
        if (db == "wiki") db_wiki <- dator:::wiki_mouse
        if (db == "msigH") db_msigH <- dator:::msigdb_H_mouse
        if (db == "msigC2") db_msigC2 <- dator:::msigdb_C2_mouse
        if (db == "msigC5") db_msigC5 <- dator:::msigdb_C5_mouse
    }


    if (verbose) {
        cutree(dclust_object$tree_row, k = k) %>% table()
        plot(dclust_object$tree_row)
        abline(h = dclust_object$tree_row$height[length(dclust_object$tree_row$height)-(k-2)],
               col = "red", lty = 2, lwd = 2)
    }

    clusterdefs <- object@metadata$limma[, contrast, ] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "feature_id")

    clusterdefs %<>%
        dplyr::left_join(., cutree(dclust_object$tree_row, k = k) %>%
                             as.data.frame() %>%
                             `colnames<-`(c("clust_dclust")) %>%
                             tibble::rownames_to_column(var = "feature_id"),
                         by = "feature_id")







    object <- add_cluster(object,
                          cluster_type = cluster_type,
                          contrastdefs = contrastdefs,
                          sig_level = sig_level,
                          k = k)

    if (!all(contrastdefs == (attr(object@metadata$clusterdefs, which = "dimnames") %$%
                              contrast))) {
        stop("'contrastedfs' don't match definition in 'clusterdefs'")
        # should not happen since it's always overwritten
        # by add_cluster
    }



    ora_list <- contrastdefs %>% lapply(function (x) {

        tmp_fdata <- fdata(object)
        gene_universe <- fdata(object)$ENTREZID[object@metadata$clusterdefs[, x, "p"] %>%
                                                    is.na() %>% `!`] %>%
            dator::tidy_keys()



        #dbres_list <- cluster_type %>% lapply(function (y) {








            cluster_col <- colnames(object@metadata$clusterdefs[, x, ]) %>%
                .[grepl(pattern = "clust_", .)] %>%
                .[y == cluster_type]

            tmp_fdata$clust_id <- object@metadata$clusterdefs[, x, cluster_col]

            no_of_groups <- tmp_fdata$clust_id %>%
                unique() %>% .[!is.na(.)] %>%
                order(decreasing = TRUE)





            ora_result

            })


        names(dbres_list) <- cluster_type
        dbres_list

    #})

    names(ora_list) <- contrastdefs

    object@metadata[[length(object@metadata) + 1]] <- ora_list
    names(object@metadata)[length(names(object@metadata))] <- paste0("ora_", db)

    return(object)

}
