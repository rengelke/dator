
#' Run over-representation analysis and add results
#'
#' @param object SummarizedExperiment
#' @param db annotation database to use: "wiki", "kegg", "go", "msigH", "msigC1", "msigC2"
#' @param cluster_type type of cluster to generate: "pval", "pval_bins", "bins"
#' @param contrastdefs contrast definition(s) for which to perform analysis
#' @param sig_level significance cutoff for p-value based cluster "pval"
#' @param k number of bins for bin cluster "bins"
#' @param ...
#'
#' @return
#' @export
#'
#' @import magrittr
#' @import dplyr
#' @import autonomics
#'
#' @examples add_ora(object, db = "wiki")
add_ora <- function(object,
                    db = c("wiki", "kegg", "go", "msigH", "msigC2", "msigC5"),
                    cluster_type = c("pval", "pval_bins", "bins"),
                    contrastdefs = NULL,
                    sig_level = 0.05,
                    k = 5,
                    ...) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    if (assertthat::has_name(object@metadata, c("contrastdefs", "limma")) == FALSE) {
        stop("Object metadata does not contain limma results")
    }
    assertive.types::assert_is_any_of(db, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(contrastdefs, classes = c("NULL", "character"))
    db <- match.arg(db, c("wiki", "kegg", "go",
                          "msigH", "msigC2", "msigC5"), several.ok = FALSE)

    if (is.null(contrastdefs)) {
        contrastdefs <- names(object@metadata$contrastdefs)
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



        dbres_list <- cluster_type %>% lapply(function (y) {

            cluster_col <- colnames(object@metadata$clusterdefs[, x, ]) %>%
                .[grepl(pattern = "clust_", .)] %>%
                .[y == cluster_type]

            tmp_fdata$clust_id <- object@metadata$clusterdefs[, x, cluster_col]

            no_of_groups <- tmp_fdata$clust_id %>%
                unique() %>% .[!is.na(.)] %>%
                sort()



            if (any(db %in% "wiki")) {

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


            if (any(db %in% "kegg")) {

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


            if (any(db %in% "go")) {

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

            if (any(db %in% "msigH")) {

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

            if (any(db %in% "msigC2")) {

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


            if (any(db %in% "msigC5")) {

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

            })


        names(dbres_list) <- cluster_type
        dbres_list

    })

    names(ora_list) <- contrastdefs

    object@metadata[[length(object@metadata) + 1]] <- ora_list
    names(object@metadata)[length(names(object@metadata))] <- paste0("ora_", db)

    return(object)

}
