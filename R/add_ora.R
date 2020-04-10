
#' Run over-representation analysis and add results
#'
#' @param object SummarizedExperiment object
#' @param db annotation database to use: "wiki", "kegg", "go_bp", "go_mf", "go_cc", "msigH", "msigC1", "msigC2"
#' @param cluster_type type of cluster to generate: "pval", "pbins", "bins", "custom"
#' @param contrasts contrasts, as defined in clustedefs of SummarizedExperiment object
#' @param sig_level p-value cutoff for significance cluster "pval"
#' @param p_cutoff p-value cutoff for over-representation analysis results
#' @param q_cutoff q-value cutoff for over-representation analysis results
#' @param simplify_go simplify GO enrichment analysis results
#' @param k number of bins for bin cluster "bins"
#' @param custom_db custom database; data frame containing "id", "gene" and "name" columns
#' @param ... currently not in use
#'
#' @return updated SummarizedExperiment object
#' @export
#'
#' @importFrom magrittr %<>% %>% %$%
#' @import autonomics
#'
#' @examples add_ora(object, db = "wiki")
add_ora <- function(object,
                    db = c("wiki", "kegg", "go_bp", "go_mf", "go_cc", "msigH", "msigC2", "msigC5", "custom"),
                    cluster_type = c("pval", "pbins", "bins", "custom"),
                    contrasts = NULL,
                    sig_level = 0.05,
                    p_cutoff = 0.05,
                    q_cutoff = 0.3,
                    simplify_go = FALSE,
                    breaks = c(0.05, 0.1, 0.2),
                    k = 5,
                    custom_db = NULL,
                    custom_cluster = NULL,
                    ...) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))
    if (assertthat::has_name(object@metadata, c("contrastdefs", "limma")) == FALSE) {
        stop("Object metadata does not contain limma results")
    }
    assertive.types::assert_is_any_of(db, classes = c("NULL", "character"))
    assertive.types::assert_is_any_of(contrasts, classes = c("NULL", "character"))
    db <- match.arg(db, c("wiki", "kegg", "go_bp", "go_mf", "go_cc",
                          "msigH", "msigC2", "msigC5", "custom"), several.ok = FALSE)


    if (db == "custom") {
        if (is.null(custom_db)) {
            stop("Error: `custom` db parameter requires custom database containing term and gene columns.")
        } else {
            custom_db %<>% `colnames<-`(c("term", "gene"))
        }
    }


    if (is.null(custom_cluster)) {
        cluster_type <- cluster_type[cluster_type != "custom"]
    }

    if (is.null(contrasts) & "custom" %in% cluster_type) {
        stop("Please define a single contrast when using `custom` cluster type")
    }
    if ("custom" %in% cluster_type) {
        if (!all(c("feature_id", "clust_custom") %in% colnames(custom_cluster))) {
            stop("`custom_cluster` object must be a data.frame containing `feature_id`
                 and `clust_custom` columns")
        }
    }


    if (is.null(contrasts)) {
        contrasts <- names(object@metadata$contrastdefs)
    }


    organism <- autonomics.annotate::infer_organism(fdata(object)$feature_uniprot[1:3],
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



    object <- dator::add_cluster(object,
                          cluster_type = c("pval", "pbins", "bins", "custom"),
                          contrasts = contrasts,
                          sig_level = sig_level,
                          breaks = breaks,
                          custom_cluster = custom_cluster,
                          k = k)

    if (!all(contrasts == (attr(object@metadata$clusterdefs, which = "names")))) {
        stop(paste0("'contrasts' ", paste0("(", paste(contrasts, collapse = ", "), ")" ),
                    " do not match 'contrasts' in cluster definitions"))
    }



    ora_list <- contrasts %>% lapply(function (x) {

        tmp_fdata <- fdata(object)
        tmp_fdata$ENTREZID_cl <- tmp_fdata$ENTREZID %>%
            dator::tidy_keys(na.rm = FALSE)
        gene_universe <- fdata(object)$ENTREZID[object@metadata$clusterdefs[[x]][, "p"] %>%
                                                    is.na() %>% `!`] %>%
            dator::tidy_keys(na.rm = TRUE)



        dbres_list <- cluster_type %>% lapply(function (y) {

            cluster_col <- colnames(object@metadata$clusterdefs[[x]]) %>%
                .[grepl(pattern = paste0("clust_", y), .)]

            tmp_fdata$clust_id <- object@metadata$clusterdefs[[x]][, cluster_col[1]]

            group_ids <- tmp_fdata$clust_id %>%
                unique() %>% .[!is.na(.)] %>%
                sort()



            if (db == "wiki") {

                db_wiki %>%
                    dplyr::filter(gene %in% tmp_fdata$ENTREZID_cl) %>%
                    dplyr::group_by(gene) %>%
                    dplyr::summarize(ENTREZID_cl = unique(gene),
                              ID = paste(wpid, collapse = "/")) %>%
                    dplyr::select(ENTREZID_cl, ID) %>%
                    dplyr::left_join(tmp_fdata, ., by = "ENTREZID_cl")

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    clusterProfiler::enricher(gene = gene_oi,
                                              universe = gene_universe,
                                              pvalueCutoff = p_cutoff,
                                              qvalueCutoff = q_cutoff,
                                              TERM2GENE = db_wiki %>%
                                                  dplyr::select(wpid, gene),
                                              TERM2NAME = db_wiki %>%
                                                  dplyr::select(wpid, name))@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }


            if (db == "kegg") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    clusterProfiler::enrichKEGG(gene = gene_oi,
                                                universe = gene_universe,
                                                organism = organism_shrt,
                                                pvalueCutoff = p_cutoff,
                                                qvalueCutoff = q_cutoff,
                                                use_internal_data = FALSE)@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }


            if (db == "go_bp") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    tmp_res <- clusterProfiler::enrichGO(gene = gene_oi,
                                                         universe = gene_universe,
                                                         OrgDb = org_db,
                                                         ont = "BP",
                                                         pAdjustMethod = "none",
                                                         pvalueCutoff = p_cutoff,
                                                         qvalueCutoff = q_cutoff,
                                                         readable = TRUE)
                    if (simplify_go) {tmp_res %<>% clusterProfiler::simplify(by = "pvalue")}

                    tmp_res@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }

            if (db == "go_mf") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    tmp_res <- clusterProfiler::enrichGO(gene = gene_oi,
                                                         universe = gene_universe,
                                                         OrgDb = org_db,
                                                         ont = "MF",
                                                         pAdjustMethod = "none",
                                                         pvalueCutoff = p_cutoff,
                                                         qvalueCutoff = q_cutoff,
                                                         readable = TRUE)
                    if (simplify_go) {tmp_res %<>% clusterProfiler::simplify(by = "pvalue")}

                    tmp_res@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }

            if (db == "go_cc") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    tmp_res <- clusterProfiler::enrichGO(gene = gene_oi,
                                                         universe = gene_universe,
                                                         OrgDb = org_db,
                                                         ont = "CC",
                                                         pAdjustMethod = "none",
                                                         pvalueCutoff = p_cutoff,
                                                         qvalueCutoff = q_cutoff,
                                                         readable = TRUE)
                    if (simplify_go) {tmp_res %<>% clusterProfiler::simplify(by = "pvalue")}

                    tmp_res@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }

            if (db == "msigH") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    clusterProfiler::enricher(gene = gene_oi,
                                              universe = gene_universe,
                                              TERM2GENE = db_msigH %>%
                                                  dplyr::select(gs_name, entrez_gene),
                                              pvalueCutoff = p_cutoff,
                                              qvalueCutoff = q_cutoff)@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }

            if (db == "msigC2") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    clusterProfiler::enricher(gene = gene_oi,
                                              universe = gene_universe,
                                              TERM2GENE = db_msigC2 %>%
                                                  dplyr::select(gs_name, entrez_gene),
                                              pvalueCutoff = p_cutoff,
                                              qvalueCutoff = q_cutoff)@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }


            if (db == "msigC5") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    clusterProfiler::enricher(gene = gene_oi,
                                              universe = gene_universe,
                                              TERM2GENE = db_msigC5 %>%
                                                  dplyr::select(gs_name, entrez_gene),
                                              pvalueCutoff = p_cutoff,
                                              qvalueCutoff = q_cutoff)@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }

            if (db == "custom") {

                ora_result <- group_ids %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys(na.rm = TRUE)

                    clusterProfiler::enricher(gene = gene_oi,
                                              universe = gene_universe,
                                              TERM2GENE = custom_db %>%
                                                  dplyr::select(term, gene),
                                              # TERM2NAME = custom_db %>%
                                              #     dplyr::select(term, name),
                                              pvalueCutoff = p_cutoff,
                                              qvalueCutoff = q_cutoff)@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 10),
                                      !!cluster_col[1] := z)

                })
                ora_result %<>% do.call(rbind, .)

            }

            ora_result

            })


        names(dbres_list) <- cluster_type
        dbres_list

    })

    names(ora_list) <- contrasts

    ora_name <- paste0("ora_", db)

    if (ora_name %in% names(object@metadata)) {
        object@metadata[[ora_name]] <- ora_list
    } else {
        object@metadata[[length(object@metadata) + 1]] <- ora_list
        names(object@metadata)[length(names(object@metadata))] <- ora_name
    }

    return(object)
}

