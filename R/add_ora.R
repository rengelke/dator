




add_ora <- function(object,
                    db = c("wiki", "kegg", "go"),
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
    db <- match.arg(db, c("wiki", "kegg", "go"), several.ok = TRUE)

    if (is.null(contrastdefs)) {
        contrastdefs <- names(object@metadata$contrastdefs)
    }

    organism <- autonomics.annotate::infer_organism(fdata(object)$feature_id[1:3], "uniprot")
    if (organism == "Homo sapiens") {
        organism <- "hsa"
        org_db <- "org.Hs.eg.db"
        if ("wiki" %in% db) wiki_data <- dator:::wiki_human
    }
    if (organism == "Mus musculus") {
        organism <- "mmu"
        org_db <- "org.Mm.eg.db"
        if ("wiki" %in% db) wiki_data <- dator:::wiki_mouse
    }

    wpid2gene <- wiki_data %>% dplyr::select(wpid, gene) #TERM2GENE
    wpid2name <- wiki_data %>% dplyr::select(wpid, name) #TERM2NAME


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
                order(decreasing = TRUE)



            if (db == "wiki") {

                ora_wiki <- no_of_groups %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys()

                    clusterProfiler::enricher(gene = gene_oi,
                                              universe = gene_universe,
                                              pvalueCutoff = 1,
                                              qvalueCutoff = 1,
                                              TERM2GENE = wpid2gene,
                                              TERM2NAME = wpid2name)@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 2),
                                      !!cluster_col := z)

                })
                ora_wiki %<>% do.call(rbind, .)

            }


            if (db == "kegg") {

                ora_kegg <- no_of_groups %>% lapply(function (z) {

                    gene_oi <- tmp_fdata %>% dplyr::filter(clust_id == z) %$%
                        ENTREZID %>%
                        dator::tidy_keys()

                    clusterProfiler::enrichKEGG(gene = gene_oi,
                                                universe = gene_universe,
                                                organism = organism,
                                                pvalueCutoff = 1,
                                                qvalueCutoff = 1,
                                                use_internal_data = FALSE)@result %>%
                        dplyr::mutate(logp = -1*log(pvalue, 2),
                                      !!cluster_col := z)

                })
                ora_kegg %<>% do.call(rbind, .)

            }


            if (db == "go") {

                ora_go <- no_of_groups %>% lapply(function (z) {

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
                ora_go %<>% do.call(rbind, .)

            }




        }





    })


        no_of_grps

        # -----------------------------------------------------------------------------------
        if (dbtype == "wiki")
        {
            wiki_ora_result <- no_of_grps %>% lapply(., function(x)
            {
                gene_oi <-
                    analysis_data %>% dplyr::filter(cluster_id == x & ENTREZID != "")
                # subset_gene_clusts <- subset(analysis_data, select = grep("clust_", names(analysis_data)))
                # gene_oi <- analysis_data %>% dplyr::filter(subset_gene_clusts == x)
                wiki_enrichr_result <- enricher(
                    gene = gene_oi$ENTREZID,
                    universe = gene_universe,
                    pvalueCutoff = 0.05,
                    TERM2GENE = annotate_func(dbtype = "wiki", species = "mouse") %>% dplyr::select(wpid, gene),
                    TERM2NAME = annotate_func(dbtype = "wiki", species = "mouse") %>% dplyr::select(wpid, name)
                )

                wiki_enrichr_result@result %>% dplyr::mutate(p_log = -1 * log(pvalue, 10)) %>%
                    dplyr::filter(p_log > 1.3) -> plot_ora_wiki
                plot_ora_wiki %>% mutate(clust_id = x)
            })
            return(wiki_ora_result)
        }
}


no_of_grps <-
    analysis_data %>% dplyr::select(starts_with("clust_")) %>% unique() %>% .[!is.na(.)]




gene_universe <- object %>% .[, c("ENTREZID")] %>% as.character()


if (db == "kegg")
{
    kk_ora_result <- no_of_grps %>% lapply(., function(x)
    {
        gene_oi <-
            analysis_data %>% dplyr::filter(cluster_id == x & ENTREZID != "")
        # second option - below code creates a vector of differnt column names that hold cluster ids.
        # alternate option - is to create  a vector with different cluster_id_names

        # colnames(analysis_data) -> colnames_storage
        # colnames_storage[grepl("clust_", colnames_storage)] -> col_name_to_loop

        # gene_oi = NULL
        # for (i in cluster_ids) {
        #   gene_oi <- analysis_data %>% dplyr::filter(i == x)
        # }

        kegg_enrichr_result <-
            enrichKEGG(
                gene = gene_oi$ENTREZID ,
                universe = gene_universe,
                organism = "mmu",
                pvalueCutoff = 0.05
            )
        kegg_enrichr_result@result %>% dplyr::mutate(p_log = -1 * log(pvalue, 10)) %>%
            dplyr::filter(p_log > 1.3) -> plot_ora_kegg
        plot_ora_kegg %>% mutate(clust_id = x)

    })
    return(kk_ora_result)
}

# -----------------------------------------------------------------------------------
if (dbtype == "wiki_enrich")
{
    wiki_ora_result <- no_of_grps %>% lapply(., function(x)
    {
        gene_oi <-
            analysis_data %>% dplyr::filter(cluster_id == x & ENTREZID != "")
        # subset_gene_clusts <- subset(analysis_data, select = grep("clust_", names(analysis_data)))
        # gene_oi <- analysis_data %>% dplyr::filter(subset_gene_clusts == x)
        wiki_enrichr_result <- enricher(
            gene = gene_oi$ENTREZID,
            universe = gene_universe,
            pvalueCutoff = 0.05,
            TERM2GENE = annotate_func(dbtype = "wiki", species = "mouse") %>% dplyr::select(wpid, gene),
            TERM2NAME = annotate_func(dbtype = "wiki", species = "mouse") %>% dplyr::select(wpid, name)
        )

        wiki_enrichr_result@result %>% dplyr::mutate(p_log = -1 * log(pvalue, 10)) %>%
            dplyr::filter(p_log > 1.3) -> plot_ora_wiki
        plot_ora_wiki %>% mutate(clust_id = x)
    })
    return(wiki_ora_result)
}
# -----------------------------------------------------------------------------------
if (dbtype == "go_enrich")
{
    go_ora_result <- no_of_grps %>% lapply(., function(x)
    {
        gene_oi <-
            analysis_data %>% dplyr::filter(cluster_id == x & ENTREZID != "")
        print(gene_oi)
        go_enrichr_result <- enrichGO(
            gene = gene_oi$ENTREZID,
            universe = gene_universe,
            OrgDb = org.Mm.eg.db,
            ont = "CC",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.05,
            readable = TRUE
        )
        go_enrichr_result@result %>% dplyr::mutate(p_log = -1 * log(pvalue, 10)) %>%
            dplyr::filter(p_log > 1.0) -> plot_ora_go
        plot_ora_go %>% mutate(clust_id = x)

    })

    return(go_ora_result)
}
}
