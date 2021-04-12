#' Add Database Annotation to SummarizedExperiment object
#'
#' @param object SummarizedExperiment object
#' @param get_keys annotation keys to select records for from the database
#'
#' @return modified SummarizedExperiment object
#' @export
#'
#' @examples add_annotation(object, get_keys = c("ENTREZID", "GOALL"))
add_annotation <- function (object,
                            get_keys = c("ENTREZID", "GOALL")) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))

    db <- match.arg(get_keys, c("ENTREZID", "GOALL", "PATH", "PFAM", "PROSITE"),
                    several.ok = TRUE)


    organism <- autonomics.annotate::infer_organism(fdata(object)$feature_uniprot[1:3],
                                                    "uniprot")

    if (organism == "Homo sapiens") {org_db <- org.Hs.eg.db::org.Hs.eg.db}
    if (organism == "Mus musculus") {org_db <- org.Mm.eg.db::org.Mm.eg.db}


    key_uniprot <- fdata(object)$feature_uniprot %>%
            stringr::str_split(., pattern = ";", simplify = TRUE) %>%
        .[, 1]


    annot_list <- list()
    for (i in seq_along(get_keys)) {

        annot_tmp <- AnnotationDbi::select(org_db, keys = key_uniprot,
                                           keytype = "UNIPROT",
                                           columns = c("UNIPROT", get_keys[i]))

        annot_list[[i]] <-annot_tmp %>%
            dplyr::group_by(UNIPROT) %>%
            dplyr::summarize( !!get_keys[i] := !!rlang::sym(get_keys[i]) %>%
                                  paste0(., collapse = "; "))
    }

    annot_list %>% purrr::reduce(dplyr::full_join, by = "UNIPROT") %>%
        dplyr::left_join(fdata(object), .,
                         by = c("feature_uniprot" = "UNIPROT")) %>%
        {. ->> fdata(object)}

    object

}
