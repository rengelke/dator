#' Add Database Annotation to SummarizedExperiment object
#'
#' @param object SummarizedExperiment object
#'
#' @return modified SummarizedExperiment object
#' @importFrom magrittr %<>% %>%
#' @export
#'
#' @examples add_fasta(object)
add_feature_sequence <- function (object) {


    assertive.types::assert_is_any_of(object, classes = c("SummarizedExperiment"))

    if (!any(colnames(autonomics.import::fdata(object)) %in% "feature_uniprot")) {
        stop("Object requires UniProt ID column (`feature_uniprot`) in feature data")
    }


    organism <- autonomics.annotate::infer_organism(autonomics.import::fdata(object)$feature_uniprot[1:3],
                                                    "uniprot")

    if (organism == "Homo sapiens") {org_db <- dator:::uniprot_human}
    if (organism == "Mus musculus") {org_db <- dator:::uniprot_mouse}

    n <- autonomics.import::fdata(object)$feature_uniprot %>% length() %>%
        `/`(95) %>% ceiling() %>%
        seq_len()

    suppressWarnings(
        uniprot_split <- autonomics.import::fdata(object)$feature_uniprot %>%
            split(f=n)
    )


    prot_seq_list <- uniprot_split %>% lapply(function (x) {
        UniProt.ws::select(org_db,
                           keys = x,
                           keytype = "UNIPROTKB",
                           columns = c("SEQUENCE"))
    })

    prot_seq <- prot_seq_list %>% do.call(rbind, .) %>%
        `colnames<-`(c("UNIPROTKB", "feature_sequence"))


    autonomics.import::fdata(object) %<>%
        dplyr::left_join(., prot_seq[!(prot_seq %>% duplicated()), ],
                         by = c("feature_uniprot" = "UNIPROTKB"))

    object
}
