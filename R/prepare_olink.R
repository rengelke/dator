

# file <- "C:/Users/rue2006/Documents/R_wd/mazloum_mouse/data/20200203_WCQO-20-001B/20200203_WCQO-20-001B.csv"
# object <- read_olink(file)
#
# filter_sample_type = c("Sample")
# filter_feature_type = c("Protein")
# filter_sample_quality = c("Warning", "Pass")
# filter_feature_quality = 1.0
# rm_features_below_lod_for_some_sample = TRUE
#
# rm_complete_feature_nondetects = TRUE
# rm_complete_sample_nondetects = TRUE
# rm_single_value_svars = FALSE
# rm_single_value_fvars = FALSE
# rm_na_svars = TRUE
# rm_na_fvars = TRUE
# log10transform = FALSE


#' Prepare Olink
#'
#' @param object SummarizedExperiment
#' @param filter_sample_type string vector: sample  types to be filtered for. Subset of c('Sample', 'Control').
#' @param filter_feature_type string vector: feature types to be filtered for. Subset of c('Protein', 'Control')
#' @param filter_sample_quality string vector: sample  qualities to be filtered for. Subset of c('Pass', 'Warning')
#' @param filter_feature_quality number: fraction 0-1 indicating missing data frequency cutoff
#' @param rm_features_below_lod_for_some_sample logical: wether to remove values below limit of detection (LOD)
#' @param rm_complete_feature_nondetects logical: wether to remove complete feature nondetects
#' @param rm_complete_sample_nondetects logical: wether to remove complete sample nondetects
#' @param rm_single_value_svars logical:  whether to remove single value svars
#' @param rm_single_value_fvars logical:  whether to remove single value fvars
#'
#' @return
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @examples prepare_olink(object)
prepare_olink <- function (object,
                           filter_sample_type = c("Sample", "Control"),
                           filter_feature_type = c("Protein", "Control"),
                           filter_sample_quality = c("Pass", "Warning"),
                           filter_feature_quality = 1.0,
                           rm_features_below_lod_for_some_sample = TRUE,
                           rm_complete_feature_nondetects = TRUE,
                           rm_complete_sample_nondetects = TRUE,
                           rm_single_value_svars = FALSE,
                           rm_single_value_fvars = FALSE) {


    c(filter_sample_type, filter_sample_quality, filter_feature_type) %>%
        assertive.types::assert_is_character()
    filter_feature_quality %>%
        assertive.types::assert_is_a_number() %>%
        assertive.numbers::assert_all_are_greater_than(0) %>%
        assertive.numbers::assert_all_are_less_than_or_equal_to(1)
    assertive.types::assert_is_logical(c(rm_features_below_lod_for_some_sample,
                                         rm_complete_feature_nondetects,
                                         rm_complete_sample_nondetects,
                                         rm_single_value_svars,
                                         rm_single_value_fvars))

    filter_sample_type %<>% match.arg(c("Sample", "Control"), several.ok = TRUE)
    filter_feature_type %<>% match.arg(c("Protein", "Control"), several.ok = TRUE)
    filter_sample_quality %<>% match.arg(c("Pass", "Warning"), several.ok = TRUE)



    # filter_sample_type
    if ("sample_id" %in% autonomics.import::svars(object)) {
        autonomics.import::sdata(object)$sample_type <- ifelse(
            autonomics.import::sdata(object)$sample_id %>%
                stringr::str_detect("^Neg Ctrl |^Ext Ctrl |^IPC "),
            "Control", "Sample")
        idx <- sdata(object)$sample_type %in% filter_sample_type
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(, idx)
    }

    #filter_feature_type
    if ("Assay" %in% autonomics.import::fvars(object)) {
        autonomics.import::fdata(object)$feature_type <- ifelse(
            autonomics.import::fdata(object)$Assay %>%
                stringr::str_detect("^Inc Ctrl|^Det Ctrl|^Ext Ctrl"),
            "Control", "Protein")
        idx <- fdata(object)$feature_type %in% filter_feature_type
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(idx, )
        autonomics.import::fdata(object) %<>% droplevels()
    }

    if ("QC_Warning" %in% autonomics.import::svars(object)) {
        idx <- autonomics.import::sdata(object)$QC_Warning %in% filter_sample_quality
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(, idx)
    }

    if ("Missing_Data_freq" %in% autonomics.import::fvars(object)) {
        idx <- autonomics.import::fdata(object)$Missing_Data_freq <= filter_feature_quality
        object %<>% magrittr::extract(idx, )
        autonomics.import::fdata(object) %<>% droplevels()
    }

    if (rm_features_below_lod_for_some_sample) {
        exprs(object) %>% sweep(., 2, fdata(object)$LOD, `-`, check.margin = FALSE) %>%
            replace(. <= 0, NA)
    }

    if (rm_complete_feature_nondetects) {
        idx <- exprs(object) %>% apply(., 1, function (x) all(x  == 0 | is.na(x)))
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(!idx, )
    }
    if (rm_complete_sample_nondetects) {
        idx <- exprs(object) %>% apply(., 2, function (x) all(x  == 0 | is.na(x)))
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(, !idx)
    }

    if (rm_single_value_svars)
        autonomics.import::sdata(object) %<>% autonomics.support::rm_single_value_columns()
    if (rm_single_value_fvars)
        autonomics.import::fdata(object) %<>% autonomics.support::rm_single_value_columns()

    object
}
