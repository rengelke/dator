prepare_olink <- function (object,
                           filter_sample_type = c("Control", "Sample"),
                           filter_sample_quality = c("Warning", "Pass"),
                           filter_feature_type = c("Control", "Protein"),
                           filter_feature_quality = 1.0,
                           rm_complete_feature_nondetects = TRUE,
                           rm_complete_sample_nondetects = TRUE,
                           rm_single_value_svars = FALSE,
                           rm_single_value_fvars = FALSE,
                           rm_na_svars = TRUE,
                           rm_na_fvars = TRUE,
                           log10transform = FALSE)

{
    SampleType <- RowCheck <- Type <- ColCheck <- NULL
    assertive.types::assert_is_character(c(filter_sample_type,
                                           filter_sample_quality, filter_feature_type))
    assertive.types::assert_is_a_number(filter_feature_quality)
    assertive.types::assert_is_logical(c(rm_complete_feature_nondetects,
                                         rm_complete_sample_nondetects, rm_single_value_svars,
                                         rm_single_value_fvars, rm_na_svars, rm_na_fvars,
                                         log10transform))

    # filter_sample_type
    if ("Assay" %in% autonomics.import::svars(object)) {
        sdata(object)$SampleType <- ifelse(
            sdata(object)$Assay %>%
                stringr::str_detect("^Neg Ctrl |^Ext Ctrl |^IPC "), "Control", "Sample")
        idx <- sdata(object)$SampleType %in% filter_sample_type
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(, idx)
    }
    if ("QC Warning" %in% autonomics.import::svars(object)) {
        idx <- sdata(object)$`QC Warning` %in% filter_sample_quality
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(, idx)
    }
    if ("OlinkID" %in% autonomics.import::fvars(object)) {
        fdata(object)$Type <- ifelse(
            fdata(object)$OlinkID %>% is.na(), "Control", "Protein")
        idx <- fdata(object)$Type %in% filter_feature_type
        idx <- idx & !is.na(idx)
        object %<>% magrittr::extract(idx, )
        fdata(object) %<>% droplevels()
    }
    if ("Missing Data freq." %in% autonomics.import::fvars(object)) {
        idx <- fdata(object)$`Missing Data freq.` <= filter_feature_quality
        object %<>% magrittr::extract(idx, )
        fdata(object) %<>% droplevels()
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
    if (rm_na_fvars)
        autonomics.import::fdata(object) %<>% autonomics.support::rm_na_columns()
    if (rm_na_svars)
        autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns()

    if (log10transform)
        exprs(object) %<>% log10()

    object
}
