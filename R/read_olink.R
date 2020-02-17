read_olink <- function (file, fid_var = "OlinkID", fname_var = "Assay")
{
    assertive.files::assert_all_are_existing_files(file)
    assertive.types::assert_is_a_string(fid_var)
    assertive.types::assert_is_a_string(fname_var)
    content <- readLines(file)
    n_row <- length(content)
    n_col <- utils::count.fields(file, quote = "", sep = ",") %>%
        max()
    f_row <- 1 + which(stringi::stri_detect_fixed(content, "NPX data"))
    s_row_h <- content %>% magrittr::extract(f_row:length(.)) %>%
        purrr::detect_index(function(y) stringi::stri_detect_regex(y,
                                                                   "QC", negate = FALSE)) %>% magrittr::add(f_row - 1)
    s_row <- content %>% magrittr::extract(f_row:length(.)) %>%
        purrr::detect_index(function(y) stringi::stri_detect_regex(y,
                                                                   "^\\,", negate = FALSE)) %>% magrittr::add(f_row - 1)
    f_col <- 1
    fid_cols <- (f_col + 1):(n_col - 4)
    fid_rows <- content %>% magrittr::extract(f_row:(s_row -
                                                         1)) %>% stringi::stri_extract_first_words() %>% magrittr::equals(fid_var) %>%
        which() %>% magrittr::add(f_row - 1)
    sid_rows <- (1 + s_row):(n_row - 3)
    sid_cols <- 1
    # object <- file %>% autonomics.import:::read_omics(
    #     fid_rows = fid_rows,
    #     fid_cols = fid_cols,
    #     sid_rows = sid_rows,
    #     sid_cols = sid_cols,
    #     expr_rows = (s_row + 1):(n_row - 3),
    #     expr_cols = (f_col + 1):(n_col - 4),
    #     fvar_rows = c(f_row:(s_row - 1), (n_row - 1):n_row),
    #     fvar_cols = f_col,
    #     fdata_rows = c(f_row:(s_row - 1), (n_row - 1):n_row),
    #     fdata_cols = (f_col + 1):(n_col - 4),
    #     svar_rows = s_row_h,
    #     svar_cols = (n_col - 3):n_col,
    #     sdata_rows = (s_row + 1):(n_row - 3),
    #     sdata_cols = (n_col - 3):n_col,
    #     transpose = T, verbose = F)

    object_raw <- read.csv(file, header = FALSE, stringsAsFactors = FALSE,
                           na.strings = c("NA", "NaN", ""))

    expr_rows <- (s_row + 1):(n_row - 3)
    expr_cols <- (f_col + 1):(n_col - 4)

    exprs1 <- object_raw[expr_rows, expr_cols] %>%
        `colnames<-`(object_raw[fid_rows, fid_cols] %>% make.names(unique = TRUE)) %>%
        `rownames<-`(object_raw[sid_rows, sid_cols]) %>%
        data.matrix()

    fvar_rows <- c(f_row:(s_row - 1), (n_row - 1):n_row)
    fvar_cols <- f_col
    fdata_rows <- c(f_row:(s_row - 1), (n_row - 1):n_row)
    fdata_cols <- (f_col + 1):(n_col - 4)


    fdata1 <- object_raw[fdata_rows, fdata_cols] %>% t() %>%
        `rownames<-`(1:nrow(.)) %>%
        `colnames<-`(object_raw[fvar_rows, fvar_cols]) %>%
        as.data.frame(stringsAsFactors = FALSE)
    fdata1$LOD %<>% as.numeric()
    fdata1$`Missing Data freq.` %<>% stringr::str_remove("%") %>%
        as.numeric() %>%
        magrittr::divide_by(100)
    rownames(fdata1) <- fdata1$OlinkID %>% make.names(unique = TRUE)

    svar_rows <- s_row_h
    svar_cols <- c(sid_cols, (n_col - 3):n_col)
    sdata_rows <- (s_row + 1):(n_row - 3)
    sdata_cols <- c(sid_cols, (n_col - 3):n_col)

    sdata1 <- object_raw[sdata_rows, sdata_cols] %>%
        `colnames<-`(object_raw[svar_rows, svar_cols]) %>%
        `rownames<-`(object_raw[sid_rows, sid_cols])

    object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs1),
                                                         rowData = fdata1,
                                                         colData = sdata1)

    # sdata(object) %>% (function(y) {
    #     y$subgroup <- y[[subgroup_var]]
    #     y %>% autonomics.support::pull_columns(c("sample_id",
    #         "subgroup"))
    # })

    fdata(object) %>% (function(y) {
        y$feature_name <- y[[fname_var]]
        y %>% autonomics.support::pull_columns(c("feature_id",
                                                 "feature_name"))
    })
    object
}
