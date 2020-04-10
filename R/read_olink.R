#' Read Olink
#'
#' Read data from Olink csv file
#'
#' @param file string: path to Olink NPX data csv file
#'
#' @return SummarizedExperiment
#' @seealso prepare_olink
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @examples read_olink('extdata/olink_mouse_expl.csv')
read_olink <- function (file) {

    # Assert
    assertive.files::assert_all_are_existing_files(file)


    # parameter currently not definable
    # using fixed fid and fname variables
    fid_var <- "OlinkID"
    fname_var <- "Assay"


    # Initial Read
    object_raw <- read.csv(file,
                           header = FALSE,
                           stringsAsFactors = FALSE,
                           blank.lines.skip = FALSE,
                           na.strings = c("NA", "NaN", ""))


    svar_rows	= object_raw[, 1] %>% grep(pattern = "Assay")
    svar_cols	= (object_raw[svar_rows, ] %>% grep(pattern = "Plate ID") %>% .[1]):
        length(object_raw[svar_rows, ])

    fid_rows = object_raw[, 1] %>% grep(pattern = "OlinkID")    #according to fid_var
    fid_cols = 2:(object_raw[svar_rows, ] %>% grep(pattern = "Plate ID") %>% .[1] %>% `-`(1))

    sid_rows = (object_raw[, 1] %>% grep(pattern = "OlinkID") %>% `+`(2)):
        (object_raw[, 1] %>% grep(pattern = "LOD") %>% .[1] %>% `-`(2))
    sid_cols = 1

    expr_rows	= sid_rows
    expr_cols	= fid_cols

    fvar_rows	= c(
        (object_raw[, 1] %>% grep(pattern = "Panel")):
            (object_raw[, 1] %>% grep(pattern = "OlinkID")),
        (object_raw[, 1] %>% grep(pattern = "LOD")):
            (object_raw[, 1] %>% grep(pattern = "Missing Data freq."))
        )
    fvar_cols	= 1

    fdata_rows = fvar_rows
    fdata_cols = fid_cols

    sdata_rows = sid_rows
    sdata_cols = svar_cols

    fname_rows = svar_rows  #according to fname_var
    fname_cols = fid_cols



    # Extract exprs

    fids1 <- object_raw[fid_rows, fid_cols] %>% as.character()

    sids1 <- object_raw[sid_rows, sid_cols] %>% as.character()

    exprs1 <- object_raw[expr_rows, expr_cols] %>%
        data.matrix() %>% t() %>%
        `rownames<-`(fids1 %>% make.names(unique = TRUE)) %>%
        `colnames<-`(sids1)


    # Extract feature annotations
    #    Leave rownames(fdata1) empty: fids1 may contain non-valid values
    #    This happens in MaxQuant files, which sometimes contain missing rows (I think after opening in excell)
    fdata1 <- data.frame(feature_id = fids1, stringsAsFactors = FALSE)
    fdata_available <- !is.null(fvar_rows) & !is.null(fvar_cols)
    if (fdata_available){
        fvars1 <- object_raw[fvar_rows, fvar_cols] %>%
            stringr::str_replace_all(" ", "_") %>%
            stringr::str_replace_all("\\.", "")
        fdata1 %<>% cbind(t(object_raw[fdata_rows, fdata_cols]), stringsAsFactors = FALSE) %>%
            `colnames<-`(c("feature_id", fvars1)) %>%
            data.frame(stringsAsFactors = FALSE, check.names = FALSE)

        if ("LOD" %in% fvars1) {
            fdata1$LOD %<>% as.numeric()
        }
        if ("Missing_Data_freq" %in% fvars1) {
            fdata1$Missing_Data_freq %<>%
                stringr::str_replace("\\%", "") %>%
                as.numeric() %>% `/`(100)
        }
    }




    # Extract sample annotations
    #    Leave rownames(sdata1) empty: sids may contain non-valid values
    #    This happens in SOMA files, where CLIENT_IDENTIFIER is not unique for calibrator and buffer samples
    sdata1 <- data.frame(sample_id = sids1, stringsAsFactors = FALSE)
    sdata_available <- !is.null(svar_rows) & !is.null(svar_cols)
    if (sdata_available){
        svars1 <- object_raw[svar_rows, svar_cols] %>%
            stringr::str_replace_all(" ", "_")
        svars1[grepl(svars1, pattern = "QC_Deviation")] %<>%
            paste0(., c("_Inc_Ctrl", "_Det_Ctrl"))
        sdata1 %<>% cbind(object_raw[sdata_rows, sdata_cols], stringsAsFactors = FALSE) %>%
            `colnames<-`(c("sample_id", svars1)) %>%
            data.frame(stringsAsFactors = FALSE, check.names = FALSE)

        if ("QC_Deviation_from_median_Inc_Ctrl" %in% svars1) {
            sdata1$QC_Deviation_from_median_Inc_Ctrl %<>% as.numeric()
        }
        if ("QC_Deviation_from_median_Inc_Ctrl" %in% svars1) {
            sdata1$QC_Deviation_from_median_Det_Ctrl %<>% as.numeric()
        }
    }


    # Rm samples with missing ids
    idx <- !is.na(sids1)
    if (any(!idx)){
        sids1  %<>% magrittr::extract(idx)
        sdata1 %<>% magrittr::extract(idx,)
        exprs1 %<>% magrittr::extract(, idx)
    }


    # Name features and samples
    fids1 %<>% make.names(unique = TRUE)
    #fids1 %<>% autonomics.support::uniquify('make.unique')  # will not remove first NA
    sids1 %<>% autonomics.support::uniquify('make.unique')
    rownames(exprs1) <- rownames(fdata1) <- fids1
    colnames(exprs1) <- rownames(sdata1) <- sids1



    object <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = exprs1))
    autonomics.import::fdata(object) <- fdata1     # fdata rownames are lost!
    autonomics.import::sdata(object) <- sdata1
    autonomics.import::analysis(object) <- list(
        nfeatures = c(all = nrow(exprs1)),
        nsamples  = c(all = ncol(exprs1)))


    # sdata
    autonomics::sdata(object) %<>% autonomics.support::pull_columns('sample_id')

    # fdata
    assertive.sets::assert_is_subset(fname_var, autonomics::fvars(object))
    assertive.sets::assert_is_subset("Uniprot_ID", autonomics::fvars(object))
    autonomics.import::fdata(object) %<>% (function (y) {
        y$feature_name <- y[[fname_var]]
        y$feature_uniprot <- y[["Uniprot_ID"]]
        y %>% autonomics.support::pull_columns(c('feature_id', 'feature_name',
                                                 'feature_uniprot'))
    })


    # Return
    object

}
