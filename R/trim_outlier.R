#' Trim outlier based on standard deviations from the median
#'
#' @param object SummarizedExperiment object with expression values
#' @param n_sigma numerical value indicating how many standard deviations avove and below median define an outlier
#'
#' @return
#' @export
#' @import magrittr
#'
#' @examples
trim_outlier <- function (object, n_sigma = 3.5) {

    exprs(object) %<>% apply(., 1, function (x) {
        ifelse(x > (median(x, na.rm = TRUE) + n_sigma * sd(x, na.rm = TRUE)),
               median(x, na.rm = TRUE) + n_sigma * sd(x, na.rm = TRUE),
               ifelse(x < (median(x, na.rm = TRUE) - n_sigma * sd(x, na.rm = TRUE)),
                      median(x, na.rm = TRUE) - n_sigma * sd(x, na.rm = TRUE),
                      x)
        )}) %>% t()
    object

}
