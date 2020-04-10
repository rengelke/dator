
#' Scaling and Centering of Matrix Objects
#'
#' @param x a numeric matrix
#'
#' @return a matrix scaled and centered by rows
#' @export
#'
#' @examples scale_rows(mat)
scale_rows <- function (x) {

    m = apply(x, 1, mean, na.rm = TRUE)
    s = apply(x, 1, sd, na.rm = TRUE)
    return((x - m)/s)
}
