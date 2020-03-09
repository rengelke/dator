

#' Clean keys
#'
#' @param x vector containing database identifiers
#' @param na.rm logical, remove NA values
#'
#' @return
#' @export
#' @import magrittr
#' @importFrom stringr str_split
#'
#' @examples tidy_keys(c("68195; 100037283" "224903", "NaN", NA, "14660", "", "67532, 100034361"))
tidy_keys <- function (x, na.rm = FALSE) {

    x %<>%
        as.character() %>%
        stringr::str_split(., ";", simplify = TRUE) %>% .[, 1] %>%
        stringr::str_split(., ",", simplify = TRUE) %>% .[, 1]

    x %<>%
        replace(.=="", NA) %>%
        replace(.=="NA", NA) %>%
        replace(.=="NaN", NA)

    if (na.rm) {
        x[!is.na(x)]
    } else {
        x
    }

}

