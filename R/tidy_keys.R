

#' Clean keys
#'
#' @param x vector containing database identifiers
#'
#' @return
#' @export
#' @import magrittr
#' @importFrom stringr str_split
#'
#' @examples
tidy_keys <- function (x) {

    x %>%
        as.character() %>%
        .[.!=""] %>%
        .[.!="NA"] %>%
        .[.!="NaN"] %>%
        stringr::str_split(., ";", simplify = TRUE) %>% .[, 1] %>%
        stringr::str_split(., ",", simplify = TRUE) %>% .[, 1]


}
