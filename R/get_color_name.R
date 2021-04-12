#' Identify closest match to a named color
#'
#' @param x a color specified as a hex string
#'
#' @return
#' @export
#' @import magrittr
#' @import colorspace
#'
#' @examples get_color_name("#1C52A4")
get_color_name <- function (x) {

    x_hcl <- x %>% colorspace::hex2RGB(.) %>%
        as(., "HSV") %>% .@coords

    x_min <- color_names$HEX %>% colorspace::hex2RGB(.) %>%
        as(., "HSV") %>% .@coords %>%
        apply(., 1, function (y) {
            ((x_hcl[, "H"] - y["H"])^2 +
                 (x_hcl[, "S"] - y["S"])^2 +
                 (x_hcl[, "V"] - y["V"])^2)^0.5 %>%
                as.numeric()

        }) %>% which.min()

    color_names[x_min, "Name"]
}

