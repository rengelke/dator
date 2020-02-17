

#' Title
#'
#' @param color in hexadecimal notation
#' @param c_max integer value indicating maximum chroma value for palette generation
#'
#' @return
#' @export color2palette
#'
#' @import magrittr
#' @import colorspace
#'
#' @examples
#'
color2palette <- function(color, c_max = 140) {

    library(magrittr)

    color <- "#E87721"

    color_hcl <- colorspace::hex2RGB(color) %>%
        as(., "polarLUV") %>% .@coords %>%
        round(0)

    h_seq <- seq(color_hcl[, "H"], color_hcl[, "H"] + 360, by = 15) %>%
        ifelse(. >= 360, . -360, .)


    color_hcl[, "L"]

    find_hcl_color(210, c_max, 62)


}
