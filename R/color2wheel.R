#' Title
#'
#' @param color color in hexadecimal notation
#' @param num integer value indicating how many colors to be generated for the color wheel
#'
#' @return
#' @export
#'
#' @examples
#'
color2wheel <- function (color, num = 12, plot = FALSE) {


    library(magrittr)

    color <- "#E87721"

    color_hcl <- colorspace::hex2RGB(color) %>%
        as(., "polarLUV") %>% .@coords %>%
        round(0)

    h_seq <- seq(color_hcl[, "H"], color_hcl[, "H"] + 359, by = 1) %>%
        ifelse(.>360, .-360, .)




    find_hcl_color(210, c_max, 62)

    list(
        complementary = c(color, find_hcl_color(h_seq[181], c_max, color_hcl[, "L"]))


    )

    # complementary 1, 181
    # split_comp_30 1, 151, 211
    # split_comp_45 1, 136, 226
    # tetradic      1, 181, 61, 121
    # tetradic_rev  1, 181, 301, 241
    # triadic       1, 121, 241
    # square        1, 91, 181, 271
    # pentagon      1, 73, 145, 217, 289


}
