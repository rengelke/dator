
#' Title
#'
#' @param h numeric value of hue coordinate
#' @param c_max numeric value of maximum chroma coordinate
#' @param l numeric value of luminescence coordinate
#'
#' @return
#' @export
#'
#' @examples
#'
find_hcl_color <- function (h, c_max, l) {
    if (is.na(colorspace::polarLUV(h, c_max, l) %>% colorspace::hex())) {
        hexi1 <- NA
        cn <- 1
        while (is.na(hexi1)) {
            hexi1 <- colorspace::polarLUV(l, c_max - cn, h) %>% colorspace::hex()
            cn = cn+1
        }
    } else {
        hexi1 <- colorspace::polarLUV(l, c_max, h) %>% hex()
    }
    hexi1
}


