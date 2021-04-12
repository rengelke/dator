

#' Winsorize based on tail percentile
#'
#' @param x numerical vector
#' @param fraction fraction of the data in the tails to be replaced
#'
#' @return winsorized vector
#' @export
#'
#' @examples winsorize1(rnorm(100, mean = 0, sd = 3), fraction = 0.05)
winsorize1 <- function (x, fraction = 0.05)
{
    if(length(fraction) != 1 || fraction < 0 ||
       fraction > 0.5) {
        stop("bad value for 'fraction'")
    }
    lim <- quantile(x, probs=c(fraction, 1 - fraction), na.rm = TRUE)
    x[ x < lim[1] ] <- lim[1]
    x[ x > lim[2] ] <- lim[2]
    x
}


#' Winsorize based on standard deviation from mean
#'
#' @param x numerical value
#' @param multiple values outside multiple of standard deviation to be replaced
#'
#' @return
#' @export
#'
#' @examples winsorize2(rnorm(100, mean = 0, sd = 3), multiple = 3)
winsorize2 <- function (x, multiple = 3)
{
    if(length(multiple) != 1 || multiple <= 0) {
        stop("bad value for 'multiple'")
    }
    med <- mean(x, na.rm = TRUE)
    y <- x - med
    sc <- mad(y, center=0, na.rm = TRUE) * multiple
    y[ y > sc ] <- sc
    y[ y < -sc ] <- -sc
    y + med
}
