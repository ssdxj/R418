#' stat_fun for \code{\link{wrapper_LbyL}} to calc r2 of lm fit
#' @export
LbyL_stat_fun_r2 <- function(x, y) {
  summary(lm(y ~ x, data = data.frame(x = x, y = y)))[["r.squared"]]
}

#' resp_fun for \code{\link{wrapper_LbyL}} to calc DVI index
#' @export
LbyL_resp_fun_DVI <- function(b1, b2) {
  b2 - b1
}

#' resp_fun for \code{\link{wrapper_LbyL}} to calc NDVI index
#' @export
LbyL_resp_fun_NDVI <- function(b1, b2) {
  (b2 - b1) / (b2 + b1)
}

#' resp_fun for \code{\link{wrapper_LbyL}} to calc SR index
#' @export
LbyL_resp_fun_SR <- function(b1, b2) {
  b2 / b1
}


# fun ---------------------------------------------------------------------
#' wrapper for doing LbyL VI
#'
#' @param spc spc
#' @param biochemphy colname of responser in SI
#' @param isSym is Symmetry? TRUE/FALSE
#' @param resp_fun function to calc VI, take param \code{b1, b2}, which are
#' reflectance vector corresponding to band1 and band2
#' @param stat_fun function to calc stat result, take param \code{x, y}, which are
#' numeric vector corresponding to predictor(VI) and responder(value of biochemphy)
#'
#' @return df
#' @export
wrapper_LbyL <- function(spc, biochemphy, resp_fun, stat_fun, isSym = TRUE) {

  # prepare data
  rep <- SI(spc)[[biochemphy]]
  wl <- wavelength(spc)
  if (isSym) {
    wlCombns <- combn(wl, 2, simplify = FALSE)
  } else {
    wlCombns <- expand.grid(wl1 = wl, wl2 = wl) %>%
      t() %>%
      as.data.frame() %>%
      map(function(x) unlist(x))
  }

  wrapper_fun <- function(obj) {
    if (obj[[1]] == obj[[2]]) {
      out <- c(obj, NA)
    } else {
      # reflectance of wl1
      b1 <- get_reflectance(spc, obj[1])
      # reflectance of wl2
      b2 <- get_reflectance(spc, obj[2])

      # main part: do calc vi
      vi <- resp_fun(b1 = b1, b2 = b2)

      # fit lm and ready out
      r2 <- stat_fun(vi, rep)

      out <- c(obj, r2)
    }
    return(out)
  }

  plan(multiprocess)
  out <- furrr::future_map(wlCombns, wrapper_fun, .progress = TRUE)

  # tidy out
  df <- do.call(rbind, out)
  colnames(df) <- c("wl1", "wl2", "r2")

  return(df)
}
