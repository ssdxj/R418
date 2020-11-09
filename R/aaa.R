library(tidyverse)
library(tidymodels)
library(kableExtra)
library(ggpubr)
library(ggsci)





# color palette-----------------------------------------------------------------

# colorblind-friendly palette
# from http://jfly.iam.u-tokyo.ac.jp/color/

#' The palette with grey:
#'
#' @return chr list
#' @export
paletteCb <- function() {
  c(
    "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
  )
}

#' The palette with black:
#'
#' @return chr list
#' @export
paletteCbb <- function() {
  c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
  )
}

#' heatmap palette
#'
#' @return function generate by colorRampPalette
#' @export
paletteHm <- function() {
  colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")), space = "Lab")
}

#' palette Jetter
#'
#' @return function generate by colorRampPalette
#' @export
paletteJetter <- function() {
  colorRampPalette(c(
    "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
    "yellow", "#FF7F00", "red", "#7F0000"
  ))
}


# ggplot themes -----------------------------------------------------------



#' theme for dotchart plot
#' @return ggplot theme
#' @export
themeDotplot <- theme_bw(14) +
    theme(
      axis.text.y = element_text(size = rel(.75)),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = rel(.75)),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.5),
      panel.grid.minor.x = element_blank()
    )


# masks -------------------------------------------------------------------

masks <- c(1350, 1460, 1790, 2000, 2400, 2510)
masks_withoutSWIR <- c(1350, 1460, 1790, 2510)



# breaks ------------------------------------------------------------------

wl_breaks <- c(350, 470, 560, 680, 800, 1350, 1460, 1790, 2000, 2200, 2500)
# wl_breaks <- c(1100, 1400, 1900, 2200)
sen2_breaks <- c(490, 560, 665, 705, 740, 783, 842, 865, 945, 1610, 2190)
n_sgolay <- 21

PlotID_ynsd <- seq(1, 12) %>% str_pad(2, 'left', '0') %>% str_pad(3, 'left', 'P')


# universal funs ----------------------------------------------------------

# `%not in%` <- function(x, table) {
#   is.na(match(x, table, nomatch = NA_integer_))
# }

`%not in%` <- Negate(`%in%`)

df2html <- function(df, digits) {
  df %>%
    mutate_if(is.numeric, round, digits) %>%
    kable(format = "html", digits = digits, escape = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover"))
}



# plugins -----------------------------------------------------------------

pls_VIP <- function(object) {
  if (object$method != "oscorespls") {
    stop(
      "Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'"
    )
  }
  if (nrow(object$Yloadings) > 1) {
    stop("Only implemented for single-response models")
  }

  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

## VIPjh returns the VIP of variable j with h components
pls_VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls") {
    stop(
      "Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'"
    )
  }
  if (nrow(object$Yloadings) > 1) {
    stop("Only implemented for single-response models")
  }

  b <- c(object$Yloadings)[1:h]
  T <- object$scores[, 1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[, 1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j, ]^2 / Wnorm2) / sum(SS))
}
