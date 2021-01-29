# global import -----------------------------------------------------------
#' @import tidyverse
#' @import ggsci
#' @import ggpubr
#' @import kableExtra
#' @importFrom hsdar SI spectra wavelength get_reflectance rededge vegindex nri
#' @importFrom caret train preProcess trainControl postResample defaultSummary getTrainPerf
#' @importFrom minpack.lm nlsLM nls.lm.control

library(tidyverse)
library(tidymodels)
library(kableExtra)
library(ggpubr)
library(ggsci)
# library(plsmod)


options(stringsAsFactors = FALSE)


# common var --------------------------------------------------------------

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

#' theme for dotchart plot
#' @return ggplot theme
#' @export
themeDotplot <- function() {
  theme_bw(14) +
    theme(
      axis.text.y = element_text(size = rel(.75)),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = rel(.75)),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.5),
      panel.grid.minor.x = element_blank()
    )
}




`%not in%` <- function(x, table) {
  is.na(match(x, table, nomatch = NA_integer_))
}


df2html <- function(df, digits) {
  df %>%
    mutate_if(is.numeric, round, digits) %>%
    kable(format = "html", digits = digits, escape = FALSE) %>%
    kable_styling(bootstrap_options = c("striped", "hover"))
}


# global vars -------------------------------------------------------------


masks <- c(1350, 1460, 1790, 2000, 2400, 2510)
masks_withoutSWIR <- c(1350, 1460, 1790, 2510)
wl_breaks <- c(350, 470, 560, 680, 800, 1350, 1460, 1790, 2000, 2200, 2500)
# wl_breaks <- c(1100, 1400, 1900, 2200)
sen2_breaks <- c(490, 560, 665, 705, 740, 783, 842, 865, 945, 1610, 2190)
`%not in%` <- Negate(`%in%`)
n_sgolay <- 21

PlotID_ynsd <- seq(1, 12) %>% str_pad(2, 'left', '0') %>% str_pad(3, 'left', 'P')




# tools -------------------------------------------------------------------

#' format digits in training param to str
#'
#' @param x input value
#'
#' @return formated value
pretty_value <- function(x) {
  # for safe
  out <- x

  if (is.numeric(x)) { # numeric
    tmp <- as.character(x)
    if (!str_detect(tmp, "\\.")) { # integer
      out <- tmp
    } else { # float
      if (nchar(str_split(tmp, "\\.")[[1]][2]) <= 2) { # digit number lt 2, ex 0.12
        out <- tmp
      } else { # digit number lg 2, ex 0.123455
        out <- sprintf("%.2e", x)
      }
    }
  }

  return(out)
}

