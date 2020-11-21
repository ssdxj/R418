# spc plot ----------------------------------------------------------------


#' plot selected spectra of a speclib, selected by index vector or a data.frame
#'
#' @param spc speclib obj
#' @param by a (numeric) index vector or data.frame for left_join
#' @param alpha for geom_line
#'
#' @return a ggplot obj
#' @export
#'
#' @examples
spc_plot_filters <- function(spc, by, alpha = 1){
  df <- spc_2df_withMask(spc)

  if(is.data.frame(by)){
    df_filtered <- dplyr::left_join(by, df)
  } else if(is.vector(by) & (is.numeric(by) | is.logical(by))){
    df_filtered <- df[by,]
  } else {
    stop('Wrong filters')
  }

  spc_plot_all(df_filtered, alpha = alpha)

}

#' plot mean, mean+sd and mean-sd spectra of ad speclib
#'
#' @param spc the speclib obj
#'
#' @return a ggplot obj
#' @export
#'
#' @examples
spc_plot_sd <- function(spc){
  df <- spc_2df_withMask(spc)
  df_spectra <- dplyr::select(df, matches('^[[:digit:]]'))
  vec_sd <- map_dbl(df_spectra, sd, na.rm = TRUE)
  vec_mean <- map_dbl(df_spectra, mean, na.rm = TRUE)
  vec_mean[is.nan(vec_mean)] <- NA

  data.frame(wl = parse_double(names(df_spectra)),
             mean = vec_mean,
             upper  =  vec_mean + vec_sd,
             lower  = vec_mean - vec_sd) %>%
    pivot_longer(-wl) %>%
    mutate(factor(name, levels = c('upper', 'mean', 'lower'), labels = c('mean+sd', 'mean', 'mean-sd'))) %>%
    ggplot(aes(x = wl, y = value,  linetype = name, color = name)) +
    geom_line() +
    theme(
      legend.position = 'top',
      legend.title = element_blank()
    )

}


#' plot all spectra of a speclib
#'
#' @param spc the speclib obj
#' @param alpha for geom_line
#'
#' @return a ggplot obj
#' @export
#'
#' @examples
spc_plot_all <- function(spc, alpha = 1/5){
  spc %>%
    spc_2df_withMask() %>%
    spc_melt() %>%
    ggplot() +
    geom_line(aes(x = wl, y = reflect, group = .ID), alpha = alpha) +
    scale_y_continuous(name = 'Reflectance', breaks = scales::pretty_breaks(n = 5)) +
    scale_x_continuous(name = 'Wavelength(nm)', breaks = wl_breaks)

}


# gglm --------------------------------------------------------------------


#' shortcut obs vs pred scatter plot for 'lm' obj
#'
#' @param fit 'lm' obj
#' @param smooth stat_smooth method
#'
#' @return p(ggplot2)
#' @export
#'
#' @examples
#' gglm(lm(y~x))
gglm <- function(fit, smooth = "lm") {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = smooth) +
    labs(title = paste(
      "Adj R2 = ", sprintf("%.3f", summary(fit)$adj.r.squared),
      "Sigma =", sprintf("%.3f", summary(fit)$sigma)
    ))
}



# LbyL --------------------------------------------------------------------

#' heatmap plot of \code{\link{wrapper_LbyL}} function result.
#'
#' @param df \code{\link{wrapper_LbyL}} function result
#'
#' @return p(ggplot2)
#' @export
LbyL_heatmap_R2 <- function(df) {
  if (!inherits(df, "data.frame")) df <- as.data.frame(df)
  ggplot(df) +
    geom_tile(aes(wl2, wl1, fill = r2)) +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(x = "Band2/(nm)", y = "Band1/(nm)", fill = "R2") +
    theme_pubr() +
    theme(
      legend.position = "right"
    )
}
