# spc plot ----------------------------------------------------------------

#' plot speclib
#'
#' @param x spc obj
#' @param FUN  ids of spectral(vector) or NULL
#'
#' @return ggplot
#'
spc_plot1 <- function(x, FUN = NULL) {
  # in case only 1 spectra
  if (nspectra(x) == 1) FUN <- 1

  if (is.null(FUN)) {
    # case no FUN
    mean_spec <- apply(x, FUN = mean, na.rm = TRUE)
    sd_spec <- apply(x, FUN = sd, na.rm = TRUE)

    p <- do.call(rbind, list(
      (spectra(mean_spec) + spectra(sd_spec))[1, ],
      spectra(mean_spec)[1, ],
      (spectra(mean_spec) - spectra(sd_spec))[1, ]
    )) %>%
      as_tibble() %>%
      set_colnames(wavelength(x)) %>%
      mutate(fun = c("mean+sd", "mean", "mean-sd")) %>%
      spc_melt() %>%
      ggplot() +
      geom_line(aes(wl, reflect, linetype = fun)) +
      scale_linetype_manual(
        breaks = c("mean+sd", "mean", "mean-sd"),
        labels = c("mean+sd", "mean", "mean-sd"),
        values = c(`mean+sd` = "dashed", `mean-sd` = "dashed", mean = "solid")
      ) +
      labs(x = "Wavelength(nm)", y = "Reflectance") +
      theme_grey() +
      theme(
        legend.title = element_blank(),
        legend.justification = c(1.5, 1.5),
        legend.position = c(1, 1)
      )
  } else {
    if (is.numeric(FUN)) {
      df <- spc_2df(x)[FUN, ]
      IDs <- paste(df$SampleDate, df$PlotID, df$SampleID)
      p <- df %>%
        mutate(IDs = IDs) %>%
        spc_melt() %>%
        ggplot() +
        geom_line(aes(wl, reflect, color = IDs)) +
        labs(x = "Wavelength(nm)", y = "Reflectance") +
        scale_color_manual(values = paletteHm()(length(FUN))) +
        theme_grey() +
        theme(
          legend.title = element_blank(),
          legend.position = "top"
        )
    } else {
      stop("not valied FUN")
    }
  }

  return(p)
}


#' Speclib plot(all or selected PlotID)
#'
#' @param spc the Speclib object
#' @param idCol characters, the colname of groups (ex. PlotID)
#' @param idFilter vector of characters(c("P01", "P02"))
#' @param ylims
#'
#' @return ggplot obj
#' @export
#'
#' @examples
spc_plot2 <- function(spc, idCol = 'PlotID',  idFilter = NULL, ylims = c(0, 0.8)){
  input <- spc_2df_meltWithMask(spc)

  if(!is.null(idFilter)) input <- dplyr::filter(input, .data[[idCol]] %in% idFilter)

  ggplot(input) +
    geom_line(aes_string(x = 'wl', y = 'reflect', group = idCol)) +
    scale_y_continuous(name = 'Reflectance', breaks = pretty_breaks(n = 9),
                       limits = ylims) +
    scale_x_continuous(name = 'Wavelength(nm)', breaks = wl_breaks) +
    theme(
      legend.position = 'top'
    )
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



# caret::train ------------------------------------------------------------

#' self version of plot.train
#'
#' @param obj caret::train obj.
#' \itemize{
#' \item svmRadialSigma
#' \item rf
#' \item keras
#' \item gassprRadial
#' \item pls
#' \item SMLR(stepwise linear regression)
#' }
#'
#' @return p(ggplot2)
#' @export
#'
train_plot <- function(obj) {
  method <- obj$method
  df <- obj$results

  if (method == "svmRadialSigma") {
    p <- df %>%
      mutate(sigma = pretty_values(sigma)) %>%
      mutate(cost = C) %>%
      ggplot(aes(cost, RMSE)) +
      geom_point(aes(color = sigma)) +
      geom_line(aes(color = sigma)) +
      scale_color_npg()
  } else if (method == "nnet") {
    p <- df %>%
      mutate(decay = pretty_values(decay)) %>%
      mutate(units = size) %>%
      ggplot(aes(size, RMSE)) +
      geom_point(aes(color = decay)) +
      geom_line(aes(color = decay)) +
      scale_color_npg()
  } else if (method == "rf") {
    p <- df %>%
      mutate(ntree = as.factor(ntree)) %>%
      ggplot(aes(mtry, RMSE)) +
      geom_point(aes(color = ntree)) +
      geom_line(aes(color = ntree)) +
      scale_color_npg()
  } else if (method == "gaussprRadial") {
    p <- df %>%
      ggplot(aes(sigma, RMSE)) +
      geom_point(color = pal_npg()(1)) +
      geom_line(color = pal_npg()(1))
  } else if (method == "pls") {
    p <- df %>%
      ggplot(aes(ncomp, RMSE)) +
      geom_point(color = pal_npg()(1)) +
      geom_line(color = pal_npg()(1)) +
      scale_x_continuous(breaks = pretty_breaks())
  } else if (method == "SMLR") {
    p <- df %>%
      ggplot(aes(nvmax, RMSE)) +
      geom_point(color = pal_npg()(1)) +
      geom_line(color = pal_npg()(1)) +
      scale_x_continuous(breaks = pretty_breaks())
  } else {
    (
      stop(sprintf("%s: not valid!!!", method))
    )
  }

  return(p)
}


#' wrapper of train_plot for list of train objs
#'
#' @param objs list of train objs
#'
#' @return list of p(ggplot)
#' @export
trains_plot <- function(objs) map(objs, train_plot)


#' scatterplot of obs vs pred from modified version of caret::train
#' (with items: add_trainDf, add_testDf, add_gof)
#'
#' @param train_obj train obj caret::train obj
#' (updated with add_trainDf, add_testDf, add_gof)
#'
#' @return p(ggplot2) or list of p
#'
#' @export
train_obsVSpred_plot <- function(train_obj) {
  trainDf <- train_obj$add_trainDf
  testDf <- train_obj$add_testDf

  # case no test
  if (is.null(testDf)) {
    df <- trainDf
    df_lim <- range(c(df$obs, df$pred))
    gof <- train_obj$add_gof %>%
      mutate(label = sprintf(
        "Rsquared=%.2f\nRMSE=%.2f\nMAE=%.2f",
        TrainRsquared, TrainRMSE, TrainMAE
      ))

    ggplot() +
      geom_point(aes(obs, pred), data = df) +
      geom_text(aes(label = label),
        x = -Inf, y = Inf, size = 5,
        hjust = -0.1, vjust = 1.2,
        data = gof
      ) +
      geom_abline(intercept = 0, slope = 1) +
      coord_equal(xlim = df_lim, ylim = df_lim) +
      ggpubr::theme_pubr()

    # case both train and test
  } else {
    df <- rbind(
      trainDf %>% mutate(key = "Calibration"),
      testDf %>% mutate(key = "Validation")
    )

    df_lim <- df %>%
      group_by(key) %>%
      summarise(amax = max(c(obs, pred)), amin = min(c(obs, pred))) %>%
      gather(tmp, value, amax, amin) %>%
      mutate(obs = value, pred = value)

    gof <- train_obj$add_gof %>%
      mutate(Calibration = sprintf(
        "Rsquared=%.2f\nRMSE=%.2f\nMAE=%.2f",
        TrainRsquared, TrainRMSE, TrainMAE
      )) %>%
      mutate(Validation = sprintf(
        "Rsquared=%.2f\nRMSE=%.2f\nMAE=%.2f",
        TestRsquared, TestRMSE, TestMAE
      )) %>%
      gather(key, label, Calibration, Validation)


    ggplot() +
      geom_point(aes(obs, pred), data = df) +
      geom_text(aes(label = label),
        x = -Inf, y = Inf, size = 5,
        hjust = -0.1, vjust = 1.2,
        data = gof
      ) +
      geom_abline(intercept = 0, slope = 1) +
      geom_blank(aes(obs, pred), data = df_lim) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      facet_grid(~key, scales = "fixed") +
      coord_equal() +
      ggpubr::theme_pubr()
  }
}



#' wrapper of \code{\link{train_obsVSpred_plot}} for list of caret::train objs
#' @param train_objs list of caret::train obj
#' (updated with add_trainDf, add_testDf, add_gof)
#'
#' @return p(ggplot) or list of p
#' @export
#'
trains_obsVSPred_plot <- function(train_objs, ...) {
  out <- list()
  gof_df <- map_df(train_objs, "add_gof", .id = "tag") %>%
    mutate(Clabel = sprintf(
      "Rsquared=%.2f\nRMSE=%.2f",
      CVRsquared, CVRMSE
    ))

  df1 <- map_df(train_objs, "add_trainDf", .id = "tag")
  limDf1 <- df1 %>%
    group_by(tag) %>%
    summarise(
      amax = max(c(obs, pred)),
      amin = min(c(obs, pred))
    ) %>%
    ungroup() %>%
    gather(key, value, amax, amin) %>%
    mutate(obs = value, pred = value)
  p1 <- ggplot(df1, aes(x = obs, y = pred)) +
    geom_point(aes_string(...)) +
    geom_text(aes(label = Clabel),
      size = 3,
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
      parse = FALSE, data = gof_df
    ) +
    geom_abline(intercept = 0, slope = 1, color = "gray50") +
    geom_blank(data = limDf1) +
    facet_wrap(~tag, scales = "fixed", labeller = label_parsed) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    labs(x = "Observed LAI", y = "Predicted LAI") +
    theme(aspect.ratio = 1)

  # incase no testplot
  out <- p1
  # testDf
  df2 <- map_df(train_objs, "add_testDf", .id = "tag")
  if (nrow(df2) > 0) {
    gof_df <- map_df(train_objs, "add_gof", .id = "tag") %>%
      mutate(Vlabel = sprintf(
        "Rsquared=%.2f\nRMSE=%.2f",
        TestRsquared, TestRMSE
      ))


    limDf2 <- group_by(df2, tag) %>%
      summarise(
        amax = max(c(obs, pred)),
        amin = min(c(obs, pred))
      ) %>%
      ungroup() %>%
      gather(key, value, amax, amin) %>%
      mutate(obs = value, pred = value)

    p2 <- ggplot(df2, aes(x = obs, y = pred)) +
      geom_point(aes_string(...)) +
      geom_text(aes(label = Vlabel),
        size = 3,
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
        parse = FALSE, data = gof_df
      ) +
      geom_abline(intercept = 0, slope = 1, color = "gray50") +
      facet_wrap(~tag, scales = "fixed", labeller = label_parsed) +
      geom_blank(data = limDf2) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(breaks = pretty_breaks()) +
      labs(x = "Observed LAI", y = "Predicted LAI") +
      theme(aspect.ratio = 1)

    out <- list(p1 = p1, p2 = p2)
  }

  return(out)
}


#' dotplot of fits gof df
#'
#' @param input multirow gof df
#' @param tagName ID col name(group tag)
#'
#' @return p(ggplot2)
#' @export
#'
trains_gof_dotplot <- function(input, tagName = "method") {
  input %>%
    gather(
      matric, value, contains("RMSE"), contains("MAE"),
      contains("Rsquared")
    ) %>%
    mutate(isTrain = str_extract(matric, "CV|Train|Test")) %>%
    mutate(matric = str_extract(matric, "RMSE|MAE|Rsquared")) %>%
    mutate(line_group = paste(isTrain, matric)) %>%
    ggplot(aes_string(tagName, "value")) +
    geom_point(aes(color = isTrain, shape = isTrain)) +
    geom_line(aes(color = isTrain, linetype = isTrain, group = line_group)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    labs(x = tagName, y = NULL) +
    facet_grid(~matric, scales = "free") +
    coord_flip() +
    themeDotplot +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )
}

#' Residual plot(facet by along ID and along obs value)
#'
#' @param train_obj train obj(updated with add_trainDf, add_testDf, add_gof)
#'
#' @return p(ggplot2)
train_residual_plot <- function(train_obj) {
  trainDf <- train_obj$add_trainDf
  testDf <- train_obj$add_testDf

  if (is.null(testDf)) {
    df1 <- data.frame(x = 1:nrow(trainDf), y = trainDf$pred - trainDf$obs, key2 = "ID")
    df2 <- data.frame(x = trainDf$obs, y = trainDf$pred - trainDf$obs, key2 = "Obs")
    df <- rbind(df1, df2)

    p <- ggplot(df, aes(x, y)) +
      geom_point() +
      geom_hline(yintercept = 0, color = "grey50") +
      facet_grid(~key2, scales = "free_x") +
      labs(x = NULL, y = "Residual")
  } else {
    df1 <- rbind(
      data.frame(
        x = 1:nrow(trainDf), y = trainDf$pred - trainDf$obs,
        key1 = "Calibration", key2 = "ID"
      ),
      data.frame(
        x = 1:nrow(testDf), y = testDf$pred - testDf$obs,
        key1 = "Validation", key2 = "ID"
      )
    )

    df2 <- rbind(
      data.frame(
        x = trainDf$obs, y = trainDf$pred - trainDf$obs,
        key1 = "Calibration", key2 = "Obs"
      ),
      data.frame(
        x = testDf$obs, y = testDf$pred - testDf$obs,
        key1 = "Validation", key2 = "Obs"
      )
    )

    df <- rbind(df1, df2)
    p <- ggplot(df, aes(x, y)) +
      geom_point(aes(color = key1, shape = key1)) +
      geom_hline(yintercept = 0, color = "grey50") +
      facet_grid(. ~ key2, scales = "free_x") +
      labs(x = NULL, y = "Residual")
  }

  return(p)
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
