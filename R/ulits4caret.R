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

#' format digits in training param to str
#'
#' @param x input value vector
#'
#' @return formated value vector
#'
#' @export
pretty_values <- function(x) map_chr(x, pretty_value)




#' caret::train wrapper
#'
#' @param method model(list) or model name
#' @param tuneGrid caret tuneGrid (df)
#' @param tuneLength caret tuneLength (numeric)
#' @param fm formula
#' @param df data
#' @param biochemphyNM y name for label
#'
#' @return modified caret::train obj
#' @export
#'
#' @examples
doCaret <- function(method, tuneGrid, tuneLength, fm, df, biochemphyNM){

  # data split
  df_inTrain <- dplyr::filter(df, inTrain == 'inTrain')
  df_Test <- dplyr::filter(df, inTrain != 'inTrain')

  # setup cross-validation
  set.seed(324)
  ctrl <- trainControl(
    method = 'repeatedcv',
    number = 5,
    repeats = 10,
    p = 0.75,
    verboseIter = FALSE,
    savePredictions = TRUE,
    selectionFunction = 'best',
    allowParallel = TRUE
  )


  # do train
  start_time <- Sys.time()
  fit <- train(
    fm, data = df_inTrain,
    method = method,
    preProcess = c('center', 'scale'),
    trControl = ctrl,
    tuneGrid = tuneGrid,
    tuneLength = tuneLength
  )
  end_time <- Sys.time()


  # add info to caret::train object

  # 1: yname, yname label(biochemphyNM)
  yname <- as.character(attr(fit$terms, 'variables'))[2]
  fit$yname <- yname
  fit$biochemphyNM <- biochemphyNM

  # 2: time_used
  fit$time_used <- end_time - start_time

  # 3. df_inTrain, df_Test
  df_inTrain$pred <- predict(fit, newdata = df_inTrain)
  df_inTrain$obs <- df_inTrain[[yname]]
  df_Test$pred <- predict(fit, newdata = df_Test)
  df_Test$obs <- df_Test[[yname]]
  fit$df_inTrain <- df_inTrain
  fit$df_Test <- df_Test


  # 4. gof(train, validation, cv)
  metrics_cv <- fit$bestTune %>% left_join(fit$results) # gof of cv
  metrics_cv <- c(RMSE_cv = metrics_cv$RMSE[1],
                  Rsquared_cv = metrics_cv$Rsquared[1],
                  MAE_cv = metrics_cv$MAE[1])
  metrics_train <- postResample(pred = fit$df_inTrain$pred, obs = fit$df_inTrain[[yname]])
  names(metrics_train) <- c('RMSE_train', 'Rsquared_train', 'MAE_train')
  metrics_test <- postResample(pred = fit$df_Test$pred, obs = fit$df_Test[[yname]])
  names(metrics_test) <- c('RMSE_test', 'Rsquared_test','MAE_test')
  fit$metrics_df <- rbind(
    as_tibble(metrics_train, rownames = 'metric'),
    as_tibble(metrics_cv, rownames = 'metric'),
    as_tibble(metrics_test, rownames = 'metric')
  )


  # 5. plots (p_train, p_test)
  # for axis limits
  dat <- c(fit$df_inTrain$pred, fit$df_Test$pred, fit$df_inTrain[[yname]], fit$df_inTrain[[yname]])
  dat_max <- max(dat)
  dat_min <- min(dat)
  dat_range <- dat_max-dat_min
  dat_limits <- c(dat_min - dat_range*0.05, dat_max + dat_range*0.05)
  fit$dat_limits <- dat_limits
  # for axis labels
  obs_label <- sprintf('Observed %s', biochemphyNM)
  pred_label <- sprintf('Predicted %s', biochemphyNM)
  text_cv <- sprintf('RMSEcv=%.2f, R2cv=%.2f', metrics_cv[1], metrics_cv[2])
  text_val  <- sprintf('RMSEval=%.2f,R2val=%.2f', metrics_test[1], metrics_test[2])

  fit$p_train <- fit$df_inTrain %>%
    ggplot() +
    geom_point(aes_string(x = yname, y = 'pred'), shape = 4) +
    geom_abline(intercept = 0, slope = 1, color = 'grey50') +
    geom_text(aes(x = -Inf, y = Inf), hjust = 0, vjust = 1.5, label = text_cv) +
    labs(x = obs_label, y = pred_label) +
    coord_equal(xlim = dat_limits, ylim = dat_limits)
  fit$p_test <- fit$df_Test %>%
    ggplot() +
    geom_point(aes_string(x = yname, y = 'pred'), shape = 4) +
    geom_abline(intercept = 0, slope = 1, color = 'grey50') +
    geom_text(aes(x = -Inf, y = Inf), hjust = 0, vjust = 1.5, label = text_val) +
    labs(x = obs_label, y = pred_label) +
    coord_equal(xlim = dat_limits, ylim = dat_limits)

  # final
  return(fit)
}



#' inTrain/Test split wrapper
#'
#' @param df
#' @param yname
#' @param p
#' @param groups
#'
#' @return df with inTrain column
#' @export
#'
#' @examples
doAddinTrain <- function(df, yname, p = 0.66, groups = 5){
  df$inTrain <- 'Test'
  inTrain_index <- createDataPartition(df[[yname]], times = 1, p = p, groups = groups)[[1]]
  df$inTrain[inTrain_index] <- 'inTrain'

  return(df)
}


#' inTrain/Test split wrapper
#'
#' @param df
#' @param yname
#' @param p
#' @param groups
#'
#' @return vector of PlotIDs fall in inTrain/Test
#' @export
#'
#' @examples
doAddinTrain_PlotID <- function(df, yname, p = 0.66, groups = 5){
  inTrain_index <- createDataPartition(df[[yname]], times = 1, p = p,
                                       groups = groups)[[1]]

  # return
  list(
    inTrain = df$PlotID[inTrain_index],
    Test = df$PlotID[-inTrain_index] # incase NA in original df
  )

}


#' Title
#'
#' @param dat
#' @param info
#'
#' @return
#' @export
#'
#' @examples
add_inTrain_col <- function(dat, info){
  dat$inTrain <- ''
  inTrain <- dat$inTrain
  PlotID <- dat$PlotID
  inTrain[PlotID %in% info$inTrain] <- 'inTrain'
  inTrain[PlotID %in% info$Test] <- 'Test'
  dat$inTrain <- inTrain

  return(dat)

}


#' all bands model formula generator
#'
#' @param df
#' @param yname
#'
#' @return
#' @export
#'
#' @examples
generate_fm <- function(df, yname){
  xnames <- names(df) %>%
    str_subset('^B_') %>%
    paste(collapse = '+')

  paste(yname, xnames, sep = '~') %>%
    as.formula()
}


# inTrain index -----------------------------------------------------------

#' get inTrain vector by ordered reponse value
#'
#' @param input numeric vector(normally the response)
#' @param step 1 out of step in Test
#' @param start one of 1:step where the first Test data is.
#'
#' @return index vector for inTrain
#' @export
get_inTrain_respOrder <- function(input, step = 4, start = 2) {
  x <- input
  x_ordered <- order(x)
  inTest <- seq(start, length(x) - 1, by = step)
  out <- x_ordered[-inTest] %>% sort()

  return(out)
}


#' get inTrain vector by ordered reponse value(within each group)
#'
#' @param input numeric vector(normally the response)
#' @param group str vector(ex: stage, year)
#' @param step 1 out of step in Test
#' @param start one of 1:step where the first Test data is.
#'
#' @return index vector for inTrain
#' @export
get_inTrain_respOrder_withGroup <- function(input, group, step = 4, start = 2) {
  idx <- 1:length(input)
  inTrain <- list()

  groupValues <- unique(group)
  for (gv in groupValues) {
    flag <- group == gv
    input_sub <- input[flag]
    idx_sub <- idx[flag]

    input_sub_order <- order(input_sub)

    input_sub_orderd <- input_sub[input_sub_order]
    idx_sub_orderd <- idx_sub[input_sub_order]

    inTrain_sub <- idx_sub_orderd[-seq(start, length(input_sub), by = step)]
    inTrain[[gv]] <- inTrain_sub
  }

  inTrain <- unlist(inTrain, use.names = FALSE) %>% sort()


  return(inTrain)
}



#' Title get inTrain index by full dataset and Test dataset
#'
#' @param input_full full spc (or SI)
#' @param input_Test  test spc (or SI)
#'
#' @return inTrain index vector
#' @export
get_inTrain_fromTestSPC <- function(input_full, input_Test) {
  if (is.speclib(input_full)) input_full <- SI(input_full)
  if (is.speclib(input_Test)) input_Test <- SI(input_Test)
  if (!is.data.frame(input_full) | !is.data.frame(input_Test)) {
    stop("input param error!!!")
  }

  input_full$index <- 1:nrow(input_full)

  out <- anti_join(input_full, input_Test, by = by) %>%
    dplyr::select(index) %>%
    unlist(use.names = FALSE)

  return(out)
}


# cross-validation index --------------------------------------------------


#' Create repeated n fold corssvalidation index
#'
#' @param x numerical vector(normally the inTrain part of response vector)
#' @param fold n folds
#' @param times n times
#'
#' @return list of cv index
#' @export
get_indexCV <- function(x, folds = 10, times = 5) {
  index <- map(3^(1:times - 3), function(seed) {
    set.seed(seed)
    createFolds(x, k = folds, returnTrain = TRUE)
  })
  names(index) <- paste("Round", 1:times, sep = "")
  index <- unlist(index, recursive = FALSE)

  return(index)
}




# post plots --------------------------------------------------------------



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
