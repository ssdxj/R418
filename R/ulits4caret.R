# train wrapper -----------------------------------------------------------

#' caret::train wrapper to fit biochemphy~index curve
#'
#' @param index the index name for ssdxj_vegindex
#' @param spc  the speclib obj
#' @param biochemphy the target variable name in SI
#'
#' @return a modified caret::obj. With x(index), y(biochemphy), obs, pred in df_inTrain and df_Test
#' @export
do_caret_VI <- function(index, spc, biochemphy, fun_indexCV = get_indexCV) {

  # calc vi
  seed <- 666
  y <- SI(spc)[[biochemphy]] %>% unlist(use.names = FALSE)
  x <- ssdxj_vegindex(index, spc)
  # bind x(index), y(biochemphy) in df(later the df_inTrain and df_Test)
  df <- bind_cols(SI(spc), data.frame(x = x, y = y))
  # bind obs(biochemphy) in df(later the df_inTrain and df_Test)
  df$obs <- y

  # inTrain/Test split
  df_inTrain <- dplyr::filter(df, inTrain == 'inTrain')
  df_Test <- dplyr::filter(df, inTrain == 'Test')
  # inTrain cv split
  indexCV <- fun_indexCV(df_inTrain[[biochemphy]]) # folds = 10, times = 5

  # flag for bear function
  isReverse <- cor(x, y, use = 'pairwise.complete.obs') < 0

  # train control
  tr <- trainControl(index = indexCV, allowParallel = FALSE)


  # place holder
  fit <- NULL

  if(!isReverse){ # positive correlation with bear
    set.seed(seed)
    fit <- tryCatch(
      caret::train(y~x, data = df_inTrain, method = modelInfo_bear(),
                   trControl = tr),
      stop = function(e) {
        msg <- sprintf("Stop in bear fit for index %s!!!", index)
        print(msg)
        print(e)
        return(NULL)
      },
      error = function(e) {
        msg <- sprintf("Error in bear fit for index %s!!!", index)
        print(msg)
        print(e)
        return(NULL)
      }
    )
    if(!is.null(fit)) fit$method <- 'bear'
  } else { # negtive corelation with bearReversed
    set.seed(seed)
    fit <- tryCatch(
      caret::train(y~x, data = df_inTrain, method = modelInfo_bearReverse(),
                   trControl = tr),
      stop = function(e) {
        msg <- sprintf("Stop in bear fit for index %s!!!", index)
        print(msg)
        print(e)
        return(NULL)
      },
      error = function(e) {
        msg <- sprintf("Error in bear fit for index %s!!!", index)
        print(msg)
        print(e)
        return(NULL)
      }
    )
    if(!is.null(fit)) fit$method <- 'bearReverse'
  }

  if(is.null(fit)){ # incase of failure
    fit <- caret::train(y~x ,data = df_inTrain, method = 'lm', trControl = tr)
    fit$method <- 'lm'
  }


  # update param
  fit$xname <- index
  fit$yname <- biochemphy
  df_inTrain$pred <- predict(fit)
  df_Test$pred <- predict(fit, newdata = df_Test)
  fit$df_inTrain <- df_inTrain
  fit$df_Test <- df_Test

  return(fit)
}


#' caret::train wrapper to fit biochemphy~bands machine learning models
#'
#' @param method the caret::train method param
#' @param tuneGrid the caret::trainControl tuneGrid param
#' @param tuneLength the caret::trainControl tuneLength param
#' @param fm the caret::train fm param
#' @param df_inTrain df inTrain
#' @param df_Test df for standalone validation
#' @param indexCV cross-validation split of df_inTrain
#' @param ylabel label for target variable
#'
#' @return
#' @export
do_caret_bands <- function(method, tuneGrid, tuneLength, fm, df_inTrain,
                           df_Test, indexCV, ylabel = NULL,
                           preProcess_fun = c('center', 'scale')){

  # setup cross-validation
  set.seed(324)
  tr <- trainControl(method = 'repeatedcv',
                     index = indexCV,
                     selectionFunction = best,
                     allowParallel = TRUE,
                     search = 'grid'
  )

  # do train
  start_time <- Sys.time()
  if(is.null(tuneGrid)){
    set.seed(666)
    fit <- caret::train(
     fm, data = df_inTrain,
      method = method,
      preProcess = preProcess_fun,
      trControl = tr,
      tuneLength = tuneLength
    )
  } else {
    set.seed(666)
    fit <- caret::train(
      fm, data = df_inTrain,
      method = method,
      preProcess = preProcess_fun,
      trControl = tr,
      tuneGrid = tuneGrid
    )
  }
  end_time <- Sys.time()


  # update param
  # 1: yname, yname label(biochemphyNM)
  yname <- as.character(fm)[2]
  fit$yname <- yname
  ylabel <- ifelse(is.null(ylabel), yname, ylabel)
  fit$ylabel <-  ylabel

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
  obs_max <- max(df_inTrain[[yname]])
  obs_min <- min(df_inTrain[[yname]])
  obs_range <- obs_max - obs_min
  metrics_cv <- fit$bestTune %>% left_join(fit$results) # gof of cv
  df_metrics <- c(RMSE_cv = metrics_cv$RMSE[1],
                  Rsquared_cv = metrics_cv$Rsquared[1],
                  MAE_cv = metrics_cv$MAE[1],
                  n_RMSE_cv = metrics_cv$RMSE[1]/obs_range)

  metrics_train <- postResample(pred = fit$df_inTrain$pred,
                                obs = fit$df_inTrain[[yname]])
  names(metrics_train) <- c('RMSE_train', 'Rsquared_train', 'MAE_train')
  df_metrics <- c(df_metrics, metrics_train, n_RMSE_train = metrics_train[[1]]/obs_range)
  metrics_test <- postResample(pred = fit$df_Test$pred, obs = fit$df_Test[[yname]])
  names(metrics_test) <- c('RMSE_test', 'Rsquared_test','MAE_test')
  df_metrics <- c(df_metrics, metrics_test, n_RMSE_test = metrics_test[[1]]/obs_range)
  fit$df_metrics <- as_tibble_row(df_metrics)

  # 5. plots (p_train, p_test)
  # for axis limits
  dat <- c(fit$df_inTrain$pred, fit$df_Test$pred,
           fit$df_inTrain[[yname]], fit$df_inTrain[[yname]])
  dat_max <- max(dat)
  dat_min <- min(dat)
  dat_range <- dat_max-dat_min
  dat_limits <- c(dat_min - dat_range*0.05, dat_max + dat_range*0.05)
  fit$dat_limits <- dat_limits

  # for axis labels
  obs_label <- sprintf('Observed %s', ylabel)
  pred_label <- sprintf('Predicted %s', ylabel)
  text_cv <- sprintf('nRMSE[cv]==%.2f~R[cv]^2==%.2f', df_metrics[['n_RMSE_cv']],
                     df_metrics[['Rsquared_cv']])
  text_val  <- sprintf('nRMSE[val]==%.2f~R[val]^2==%.2f',
                       df_metrics[['n_RMSE_test']],
                       df_metrics[['Rsquared_test']])

  df_label <- data.frame(inTrain = c('inTrain', 'Test'),
                         label = c(text_cv, text_val))

  if(df_inTrain[['stageF']] %>% unique() %>% length() > 1){
    gp <- geom_point(aes(x = obs, y = pred, color = stageF, shape = stageF))
  } else {
    gp <- geom_point(aes(x = obs, y = pred), shape = 4)
  }

  fit$p <- rbind(df_inTrain, df_Test) %>%
    ggplot() +
    gp +
    geom_abline(intercept = 0, slope = 1, color = 'grey50') +
    geom_text(aes(x = -Inf, y = Inf, label = label),
              hjust = 0, vjust = 1.5, data = df_label,
              parse = TRUE) +
    labs(x = obs_label, y = pred_label) +
    facet_grid(~inTrain) +
    coord_equal(xlim = dat_limits, ylim = dat_limits)

  # final
  return(fit)
}



# inTrain/Test split ------------------------------------------------------



#' inTrain/Test split wrapper -> add inTrain column for a df
#'
#' @param df the target df
#' @param yname the column name will used for stratification
#' @param p the percentage used for standalone validation
#' @param groups n groups for stratification
#'
#' @return df with addictional inTrain column
#' @export
createDataPartition_2col <- function(df, yname, how = 'random',
                                     p = 0.75, groups = 5, step = 4){

  if(how == 'random'){
    df$inTrain <- 'Test'
    inTrain_index <- createDataPartition(df[[yname]], times = 1, p = p,
                                         groups = groups)[[1]]
    df$inTrain[inTrain_index] <- 'inTrain'
  } else if(how == 'order'){
    y <- df[[yname]]
    inTrain <- rep('inTrain', times = length(y))
    y_order <- order(y)
    inTrain[y_order[seq(2, length(y), by = step)]] <- 'Test'
    df$inTrain <- inTrain
  } else {
    stop(sprintf('how = %s is valid in createDataPartition_2col!!!', how))
  }

  return(df)
}


inTrain_split_by_order <- function(y, step){

}


#' inTrain/Test split wrapper -> return list of inTrain/Test PlotID
#'
#' @param df the target df
#' @param yname the column name will used for stratification
#' @param p the percentaget used for standalone validation
#' @param groups n groups for stratification
#'
#' @return list of vector of PlotID fall in inTrain/Test
#' @export
#'
#' @examples
createDataPartition_2PlotIDdict <- function(df, yname, p = 0.66, groups = 5){
  inTrain_index <- createDataPartition(df[[yname]], times = 1, p = p,
                                       groups = groups)[[1]]

  # return
  list(
    inTrain = df$PlotID[inTrain_index],
    Test = df$PlotID[-inTrain_index] # incase NA in original df
  )

}

# fm generaters -----------------------------------------------------------


#' all bands model formula generator
#'
#' @param df the df with columns [yname, B_400, B_401, ...]
#' @param yname the target variable column name
#'
#' @return a formula object
#' @export
generate_fm <- function(df, yname, model = 'B_'){
  if(model == 'B_'){
    xnames <- names(df) %>%
      str_subset('^B_') %>%
      paste(collapse = '+')

    fm <- paste(yname, xnames, sep = '~')
  } else if(model == 'all'){
    df_sub <- dplyr::select(df, -one_of(yname))
    xnames <- names(df_sub) %>% paste(collapse = '+')
    fm <- paste(yname, xnames, sep = '~')
  }

  return(as.formula(fm))
}




# cross-validation index --------------------------------------------------


#' Create repeated n fold corss-validation index
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




# post fit wrappers -------------------------------------------------------


#' collect gof from do_caret_VI and do_caret_bands result
#'
#' @param fit the do_caret_VI and do_caret_bands result
#'
#' @return tibble
#' @export
collect_caret_metrics <- function(fit){

  # for normalized RMSE
  obs_range_inTrain <- max(fit$df_inTrain$obs )-min(fit$df_inTrain$obs)
  obs_range_Test <- max(fit$df_Test$obs)-min(fit$df_Test$obs)

  # CV
  metrics_cv <- fit$bestTune %>% left_join(fit$results) # gof of cv

  df_metrics <- c(RMSE_cv = metrics_cv$RMSE[1],
                  Rsquared_cv = metrics_cv$Rsquared[1],
                  MAE_cv = metrics_cv$MAE[1],
                  n_RMSE_cv = metrics_cv$RMSE[1]/obs_range_inTrain)

  # train
  metrics_train <- postResample(pred = fit$df_inTrain$pred,
                                obs = fit$df_inTrain$obs)
  names(metrics_train) <- c('RMSE_train', 'Rsquared_train', 'MAE_train')
  df_metrics <- c(df_metrics, metrics_train,
                  n_RMSE_train = metrics_train[[1]]/obs_range_inTrain)

  # test
  metrics_test <- postResample(pred = fit$df_Test$pred, obs = fit$df_Test$obs)
  names(metrics_test) <- c('RMSE_test', 'Rsquared_test','MAE_test')
  df_metrics <- c(df_metrics, metrics_test,
                  n_RMSE_test = metrics_test[[1]]/obs_range_Test)

  as_tibble_row(df_metrics)
}


#' collect gof from do_caret_VI and do_caret_bands result(full stages fit)
#'
#' @param fit the do_caret_VI and do_caret_bands result
#'
#' @return tibble
#' @export
collect_caret_metrics_fullStage <- function(fit){
  foo <- function(tag_stage){
    df_inTrain <- dplyr::filter(fit$df_inTrain, stageF == tag_stage)
    df_Test <- dplyr::filter(fit$df_Test, stageF == tag_stage)

    obs_range_inTrain <- max(df_inTrain$obs )-min(df_inTrain$obs)
    obs_range_Test <- max(df_Test$obs)-min(df_Test$obs)

    # train
    df_metrics <- postResample(pred = df_inTrain$pred, obs = df_inTrain$obs)
    names(df_metrics) <- c('RMSE_train', 'Rsquared_train', 'MAE_train')
    df_metrics <- c(df_metrics,  n_RMSE_train = df_metrics[[1]]/obs_range_inTrain)

    # test
    metrics_test <- postResample(pred = df_Test$pred, obs = df_Test$obs)
    names(metrics_test) <- c('RMSE_test', 'Rsquared_test','MAE_test')
    df_metrics <- c(df_metrics, metrics_test,
                    n_RMSE_test = metrics_test[[1]]/obs_range_Test)

    df_metrics[['n_RMSE_test']]
  }

  tags_list <- c("Vegetative",   "Reproductive", "Ripening")
  names(tags_list) <- tags_list
  map_df(tags_list, foo, .id = 'stage')


}




#' generate df_curve from do_caret_VI result
#'
#' @param fit the do_caret_VI fit
#'
#' @return df
#' @export
collect_caret_VIcurve <- function(fit){

  # extract data
  xname <- fit$xname
  yname <- fit$yname
  fit_method <- fit$fit_method
  df_inTrain <- fit$df_inTrain

  # define new x vector
  x_org <- df_inTrain$x
  x_org_max <- max(x_org, na.rm = TRUE)
  x_org_min <- min(x_org, na.rm = TRUE)
  x_sim <- seq(x_org_min, x_org_max, length.out = 100)



  # calc new y vector
  coefs <- coef(fit$finalModel)
  if(fit_method == 'lm'){
    intercept <- coefs[1]
    slope <-coefs[2]
    y_sim <- intercept + slope * x_sim
  } else if(fit_method == 'bear'){
    VI_inf <- coefs[['VI_inf']]
    VI_bare <- coefs[['VI_bare']]
    K_VI <- coefs[['K_VI']]
    y_sim <- -1/(K_VI) * log((VI_inf - x_sim)/(VI_inf - VI_bare))

  } else if(fit_method == 'bearReverse'){
    VI_inf <- coefs[['VI_inf']]
    VI_bare <- coefs[['VI_bare']]
    K_VI <- coefs[['K_VI']]
    y_sim <- -1/(K_VI) * log((x_sim - VI_bare)/(VI_inf - VI_bare))


  } else if(fit_method == 'exp'){
    a <- coefs[['a']]
    b <- coefs[['b']]
    y_sim <- exp(a + b * x_sim)
  } else {
    stop(sprintf('Not valida curve: %s', fit_method))
  }


  # generate df_curve
  data.frame(x = x_sim, y = y_sim)
}


#' obs vs pred plot from do_caret_VI and do_caret_bands result
#'
#' @param fit the do_caret_VI and do_caret_bands result
#' @param RMSE_precision n digits for RMSE annotate
#' @param ylabel target variable label
#' @param isNormalzed control wether the RMSE or nRMSE annotate
#' @param isCV control wether the cv or train gof annotate
#'
#' @return
#' @export
collect_caret_ObsVsPred <- function(fit, RMSE_precision = 2, ylabel = NULL,
                                    isNormalzed = FALSE, isCV = FALSE){
  # for axis limits
  dat_limits <- fit$dat_limits

  # for axis labels
  ylabel <- ifelse(is.null(ylabel), fit$ylabel, ylabel)
  obs_label <- sprintf('Observed %s', ylabel)
  pred_label <- sprintf('Predicted %s', ylabel)



  # for axes annotate
  RMSE_precision <- as.character(RMSE_precision)
  df_metrics <- fit$df_metrics
  if(isNormalzed){
    text_val  <- sprintf(
      str_replace('nRMSE[val]==%.aaaf~R[val]^2==%.2f', 'aaa', RMSE_precision),
      df_metrics[['n_RMSE_test']],
      df_metrics[['Rsquared_test']]
    )

    if(isCV){
      text_cal <- sprintf(
        str_replace('nRMSE[cv]==%.aaaf~R[cv]^2==%.2f', 'aaa', RMSE_precision),
        df_metrics[['n_RMSE_cv']],
        df_metrics[['Rsquared_cv']])

    } else {
      text_cal <- sprintf(
        str_replace('nRMSE[cal]==%.aaaf~R[cal]^2==%.2f', 'aaa', RMSE_precision),
        df_metrics[['n_RMSE_train']],
        df_metrics[['Rsquared_train']])
    }
  } else {
    text_val  <- sprintf(
      str_replace('RMSE[val]==%.aaaf~R[val]^2==%.2f', 'aaa', RMSE_precision),
      df_metrics[['RMSE_test']],
      df_metrics[['Rsquared_test']])

    if(isCV){
      text_cal <- sprintf(
        str_replace('RMSE[cv]==%.aaaf~R[cv]^2==%.2f', 'aaa', RMSE_precision),
        df_metrics[['RMSE_cv']],
        df_metrics[['Rsquared_cv']])

    } else {
      text_cal <- sprintf(
        str_replace('RMSE[cal]==%.aaaf~R[cal]^2==%.2f', 'aaa', RMSE_precision),
        df_metrics[['RMSE_train']],
        df_metrics[['Rsquared_train']])
    }
  }

  df_label <- data.frame(inTrain = c('inTrain', 'Test'),
                         label = c(text_cal, text_val))


  # color and shape
  if(df_inTrain[['stageF']] %>% unique() %>% length() > 1){
    gp <- geom_point(aes(x = obs, y = pred, color = stageF, shape = stageF))
  } else {
    gp <- geom_point(aes(x = obs, y = pred), shape = 4)
  }

  rbind(fit$df_inTrain, fit$df_Test) %>%
    ggplot() +
    gp +
    geom_abline(intercept = 0, slope = 1, color = 'grey50') +
    geom_text(aes(x = -Inf, y = Inf, label = label),
              hjust = 0, vjust = 1.5, data = df_label,
              parse = TRUE) +
    labs(x = obs_label, y = pred_label) +
    facet_grid(~inTrain, labeller = labeller(
      inTrain = c(inTrain = 'Calibration', Test = 'Validation'))) +
    coord_equal(xlim = dat_limits, ylim = dat_limits)


}
