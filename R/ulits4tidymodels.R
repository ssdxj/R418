# model spec --------------------------------------------------------------


## define model: randomForest
spec_rf <- rand_forest(mtry = tune(),  trees = 500,  min_n = tune()) %>%
  set_engine('ranger') %>%
  set_mode('regression')

### define model: support vector machine
spec_svm <- svm_rbf(cost = tune(),
                    rbf_sigma = tune(),
                    margin = tune()) %>%
  set_engine('kernlab') %>%
  set_mode('regression')

### define model: pls
spec_pls <- pls(num_comp = tune(),
                predictor_prop = tune()) %>%
  set_engine("mixOmics") %>%
  set_mode("regression")



# ML routine --------------------------------------------------------------

#' tidymodel routine for ML fitting
#'
#' @param spec model_spec obj
#' @param res (data) recipe obj
#' @param df_training df for training
#' @param df_testing df for stand-alone validation
#' @param df_folds info for cross-validation
#'
#' @return   list(cv, df_training, df_testing, df_metrics)
#' @export
#' @examples
tidymodels_ML <- function(spec, res, df_training, df_testing, df_folds){

  ## define workflow
  rf_wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(spec)

  ## cross-validaiotn
  registerDoParallel(detectCores())
  rf_cv <- tune_grid(
    rf_wf,
    resamples = df_folds,
    grid = 20,
    control = control_grid(save_pred = TRUE)
  )

  ## finalize
  rf_bestParam <- select_best(rf_cv, metric = 'rmse')
  rf_wf_final <- finalize_workflow(rf_wf, rf_bestParam)
  rf_model <- fit(rf_wf_final, df_training)

  ## df for plot
  rf_df_training <- predict(rf_model, df_training) %>%
    bind_cols(df_training)
  rf_df_testing <- predict(rf_model, df_testing) %>%
    bind_cols(df_testing)

  rf_df_training %>%
    metrics({{biochemphy}}, .pred) %>%
    mutate()


  rf_metrics_cv <- left_join(rf_bestParam, collect_metrics(rf_cv)) %>%
    mutate(.estimate = mean) %>%
    dplyr::select(.metric, .estimator, .estimate) %>%
    mutate(group = 'cv')
  rf_metrics_train <- metrics(rf_df_training, {{biochemphy}}, .pred) %>%
    mutate(group = 'training')
  rf_metrics_test <- metrics(rf_df_testing, {{biochemphy}}, .pred) %>%
    mutate(group = 'testing')

  rf_df_metrics <- bind_rows(rf_metrics_cv, rf_metrics_train, rf_metrics_test) %>%
    dplyr::filter(.metric != 'mae')


  ## output
  list(cv = rf_cv, df_training = rf_df_training, df_testing = rf_df_testing,
       df_metrics = rf_df_metrics)
}


#' tidymodel routine for fitting a linear VI vs biochemphy relationship
#'
#' @param spc the speclib obj
#' @param indexName the index name
#' @param biochemphy the biochemphy name
#'
#' @return a tidymodel obj
#' @export
tidymodels_linearReg <- function(spc, indexName, biochemphy){

  df <- SI(spc)
  df[[indexName]] <- ssdxj_vegindex(indexName, spc)
  df$x <- df[[indexName]]
  df$y <- df[[biochemphy]]

  df_training <- dplyr::filter(df, inTrain == 'inTrain')
  df_testing <- dplyr::filter(df, inTrain == 'Test')

  set.seed(666)
  df_folds <- vfold_cv(df_training, v = 10, repeats = 5, strata = y, breaks = 7)

  spec <- linear_reg() %>%
    set_engine("lm") %>%
    set_mode("regression")

  rec <- recipe(y~x, data = df_training)

  # ## define workflow
  rf_wf <- workflow() %>%
    add_recipe(rec) %>%
    add_model(spec)

  ## cross-validaiotn
  # registerDoParallel(detectCores())
  cl <- makeCluster(detectCores())
  rf_cv <- fit_resamples(
    object = spec,
    preprocessor = rec,
    resamples = df_folds,
    control = control_grid(save_pred = TRUE)
  )
  stopCluster(cl)


  ## finalize
  rf_model <- fit(rf_wf, df_training)

  ## df for plot
  rf_df_training <- predict(rf_model, df_training) %>%
    bind_cols(df_training)
  rf_df_testing <- predict(rf_model, df_testing) %>%
    bind_cols(df_testing)


  rf_metrics_cv <- rf_cv %>%
    dplyr::select(id, id2, .metrics) %>%
    unnest(.metrics) %>%
    group_by(.metric, .estimator) %>%
    summarise(.estimate = mean(.estimate, na.rm = TRUE), .groups = 'drop') %>%
    mutate(group = 'cv')
  rf_metrics_train <- metrics(rf_df_training, {{biochemphy}}, .pred) %>%
    mutate(group = 'training')
  rf_metrics_test <- metrics(rf_df_testing, {{biochemphy}}, .pred) %>%
    mutate(group = 'testing')

  rf_df_metrics <- bind_rows(rf_metrics_cv, rf_metrics_train, rf_metrics_test) %>%
    dplyr::filter(.metric != 'mae')


  rf_model$cv <- rf_cv
  rf_model$df_metrics <- rf_df_metrics
  rf_model$df_training <- rf_df_training
  rf_model$df_testing <- rf_df_testing
  rf_model$xname <- indexName
  rf_model$yname <- biochemphy
  rf_model$func <- 'lm'

  rf_model$p_pre <- spc_index_scatter(spc, indexName, biochemphy)
  rf_model$p_curve <- workflow_curve_plot(rf_model)
  rf_model$p_post <- workflow_post_plot(rf_model)

  return(rf_model)

}





# inTrain/Test split ------------------------------------------------------

.spc_get_notNA <- function(spc, biochemphy){
  spc_2df(spc) %>%
    drop_na(all_of(biochemphy)) %>%
    spc_fromDf()
}


#' add a inTrain column into spc's SI
#'
#' @param spc the speclib obj
#' @param biochemphy the target var
#' @param p for createDataPartion
#' @param groups for createDataPartion
#' @param seed random seed
#'
#' @return new speclib lib
#' @export
#'
#' @examples
spc_add_inTrain <- function(spc, biochemphy, p = 2/3, groups = 7, seed = 324){
  spc <-spc_get_notNA(spc, biochemphy)
  df <- SI(spc)

  set.seed(324)
  inTrain <- createDataPartition(df[[biochemphy]], times = 1, p = p,
                                 list = TRUE, groups = groups)
  inTrain_col <- rep('Test', times = nrow(df))
  inTrain_col[inTrain[[1]]] <- 'inTrain'
  df$inTrain <- inTrain_col

  SI(spc) <- df
  return(spc)
}

#'  add a inTrain column into spc's SI considering the sdate group
#'
#' @param spc the speclib obj
#' @param biochemphy the target var
#' @param p for createDataPartion
#' @param groups for createDataPartion
#' @param seed random seed
#'
#' @return
#' @export
#'
#' @examples
spc_add_inTrain_sdates <- function(spc, biochemphy, p = 2/3, groups = 7,
                                   seed = 324){
  foo <- function(df){
    set.seed(324)
    x <- df[[biochemphy]]
    inTrain <- createDataPartition(x, times = 1,  p = p, list = TRUE,
                                   groups = groups)
    inTrain_col <- rep('Test', times = nrow(df))
    inTrain_col[inTrain[[1]]] <- 'inTrain'

    return(inTrain_col)

  }

  spc <- .spc_get_notNA(spc, biochemphy)
  inTrain <- SI(spc) %>%
    dplyr::select(any_of(c('SampleDate', biochemphy))) %>%
    nest_by(SampleDate) %>%
    mutate(data = list(foo(data))) %>%
    unnest(data)

  df <- SI(spc)
  df$inTrain <- inTrain$data
  SI(spc) <- df
  return(spc)
}

#' a boxplot for check the inTrain/Test split
#'
#' @param spc the speclib
#' @param biochemphy the target var
#'
#' @return
#' @export
#'
#' @examples
spc_test_inTrain <- function(spc, biochemphy){
  SI(spc) %>%
    dplyr::select(any_of(c('SampleDate', biochemphy, 'inTrain'))) %>%
    ggplot() +
    geom_boxplot(aes_string(x = 'inTrain', y = biochemphy)) +
    geom_jitter(aes_string(x = 'inTrain', y = biochemphy),
                width = 0.2, height = 0, alpha = 1/3, color = 'grey50') +
    facet_wrap(~SampleDate, scales = 'free')
}



# post plots --------------------------------------------------------------

#' plot the fitted curve
#'
#' @param fitResult tidymodel obj
#'
#' @return a ggplot obj
#' @export
#'
#' @examples
workflow_curve_plot <- function(fitResult){
  # extract data
  xname <- fitResult$xname
  yname <- fitResult$yname
  func <- fitResult$func
  df_training <- fitResult$df_training

  # define new x vector
  x_org <- df_training[[xname]]
  x_org_max <- max(x_org, na.rm = TRUE)
  x_org_min <- min(x_org, na.rm = TRUE)
  x_sim <- seq(x_org_min, x_org_max, length.out = 100)


  # calc new y vector
  if(func == 'lm'){
    coefs <- fitResult$fit$fit$fit$coefficients
  intercept <- coefs[1]
    slope <-coefs[2]
    y_sim <- intercept + slope * x_sim
  } else {
    stop(func)
  }

  # generate label
rsq <- yardstick::rsq_vec(df_training[[yname]], df_training[['.pred']])
  rmse <- yardstick::rmse_vec(df_training[[yname]], df_training[['.pred']])
  label <- sprintf('RMSE[cal]==%.2f~R[cal]^2==%.2f', rmse, rsq)

  df_sim <- data.frame(x = x_sim, y = y_sim)

  if('stageF' %in% names(df_training)){

    p_curve <- ggplot(df_training, aes(x=x, y=y)) +
      geom_point(aes(color = stageF, shape = stageF)) +
      geom_line(data = df_sim) +
      geom_label(aes(x = -Inf, y = Inf), label = label, hjust = -0.1, vjust = 1,
                 parse = TRUE) +
      labs(x = xname, y = yname) +
      theme(legend.position = 'top',
            legend.title = element_blank())

  } else {
    p_curve <- ggplot(df_training, aes(x=x, y=y)) +
      geom_point() +
      geom_line(data = df_sim) +
      labs(x = xname, y = yname) +
      geom_label(aes(x = -Inf, y = Inf), label = label, hjust = -0.1, vjust = 1,
                 parse = TRUE)


  }

  p_curve
}

#' obs vs pred scatter plot
#'
#' @param fitResult the tidymodel obj
#'
#' @return a ggplot obj
#' @export
#'
#' @examples
workflow_post_plot <- function(fitResult){
  xname <- fitResult$xname
  yname <- fitResult$yname
  df_metrics <- fitResult$df_metrics

  df_trainging <- fitResult$df_training
  df_trainging$group <- 'Calibration'
  df_testing <- fitResult$df_testing
  df_testing$group <- 'Validation'
  df <- rbind(df_trainging, df_testing)

  tmp <- c(df_trainging$y, df_trainging$.pred, df_testing$x, df_testing$y)
  lims <- c(floor(min(tmp)), ceiling(max(tmp)))
  cv_rmse <- dplyr::filter(df_metrics, group == 'cv', .metric == 'rmse')$.estimate[1]
  cv_rsq <- dplyr::filter(df_metrics, group == 'cv', .metric == 'rsq')$.estimate[1]
  cv_label <- sprintf('RMSE[cv]==%.2f~R[cv]^2==%.2f', cv_rmse, cv_rsq)

  val_rmse <- dplyr::filter(df_metrics, group == 'testing', .metric == 'rmse')$.estimate[1]
  val_rsq <- dplyr::filter(df_metrics, group == 'testing', .metric == 'rsq')$.estimate[1]
  val_label <- sprintf('RMSE[val]==%.2f~R[val]^2==%.2f', val_rmse, val_rsq)
  df_label <- data.frame(group = c('Calibration', 'Validation'),
                         label = c(cv_label, val_label))

  ggplot() +
    geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 1) +
    geom_point(aes(x = y, y = .pred, color = stageF, shape = stageF),
               data = df) +
    coord_equal(xlim = lims, ylim = lims) +
    geom_label(aes(x = -Inf, y = Inf, label = label),
               data = df_label, hjust = -0.1, vjust = 1,
                 parse = TRUE) +
    facet_grid(~group) +
    labs(x = sprintf('Observed %s', yname),
         y = sprintf('Predicted %s', yname)) +
    theme(legend.position = 'top',
          legend.title = element_blank())
}


# pre plots ---------------------------------------------------------------


#' scatter plot between biochemphy and index from spc and indexName
#'
#' @param spc  the speclib obj
#' @param indexName the name of index
#' @param biochemphy the name of biochemphy
#'
#' @return a ggplot obj
#' @export
#'
#' @examples
spc_index_scatter <- function(spc, indexName, biochemphy){
  df <- SI(spc)
  df$x <- ssdxj_vegindex(indexName, spc)
  df$y <- df[[biochemphy]]

  p_scatter <- df %>%
    ggplot() +
    geom_point(aes(x = x, y = y,  shape = year, color = stageF))  +
    geom_smooth(aes(x = x, y = y,  color = stageF, linetype = stageF), se = FALSE) +
    geom_smooth(aes(x = x, y = y), color = 'grey50', size = 1, se = FALSE) +
    labs(x = index2label(indexName),
         y = biochemphy2label(biochemphy),
         title = sprintf('Scatter between %s and %s', biochemphy, indexName)) +
    theme(legend.position = 'top')

  p_scatter_facet <- p_scatter + facet_grid(~stageF)

  list(p_scatter = p_scatter, p_scatter_facet = p_scatter_facet)
}




# index test --------------------------------------------------------------

#' A quick check of a new  two band index by give bands and function
#'
#' @param spc the speclib obj
#' @param wl1 band1
#' @param wl2 band2
#' @param fun define of index
#'
#' @return
#' @export
#'
#' @examples
spc_2bandVI_test <- function(spc, wl1, wl2, fun = function(r1, r2) {r1/r2}){

  r1 <- get_reflectance(spc, wl1)
  r2 <- get_reflectance(spc, wl2)

  df <- mutate(SI(spc),
               index = fun(r1, r2), r1 = r1, r2 = r2)


  ggplot() +
    geom_point(aes(x = LCC, y = index, color = stageF, shape = year), data = df) +
    theme(legend.position = 'top')
}
