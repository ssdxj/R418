# help funs ---------------------------------------------------------------


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



# update train obj --------------------------------------------------------


#' update caret::train obj with add_curveDf, to simulate the trained function (curve)
#' where x is \code{x <- seq(min(df$vi), max(df$vi), length.out = 1000)}, while y is
#' calculated by parse the coeffciences.
#' handled method include: \itemize{
#'   \item method == 'lm'
#'   \item method == 'lm2'
#'   \item method == 'bear'
#' }
#'
#'
#' @param trian_obj: caret::train obj
#'
#' @return df
#' @export
train_add_curveDf <- function(train_obj) {
  md <- train_obj$method
  if (md %not in% c("lm", "lm2", "exp", "bear")) return(NULL)

  # train data
  df <- train_obj$add_trainDf
  # coefs
  coefs <- coef(train_obj$finalModel)

  # simulated obs var
  x <- seq(min(df$vi), max(df$vi), length.out = 1000)
  # simulated pred var
  y <- vector(mode = "double", length = length(x))

  if (md == "lm") {
    # get fitted coefs
    a <- ifelse("(Intercept)" %in% names(coefs), coefs["(Intercept)"], 0)
    b <- coefs["x"]
    # calc y
    y <- a + b * x
  } else if (md == "lm2") {
    a <- ifelse("(Intercept)" %in% names(coefs), coefs["(Intercept)"], 0)
    b <- coefs["x"]
    bb <- coefs["I(x^2)"]
    y <- a + b * x + bb * x^2
  } else if (md == "exp") {
    a <- coefs["a"]
    b <- coefs["b"]
    y <- exp(a + b * x)
    out <- data.frame(x, y)
  } else if (md == "bear") {
    # get fitted coefs
    VI_inf <- coefs["VI_inf"]
    VI_bare <- coefs["VI_bare"]
    K_VI <- coefs["K_VI"]
    # calc y
    y <- -1 / K_VI * log((VI_inf - x) / (VI_inf - VI_bare))
    y <- predict(train_obj, newdata = data.frame(x = x))
    out <- data.frame(x, y)
  } else {
    out <- NULL
  }


  return(out)
}

#' get train and test gof from fit obj (depend add_trainDf and add_testDf in fit obj)
#'
#' @param  train_obj caret::train (caret::rfe) obj
#' (updated version with add_trainDf and add_testDf)
#'
#' @return df
#' @export
train_add_gof <- function(train_obj) {
  perfNames <- train_obj$perfNames

  if (is.null(train_obj$add_trainDf)) return("NA")

  # train gof
  out1 <- defaultSummary(train_obj$add_trainDf)
  names(out1) <- paste("Train", names(out1), sep = "")

  # Test gof
  if (!is.null(train_obj$add_testDf)) {
    out2 <- defaultSummary(train_obj$add_testDf)
    names(out2) <- paste("Test", names(out2), sep = "")
  } else {
    out2 <- NULL
  }

  # Train CV gof
  if(inherits(train_obj, 'train')){
    bestPerf <- train_obj$bestTune
    if (!is.data.frame(bestPerf)) bestPerf <- as.data.frame(bestPerf)
    colnames(bestPerf) <- gsub("^\\.", "", colnames(bestPerf))
    out3 <- merge(train_obj$results, bestPerf) %>% dplyr::select(one_of(perfNames))
    colnames(out3) <- paste("CV", colnames(out3), sep = "")
  } else if(inherits(train_obj, 'rfe')){
    optsize <- train_obj$optsize
    out3 <- dplyr::filter(train_obj$results, Variables == optsize) %>%
      dplyr::select(one_of(perfNames))
    colnames(out3) <- paste("CV", colnames(out3), sep = "")
  }

  if (!is.null(out2)) {
    out <- c(out1, out2, out3) %>% map_df(cbind)
  } else {
    out <- c(out1, out3) %>% map_df(cbind)
  }

  return(out)
}

#' get bear coefs from caret::train obj
#'
#' @param trainObj
#'
#' @return charactor
#' @export
train_add_coefs <- function(trainObj) {
  param <- coef(trainObj$finalModel) %>% round(2)
  param_nm <- names(param)
  param_str <- map2_chr(param_nm, param, ~sprintf("%s=%s", .x, .y))
  add_param <- paste(param_str, collapse = ",")
  return(add_param)
}


#' get (ML) model param
#'
#' @param train_obj caret::train obj
#'
#' @return str
#' @export
train_add_param <- function(trainObj) {
  param <- trainObj$bestTune
  param_nm <- colnames(param)
  param_value <- param[1, ] %>% map(pretty_values)
  param_str <- map2_chr(param_nm, param_value, ~sprintf("%s=%s", .x, .y))
  add_param <- paste(param_str, collapse = ",")
  return(add_param)
}


#' update train obj (add_trainDf, add_testDf, add_gof)
#'
#' @param train_obj caret::train obj
#' @param spc_inTrain Speclib obj where Calibration dataset from
#' @param spc_Test  Speclib obj where Validation dataset from
#' @param newdata in predict(fit, newdata = newdata)
#' @param newdata_inTrain in predict(fit, newdata = newdata)
#' @param biochemphy y name in SI
#'
#' @return updated caret::train Obj
#' @export
train_update <- function(train_obj, spc_inTrain, spc_Test, newdata, biochemphy,
                         newdata_inTrain = NULL) {

  # for special case, need pass in newdata_inTrain by hand
  if (is.null(newdata_inTrain)) {
    predValue <- predict(train_obj)
  } else {
    predValue <- predict(train_obj, newdata_inTrain)
  }
  if (is.matrix(predValue)) predValue <- predValue[, 1]

  # add add_trainDf
  train_obj$add_trainDf <- SI(spc_inTrain) %>%
    mutate(pred = predValue) %>%
    mutate_(obs = biochemphy)


  # add add_testDf
  if (!is.null(spc_Test)) {
    predValue <- predict(train_obj, newdata)
    if (is.matrix(predValue)) predValue <- predValue[, 1]

    train_obj$add_testDf <- SI(spc_Test) %>%
      mutate(pred = predValue) %>%
      mutate_(obs = biochemphy)
  }

  # calc gof
  gof <- train_add_gof(train_obj)
  train_obj$add_gof <- gof


  # calc param
  if(inherits(train_obj, 'train')){
    if (train_obj$method == "bear") {
      param <- train_add_coefs(train_obj)
    } else {
      param <- train_add_param(train_obj)
    }
    train_obj$param <- param
    train_obj$add_gof2 <- c(gof, param = param) %>% as.data.frame()

    # calc curveDf
    curveDf <- train_add_curveDf(train_obj)
    train_obj$add_curveDf <- curveDf
  }



  return(train_obj)
}
