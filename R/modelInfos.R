#' modelInfo for caret doing bear's law VI curve fitting
#'
#' @return the modelInfo(list)
#' @export
modelInfo_bear <- function() {
  modelInfo_bear <- list()
  modelInfo_bear$label <- "VI bear law fit by nlsLM"
  modelInfo_bear$library <- "minpack.lm"
  modelInfo_bear$loop <- NULL
  modelInfo_bear$type <- "Regression"
  modelInfo_bear$prob <- NULL
  modelInfo_bear$parameters <- data.frame(
    parameter = c("parameter"),
    class = c("character"),
    label = c("Parameter")
  )

  # param grid
  modelInfo_bear$grid <- function(x, y, len = NULL, search = "grid") data.frame(parameter = "none")

  # fit
  modelInfo_bear$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    dat <- data.frame(x = x)
    dat$.outcome <- y

    x_upper <- max(x)
    x_lower <- min(x)
    x_sd <- sd(x)

    nlsLM(.outcome ~ -1 / (K_VI) * log((VI_inf - x) / (VI_inf - VI_bare)),
      data = dat,
      trace = FALSE,
      control = nls.lm.control(maxiter = 1024),
      start = c(VI_inf = x_upper + x_sd * 3, VI_bare = x_lower - x_sd * 3, K_VI = 0.2),
      lower = c(VI_inf = x_upper, VI_bare = x_lower - x_sd * 10, K_VI = -1),
      upper = c(VI_inf = x_upper + x_sd * 10, VI_bare = x_lower, K_VI = 3)
    )
  }

  # predict
  modelInfo_bear$predict <- function(modelFit, newdata, submodels = NULL) {
    if (!is.data.frame(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    predict(modelFit, newdata)
  }

  modelInfo_bear$tags <- c("Linear Regression", "Robust Model", "Accepts Case Weights")
  modelInfo_bear$sort <- function(x) x
  modelInfo_bear$prob <- function(x) x

  return(modelInfo_bear)
}



#' modelInfo for caret doing bear's law VI curve fitting
#'
#' @return the modelInfo(list)
#' @export
modelInfo_bearReverse <- function() {
  modelInfo_bear <- list()
  modelInfo_bear$label <- "VI bear law fit by nlsLM"
  modelInfo_bear$library <- "minpack.lm"
  modelInfo_bear$loop <- NULL
  modelInfo_bear$type <- "Regression"
  modelInfo_bear$prob <- NULL
  modelInfo_bear$parameters <- data.frame(
    parameter = c("parameter"),
    class = c("character"),
    label = c("Parameter")
  )

  # param grid
  modelInfo_bear$grid <- function(x, y, len = NULL, search = "grid") data.frame(parameter = "none")

  # fit
  modelInfo_bear$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    dat <- data.frame(x = x)
    dat$.outcome <- y

    x_upper <- max(x)
    x_lower <- min(x)
    x_sd <- sd(x)

    nlsLM(.outcome ~ -1 / (K_VI) * log((x - VI_bare) / (VI_inf - VI_bare)),
      data = dat,
      trace = FALSE,
      control = nls.lm.control(maxiter = 1024),
      start = c(VI_inf = x_upper + x_sd * 3, VI_bare = x_lower - x_sd * 3, K_VI = 0.2),
      lower = c(VI_inf = x_upper, VI_bare = x_lower - x_sd * 10, K_VI = -1),
      upper = c(VI_inf = x_upper + x_sd * 10, VI_bare = x_lower, K_VI = 3)
    )
  }

  # predict
  modelInfo_bear$predict <- function(modelFit, newdata, submodels = NULL) {
    if (!is.data.frame(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    predict(modelFit, newdata)
  }

  modelInfo_bear$tags <- c("Linear Regression", "Robust Model", "Accepts Case Weights")
  modelInfo_bear$sort <- function(x) x
  modelInfo_bear$prob <- function(x) x

  return(modelInfo_bear)
}



#' modelInfo for caret doing exp VI curve fitting
#'
#' @return the modelInfo(list)
#' @export
modelInfo_exp <- function() {
  modelInfo_exp <- getModelInfo("rlm")$rlm
  modelInfo_exp$label <- "exp nlsLM"
  modelInfo_exp$library <- "minpack.lm"
  modelInfo_exp$type <- "Regression"
  modelInfo_exp$parameters <- data.frame(
    parameter = c("parameter"),
    class = c("character"),
    label = c("parameter")
  )
  modelInfo_exp$grid <- function(x, y, len = NULL, search = "grid") {
    data.frame(parameter = "none")
  }
  modelInfo_exp$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    df <- data.frame(x = x, y = y)
    minpack.lm::nlsLM(y ~ exp(a + b * x), data = df, start = list(a = 1, b = 1), na.action = na.omit)
  }
  modelInfo_exp$predict <- function(modelFit, newdata, submodels = NULL) {
    newdata <- data.frame(x = newdata)
    predict(modelFit, newdata)
  }

  return(modelInfo_exp)
}



#' GauPro Gasussian Process modelinfo for caret
#'
#' @return
#' @export
#'
#' @examples
modelInfo_GP <- function() {
  modelInfo_GP <- list()
  modelInfo_GP$label <- "GauPro Gaussian Process"
  modelInfo_GP$library <- "GauPro"
  modelInfo_GP$loop <- NULL
  modelInfo_GP$type <- "Regression"
  modelInfo_GP$prob <- NULL
  modelInfo_GP$parameters <- data.frame(
    parameter = c("parameter"),
    class = c("character"),
    label = c("Parameter")
  )

  # param grid
  modelInfo_GP$grid <- function(x, y, len = NULL, search = "grid") data.frame(parameter = "none")

  # fit
  modelInfo_GP$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    dat <- data.frame(x = x)
    dat$.outcome <- y

    x_upper <- max(x)
    x_lower <- min(x)
    x_sd <- sd(x)

    GauPro(X = x, Z = y, parallel = TRUE)
  }

  # predict
  modelInfo_GP$predict <- function(modelFit, newdata, submodels = NULL) {
    modelFit$predict(newdata)
  }

  modelInfo_GP$tags <- c("Linear Regression", "Robust Model", "Accepts Case Weights")
  modelInfo_GP$sort <- function(x) x
  modelInfo_GP$prob <- function(modelFit, newdata, summodels = NULL) {
    modelFit$predict(newdata, se = TRUE)
  }

  return(modelInfo_GP)
}


#' modelInfo fore caret doing rf tuning both mtry and ntree
#'
#' @return the modelInfo(list)
#' @export
modelInfo_rf_mtry_ntree <- function() {
  modelInfo_rf <- getModelInfo("rf")$rf
  modelInfo_rf$grid <- function(x, y, len = NULL, search = "grid") {
    expand.grid(
      mtry = caret::var_seq(p = ncol(x), classification = is.factor(y), len = len),
      ntree = c(seq(50, 500, by = 50))
    )
  }
  modelInfo_rf$parameters <- data.frame(
    parameters = c("mtry", "ntree"),
    class = c("numeric", "numeric"),
    label = c("mtry", "ntree")
  )
  modelInfo_rf$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    randomForest::randomForest(x, y, mtry = param$mtry, ntree = param$ntree, ...)
  }
  modelInfo_rf$sort <- function(x) {
    x[order(x$ntree, x$mtry), ]
  }

  return(modelInfo_rf)
}
