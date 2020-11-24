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
