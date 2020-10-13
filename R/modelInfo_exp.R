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
