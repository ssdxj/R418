#' Title
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
