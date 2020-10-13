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
