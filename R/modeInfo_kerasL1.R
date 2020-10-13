#' modelInfo for Keras L1
#'
#' @return modelInfo
#' @export
modelInfo_kerasL1 <- function() {
  modelInfo <- caret::getModelInfo('mlpKerasDecay')$mlpKerasDecay
  # define sort ---------------------------------------------------------------------------------
  modelInfo$sort <- function(x) x[order(x$units01), ]

  # define parameters ---------------------------------------------------------------------------
  modelInfo$parameters <- data.frame(
    parameter = c("units01", "epochs"),
    class = c("numeric", "numeric"),
    label = c("#Hidden Units in Layer1", "epochs")
  )

  modelInfo$grid <- function(x, y, len = NULL, search = "grid") {
      expand.grid(units01 = (1:len)*2-1,
                  epochs = 200)
  }


  modelInfo$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
      require(dplyr)
      K <- keras::backend()
      K$clear_session()
      if (!is.matrix(x)) x <- as.matrix(x)

      model <- keras::keras_model_sequential()
      model %>%
        keras::layer_dense(
          units = param$units01,
          activation = "relu",
          kernel_initializer = keras::initializer_glorot_uniform(),
          kernel_regularizer = keras::regularizer_l2(),
          input_shape = ncol(x)
        )

      model %>%
        keras::layer_dense(
          units = 1,
          kernel_regularizer = keras::regularizer_l2(),
          activation = "linear"
          )

      model %>%
        keras::compile(
          loss = "mean_squared_error",
          optimizer = "adam"
          )


      model %>%
        keras::fit(
          x = x,
          y = y,
          batch_size = nrow(x),
          epochs = param$epochs,
          verbose = 0
          )

      if (last)  model <- keras::serialize_model(model)
      list(object = model)
  }

  return(modelInfo)
}
