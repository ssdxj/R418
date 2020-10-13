#' modelInfo for Keras
#'
#' @return modelInfo
#' @export
modelInfo_kerasL2 <- function() {
  modelInfo <- list()
  modelInfo$label <-
    "keras with one hidden layer tuning units"
  modelInfo$library <- "keras"
  modelInfo$loop <- NULL
  modelInfo$type <- "Regression"
  modelInfo$prob <- function() {
  }
  modelInfo$varImp <- NULL
  modelInfo$check <- function(pkg) {
    testmod <- try(keras::keras_model_sequential(), silent = TRUE)
    if (inherits(testmod, "try-error")) {
      stop(
        "Could not start a sequential model. ",
        "`tensorflow` might not be installed. ",
        "See `?install_tensorflow`.",
        call. = FALSE
      )
    }
    TRUE
  }


  # define sort ---------------------------------------------------------------------------------
  modelInfo$sort <-
    function(x)
      x[order(x$units01, x$units02), ]

  # define parameters ---------------------------------------------------------------------------
  modelInfo$parameters <- data.frame(
    parameter = c("units01", "units02", "epochs"),
    class = c("numeric", "numeric", "numeric"),
    label = c("#Hidden Units in Layer1", "#Hidden Units in Layer1", "#Epochs")
  )

  modelInfo$grid <-
    function(x, y, len = NULL, search = "grid") {
      out <- expand.grid(
        units01 = ((1:len) * 2) - 1,
        units02 = c(0, ((1:len) * 2) - 1),
        epochs = 100
      )

      out
    }


  modelInfo$fit <-
    function(x,
                 y,
                 wts,
                 param,
                 lev,
                 last,
                 classProbs,
                 ...) {
      require(dplyr)
      K <- keras::backend()
      K$clear_session()
      if (!is.matrix(x)) {
        x <- as.matrix(x)
      }
      y <- unlist(y)


      model <- keras::keras_model_sequential()
      model %>%
        keras::layer_dense(
          units = param$units01,
          kernel_initializer = keras::initializer_glorot_normal(seed = 666),
          activation = "relu",
          input_shape = ncol(x)
        )

      if (param$units02 != 0) {
        model %>%
          keras::layer_dense(
            units = param$units02,
            kernel_initializer = keras::initializer_glorot_normal(seed = 666),
            activation = "relu"
          )
      }


      model %>%
        keras::layer_dense(units = 1, activation = "linear")

      model %>%
        keras::compile(loss = "mse", optimizer = "adam")


      model %>%
        keras::fit(
          x = x,
          y = y,
          batch_size = floor(nrow(x) / 3),
          epochs = param$epochs,
          verbose = 0
        )

      if (last) {
        model <- keras::serialize_model(model)
      }
      list(object = model)
    }


  # define predict ------------------------------------------------------------------------------
  modelInfo$predict <-
    function(modelFit, newdata, submodels = NULL) {
      if (inherits(modelFit$object, "raw")) {
        modelFit$object <- keras::unserialize_model(modelFit$object)
      }

      newdata <- as.matrix(newdata)

      out <- predict(modelFit$object, newdata)

      ## check for model type
      if (ncol(out) == 1) {
        out <- out[, 1]
      } else {
        stop("caret predict error!!!")
      }

      out
    }

  return(modelInfo)
}
