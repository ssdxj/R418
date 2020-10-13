#' wrapper for do predict of SAE models and turn result to spc
#'
#' @param model keras(SAE) model
#' @param spc input spc obj
#' @param wl
#' \describe{
#'  \item{NULL}{for encoder model, use 1:n channel for wavelength}
#'  \item{nmerical vector}{for decoder model}
#' }
#'
#' @return spc obj
#' @export
predict_keras_spc <- function(model, spc, wl = NULL) {
  # spc to matrix
  newdata <- spectra(spc)

  # do predict
  pred <- predict(model, newdata)

  # if encoded use 1:n channel, else decoded use give wl
  if (is.null(wl)) wl <- 1:ncol(pred)

  # matrix to spc
  new_spc <- speclib(pred, wl)

  # give attri back
  attri <- SI(spc)
  if (!is.null(attri)) SI(new_spc) <- attri

  # final
  return(new_spc)
}


# custom layers -------------------------------------------------------------------------------


# Keras tied layer class
DenseTransposeTiedLayer <- R6::R6Class("KerasLayer",
  inherit = KerasLayer,
  public = list(
    output_dim = NULL,
    tied_to = NULL,
    tie_weights = NULL,
    initialize = function(output_dim, tied_to) {
      self$output_dim <- output_dim
      self$tied_to <- tied_to
    },
    call = function(x, mask = NULL) {
      self$tie_weights <- self$tied_to$weights
      output <- k_dot(x, k_transpose(self$tie_weights[[1]]))
      if (!is.null(self$activation)) {
        output <- self$activation(output)
      }
      return(output)
    },
    compute_output_shape = function(input_shape) {
      list(input_shape[[1]], self$output_dim)
    }
  )
)

#' warpper to generate keras tied layer
#'
#' @param object tensor obj
#' @param output_dim output dim
#' @param tied_to  layer tie to
#'
#' @return layer
#' @export
#'
layer_DenseTransposeTiedLayer <- function(object, output_dim, tied_to) {
  create_layer(DenseTransposeTiedLayer, object, list(
    output_dim = as.integer(output_dim),
    tied_to = tied_to
  ))
}
