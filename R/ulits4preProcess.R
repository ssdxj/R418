# inTrain index -----------------------------------------------------------

#' get inTrain vector by ordered reponse value
#'
#' @param input numeric vector(normally the response)
#' @param step 1 out of step in Test
#' @param start one of 1:step where the first Test data is.
#'
#' @return index vector for inTrain
#' @export
get_inTrain_respOrder <- function(input, step = 4, start = 2) {
  x <- input
  x_ordered <- order(x)
  inTest <- seq(start, length(x) - 1, by = step)
  out <- x_ordered[-inTest] %>% sort()

  return(out)
}


#' get inTrain vector by ordered reponse value(within each group)
#'
#' @param input numeric vector(normally the response)
#' @param group str vector(ex: stage, year)
#' @param step 1 out of step in Test
#' @param start one of 1:step where the first Test data is.
#'
#' @return index vector for inTrain
#' @export
get_inTrain_respOrder_withGroup <- function(input, group, step = 4, start = 2) {
  idx <- 1:length(input)
  inTrain <- list()

  groupValues <- unique(group)
  for (gv in groupValues) {
    flag <- group == gv
    input_sub <- input[flag]
    idx_sub <- idx[flag]

    input_sub_order <- order(input_sub)

    input_sub_orderd <- input_sub[input_sub_order]
    idx_sub_orderd <- idx_sub[input_sub_order]

    inTrain_sub <- idx_sub_orderd[-seq(start, length(input_sub), by = step)]
    inTrain[[gv]] <- inTrain_sub
  }

  inTrain <- unlist(inTrain, use.names = FALSE) %>% sort()


  return(inTrain)
}



#' Title get inTrain index by full dataset and Test dataset
#'
#' @param input_full full spc (or SI)
#' @param input_Test  test spc (or SI)
#'
#' @return inTrain index vector
#' @export
get_inTrain_fromTestSPC <- function(input_full, input_Test) {
  if (is.speclib(input_full)) input_full <- SI(input_full)
  if (is.speclib(input_Test)) input_Test <- SI(input_Test)
  if (!is.data.frame(input_full) | !is.data.frame(input_Test)) {
    stop("input param error!!!")
  }

  input_full$index <- 1:nrow(input_full)

  out <- anti_join(input_full, input_Test, by = by) %>%
    dplyr::select(index) %>%
    unlist(use.names = FALSE)

  return(out)
}


# cross-validation index --------------------------------------------------


#' Create repeated n fold corssvalidation index
#'
#' @param x numerical vector(normally the inTrain part of response vector)
#' @param fold n folds
#' @param times n times
#'
#' @return list of cv index
#' @export
get_indexCV <- function(x, folds = 10, times = 5) {
  index <- map(3^(1:times - 3), function(seed) {
    set.seed(seed)
    createFolds(x, k = folds, returnTrain = TRUE)
  })
  names(index) <- paste("Round", 1:times, sep = "")
  index <- unlist(index, recursive = FALSE)

  return(index)
}


# prepare for train -------------------------------------------------------

#' prepare obj (list(inTrain, spc_full, indexCV)) for training workflow
#' \describe{
#' \item{inTrain}{\code{\link[G407]{get_inTrain_respOrder_withGroup}}}
#' \item{spc_full}{the spc}
#' \item{indexCV}{\code{\link[G407]{get_indexCV}}}
#' }
#'
#' @param spc ful spc
#' @param biochemphy name of response in SI
#' @param group  name of group in SI
#' @param times for cv
#' @param folds for cv
#' @param isSplit whether inTrain/Test split
#'
#' @return list(inTrain, spc_full, indexCV)
#'
#' @export
prepare_obj4wf <- function(spc, biochemphy, group, folds, times, isSplit) {
  if (isSplit) {
    # get inTrain index
    y <- SI(spc)[[biochemphy]]
    group <- SI(spc)[[group]]
    inTrain <- get_inTrain_respOrder_withGroup(y, group)

    # CV index
    y <- y[inTrain]
    indexCV <- get_indexCV(y, folds = folds, times = times)

    # out
    out <- list(
      inTrain = inTrain,
      spc_full = spc,
      indexCV = indexCV
    )
  } else {
    y <- SI(spc)[[biochemphy]]
    indexCV <- get_indexCV(y, folds = folds, times = times)
    out <- list(
      spc_full = spc,
      indexCV = indexCV
    )
  }

  return(out)
}
