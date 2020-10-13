#' Hyperspectral raster calculation
#'
#' @param fpath_in input raster file path(in ENVI format)
#' @param fpath_out output file path(in .tif format)
#' @param wl input raster file wavelength
#' @param nbands nbands for output
#' @param fun do function, take hsdar::Speclib as first input param, each spectrum
#' is one pixel of fpath_in
#' @param format parameter for writeStart to control output format, default raster.
#' @param noData value for NA
#' @param ... futher param for fun
#'
#' @return Hyperspectral image obj(save output to fpath_out inviable)
#' @export

raster_calc <- function(fpath_in, fpath_out, wl, nbands, fun, format = "raster",
                        noData = -9999, ...) {
  # incase
  library(hsdar)

  # input raster obj
  ra <- HyperSpecRaster(fpath_in, wl)
  NAvalue(ra) <- noData

  # output handler
  res <- writeStart(
    x = ra, filename = path_out, format = "raster", nl = nbands,
    overwrite = TRUE
  )

  # brick size
  tr <- blockSize(ra)

  # handle progress bar
  pb <- txtProgressBar(min = 0, max = tr$n, style = 3)

  # start loop
  for (i in 1:tr$n) {
    # subtract image block into spc
    v <- getValuesBlock(ra, row = tr$row[i], nrows = tr$nrows[i])

    # main routine
    out <- fun(v, ...)

    # write to out file
    res <- raster::writeValues(res, out, tr$row[i])

    # handle progress bar
    setTxtProgressBar(pb, i)
  }

  res <- writeStop(res)

  # handle progress bar
  close(pb)

  # return result
  return(res)
}
