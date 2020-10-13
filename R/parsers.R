# file parser -------------------------------------------------------------

#' parser of ForceA
#'
#' @param fpath file path
#'
#' @return parsed df
#' @export
parser_ForceA <- function(fpath) {
  # read
  data <- readLines(fPath, skipNul = TRUE)

  # parse
  data_sub <- grep("^\\d{4}/\\d{2}/\\d{2};", data, value = TRUE)
  data_sub_parsed <- lapply(data_sub, function(x) {
    unlist(str_split(x, ";"))
  })
  data_sub_parsed_df <- do.call(rbind.data.frame, data_sub_parsed)
  colnames(data_sub_parsed_df) <- c(
    "SampleDate", "Time", "longitude", "latitude", "altitude", "sat_qual",
    "temp", "group", "measure", "side", "Chl", "Flav", "Anth", "NBI", "Calib"
  )

  # select
  out <- data_sub_parsed_df[, -c(1, 3:6, 15)]

  # rename
  fname_old <- tools::file_path_sans_ext(fPath)
  fname_new <- paste(fname_old, "_parsed.csv", sep = "")

  # write
  write.csv(out, fname_new, row.names = FALSE)

  return(out)
}

#' parser of SunScan file
#'
#' @param  fpath fpath of SunScan file
#'
#' @return parsed df
#' @export
parser_sunScan <- function(fpath) {
  # read
  data <- readLines(fPath)

  # parse
  data_sub <- grep("^\\d{2}:\\d{2}:\\d{2}\\t\\d{1,2}", data, value = TRUE)
  data_sub_split <- lapply(data_sub, function(x) {
    unlist(str_split(x, "\t"))
  })
  df <- do.call(rbind.data.frame, data_sub_split)
  names(df) <- c(
    "Time", "ID1", "ID2", "Transmitted", "Spread", "Incident",
    "Beam", "ZenithAngle", "LAI", "Note"
  )

  # select
  out <- df[, -c(6, 7, 10)]

  # rename
  fname_old <- tools::file_path_sans_ext(fPath)
  fpath_new <- paste(fname_old, "_parsed.csv", sep = "")

  # write
  write.csv(out, fpath_new, row.names = FALSE)

  return(out)
}

#' parser of single HR1024i sig file
#'
#' @param fpath: fpath of sig file. The Field, SampleDate, isCanopy, PlotID,
#' SampleID and RepID are extracted from the fname.
#' So fname must in format 'YNSD_20140712_C_01_01_00_moc_resamp.sig'
#' which is 'FieldID_SampleDate_isCanopy_PlotID_SampleID_RepID_others.sig'
#' @return parsed content(list)
#' @export
#'
parser_sig <- function(fpath) {
  # list for store result
  out <- list()

  # read file into string vector
  sig_raw <- readLines(fpath)

  # grep the meta part
  sig_meta <- grep("^[[:alpha:]]+", sig_raw, value = TRUE)
  # grep the data part
  sig_main <- grep("^[[:digit:]]+", sig_raw, value = TRUE)

  # handle the data part
  sig_mainTable <- read.table(text = sig_main)

  # get wl
  out$wl <- sig_mainTable[[1]]
  # get ref
  out$ref <- sig_mainTable[[4]] / 100.0

  # handle the meta part
  # name from file
  out$name1 <- tools::file_path_sans_ext(toupper(basename(fpath)))

  # longitude
  char_greped <- grep("^longitude", sig_meta, value = TRUE)
  longitude <- trimws(strsplit(strsplit(char_greped, "=")[[1]][2], ",")[[1]][2])

  # latitude
  char_greped <- grep("^latitude", sig_meta, value = TRUE)
  latitude <- trimws(strsplit(strsplit(char_greped, "=")[[1]][2], ",")[[1]][2])

  # case not latitude and longitude for in indoor svc
  # if(longitude == '' | latitude == '') longitude = latitude = 0
  if (longitude == "" | latitude == "") {
    longitude <- latitude <- NA
  } else {
    # cheap check longitude and latitude length
    # some files have false code in this region
    if (nchar(longitude) != 11 | nchar(latitude) != 10) {
      longitude <- NA
      latitude <- NA
    } else {
      longitude_p1 <- as.double(str_sub(longitude, 1, 3))
      longitude_p2 <- as.double(str_sub(longitude, 4, 10))
      longitude_p3 <- ifelse(str_sub(longitude, 11, 11) == "E", 1, -1)
      longitude <- (longitude_p1 + longitude_p2 / 60) * longitude_p3

      latitude_p1 <- as.double(str_sub(latitude, 1, 2))
      latitude_p2 <- as.double(str_sub(latitude, 3, 9))
      latitude_p3 <- ifelse(str_sub(latitude, 10, 10) == "N", 1, -1)
      latitude <- (latitude_p1 + latitude_p2 / 60) * latitude_p3
    }
  }

  out$longitude <- longitude
  out$latitude <- latitude

  # is moc
  char.greped <- grep("^factors", sig_meta, value = TRUE)
  # out$isMOC <- ifelse(grep('Overlap: Remove', char.greped), 'MOC', 'NOMOC')

  # calc FieldID from file name
  chr_splited <- unlist(strsplit(out$name1, "_"))
  out$FieldID <- chr_splited[1]
  out$SampleDate <- chr_splited[2]
  out$PlotID <- paste("P", chr_splited[4], sep = "")
  out$SampleID <- paste("S", chr_splited[5], sep = "")
  out$RepID <- paste("R", chr_splited[6], sep = "")

  # nbands
  out$nbands <- nrow(sig_mainTable)

  return(out)
}

#' Wrapper of \code{\link{parser_sig}}, to parse all sig files in a dir.
#'
#' @param dir: dir contains sig files
#' @return list of parsed content(list)
#' @references parser_sig()
#' @export
#'
parser_sigs <- function(dir) {
  # walk the input dir
  fpaths <- list.files(dir, pattern = "*.sig$", full.names = TRUE)
  nfiles <- length(fpaths)

  # for store out
  out <- list()
  for (i in 1:nfiles) {
    sig_parsed <- parser_sig(fpaths[i])
    out[[sig_parsed$name1]] <- sig_parsed
  }

  return(out)
}

#' generate spc obj form parsed contents
#'
#' @param dat result from parser_sigs
#' @return df df of parsed result
#' @export
#'
spc_from_sigs <- function(dat) {
  nspectral <- length(dat)

  # check bands match
  wl <- lapply(dat, function(x) {
    x$wl
  })
  for (i in 1:nspectral) {
    if (sum((wl[[i]] - wl[[1]])^2) != 0) {
      stop(paste(
        "The wavelength of :", names(dat)[i],
        "does not match the first one!!!"
      ))
    }
  }

  # determin wavelength(use the first one)
  wl <- wl[[1]]

  # determin reflectance data.frame
  # print("rbinding reflectance ...")
  ref_list <- lapply(dat, function(x) {
    x$ref
  })
  ref_df <- do.call(rbind, ref_list)

  # determin meta data.frame
  meta_list <- lapply(dat, function(x) {
    x[which(names(x) %not in% c("wl", "ref"))]
  })
  meta_df <- do.call(rbind.data.frame, meta_list)

  # create speclib
  # spc <- createspeclib(spectra = ref.df, wavelength = wl)
  spc <- speclib(ref_df, wl)
  SI(spc) <- meta_df

  return(spc)
}
