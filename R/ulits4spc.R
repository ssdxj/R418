


#' ref to log(1/ref)
#'
#' @param spc
#'
#' @return
#' @export
#'
#' @examples
spc_2logReverse <- function(spc){

  ref <- spectra(spc)
  ref <- log(1/ref)
  spectra(spc) <- ref
  return(spc)
}

#' genearate Speclib obj from: csv file or data.frame
#'
#' @param input csv file or data.frame
#' @param idCol colname of id col ("PlotID")
#'
#' @return Speclib obj
#' @export
#'
#' @examples
spc_generator <- function(input, idCol = 'PlotID') {
  input_type <- class(input)
  if (inherits(input, 'character'))
    input <- read_csv(input)

  if (inherits(input, 'data.frame')) {
    df_meta <- dplyr::select(input, matches('^[[:alpha:]]')) %>%
      mutate(ID = .data[[idCol]]) %>%
      column_to_rownames('ID')

    df_spectra <- dplyr::select(input, matches('^[[:digit:]]'))
    wl <- names(df_spectra) %>% parse_number()

    spc <- speclib(spectra = as.matrix(df_spectra), wavelength = wl)
    SI(spc) <- df_meta

    idSpeclib(spc) <- rownames(df_meta)



    return(spc)
  } else{
    msg = sprintf('Input error: %s !!!', input_save)
    stop(msg)
  }
}


#' spc or spc_df to long df
#'
#' @param input speclib obj or df
#'
#' @return df (long df with "wl" and "reflect" column)
#' @export
#'
#' @examples
spc_2df_melt <- function(input){

  if (inherits(input, 'Speclib')) input <- spc_2df(input)

  input %>%
    pivot_longer(matches('^[[:digit:]]'), names_to = 'wl', values_to = 'reflect') %>%
    mutate(wl = parse_number(wl))
}


#' spc to long df with masked bands as NA for ggplot
#'
#' @param spc
#'
#' @return df (long df with "wl" and "reflect" column)
#' @export
#'
#' @examples
spc_2df_meltWithMask <- function(spc){
  if(!inherits(spc, 'Speclib')) stop('A Speclib is needed!!!')

  df <- spc_2df_melt(spc)
  masks_df <- mask(spc)

  if (is.null(masks_df)) {
    return (df )
  } else {
    wl <- df$wl
    reflect <- df$reflect
    for (i in 1:nrow(masks_df)) {
      reflect[wl >= masks_df$lb[i] & wl <=masks_df$ub[i] ] <- NA
    }

    df$reflect <- reflect
    return(df)
  }
}


#' spc bandwise cor
#'
#' @param spc the speclib object
#' @param yname meta colname
#'
#' @return df data.frame with colnames c(wl, colLabel)
#' @export
#'
#' @examples
spc_bandwiseCor <- function(spc, yname, colLabel){
  if(!inherits(spc, 'Speclib')) stop('A Speclib is needed!!!')

  spc_2df_meltWithMask(spc) %>%
    as_tibble() %>%
    dplyr::select(one_of(c(yname, 'wl', 'reflect'))) %>%
    # drop_na(reflect) %>%
    drop_na(.data[[yname]]) %>%
    group_by(wl) %>%
    summarise(!!(colLabel) := cor(.data[[yname]], .data$reflect, use = 'na.or.complete'),
              .groups = 'drop')
}



#' spc average
#'
#' @param spc
#'
#' @return
#' @export
#'
#' @examples
spc_reduceRepID <- function(spc){
  spc_2df_melt(spc) %>%
    group_by(FieldID, SampleDate, PlotID, SampleID, wl) %>%
    summarise(reflect = mean(reflect), .groups = 'drop') %>%
    pivot_wider(names_from = wl, values_from = reflect) %>%
    spc_fromDf()
}






#' Speclib obj to wide df
#'
#' @param spc the Speclib obj
#'
#' @return df (wide df)
#' @export
#'
#' @examples
spc_2df <- function(spc) {
  if(!inherits(spc, 'Speclib')) stop('A Speclib is needed!!!')
  df_meta <- SI(spc)
  df_spectra <- spectra(spc)
  colnames(df_spectra) <- wavelength(spc)
  df_spectra <- as_tibble(df_spectra)

  as_tibble(cbind(df_meta, df_spectra, stringsAsFactors = FALSE), .name_repaire = 'minimal')

}

#' Speclib obj to wide df (reflect col with 'B_' prefix)
#'
#' @param spc the Speclib obj
#'
#' @return df (wide df)
#' @export
#'
#' @examples
spc_2dfB <- function(spc){
  if(!inherits(spc, 'Speclib')) stop('A Speclib is needed!!!')
  df_meta <- SI(spc)
  df_spectra <- spectra(spc)
  colnames(df_spectra) <- paste('B', wavelength(spc), sep = '_')
  df_spectra <- as_tibble(df_spectra)

  as_tibble(cbind(df_meta, df_spectra), .name_repaire = 'minimal')
}


spc_2df4plot <- function(spc) {
  out <- NULL
  attri <- SI(spc)
  ref <- spectra(spc)

  # incase no attri
  if (ncol(attri) == 0) {
    out <- as_tibble(ref)
  } else {
    out <- as_tibble(cbind(attri, ref))
  }

  # handle colnames
  names(out) <- c(names(attri), hsdar::wavelength(spc))

  masked <- mask(spc)
  if(!is.null(maskVector)){
    for(rows in 1:nrow(masked)){
      lb <- masked$lb[rows]
      ub <- masked$ub[rows]

      NA_wl <- seq(lb+1, ub-1, by = spc@fwhm)
      names(NA_wl) <- NA_wl
      NA_ref <- map_df(NA_wl, ~rep(NA, times = nspectra(spc)))
      out <- c(out, NA_ref)
    }
  }

  return(out)
}


#' generate 'Speclib' obj from data.frame, reflectance selected by matches
#' colnames matches('^(\\d)+(\\.\\d+)?$'),wl determined by colnames of spectral data
#'
#' @param df data.frame of data
#' @param bands_reg with default '^(\\d)+(\\.\\d+)?$'
#'
#' @return Speclib obj
#' @export
spc_fromDf <- function(df, bands_reg = "^(\\d)+(\\.\\d+)?$") {
  # do select
  spectra <- dplyr::select(df, matches(bands_reg))
  attri <- dplyr::select(df, -matches(bands_reg)) %>% data.frame()

  # handle spc
  wl <- as.double(names(spectra))
  spc <- speclib(as.matrix(spectra), wl)
  SI(spc) <- as.data.frame(attri)

  return(spc)
}


#' rbind spc
#'
#' @param ... spc
#'
#' @return spc obj
#' @export
spc_rbind <- function(...) {
  out <- NULL

  spc_list <- list(...)
  wl_list <- map(spc_list, wavelength)
  wl_check <- map_lgl(wl_list, all.equal, wl_list[[1]]) %>% all()

  if (wl_check) {
    spc_df_list <- map(spc_list, spc_2df)
    spc_df <- do.call(rbind, spc_df_list)
    out <- spc_fromDf(spc_df)
    return(out)
  } else {
    stop("wavelength not match, stop!!")
  }
}


#' get inTrain spc by index vector
#'
#' @param spc: spc to be filtered
#' @param inTrain: inTrain index vector
#'
#' @return spc obj
#' @export
#'
spc_inTrain <- function(spc, inTrain) {
  spc_fromDf(spc_2df(spc)[inTrain, ])
}

#' melt the reflectance part of spc df
#'
#' @param x spc or spc_df
#' @param band_reg with default '^(\\d)+(\\.\\d+)?$'
#'
#' @return tibble
#' @export
#'
spc_melt <- function(x, band_reg = "^(\\d)+(\\.\\d+)?$") {
  if (is.speclib(x)) x <- spc_2df(x)

  x %>%
    gather(key = "wl", value = "reflect", matches(band_reg), convert = TRUE) %>%
    as_tibble()
}


#' Self use function. Average spc SampleID or PlotID or Treatment.
#'
#' @param spc: spc obj
#' @param by: SampleID/PlotID/Treatment/SampleDate
#' \enumerate{
#'   \item SampleID: group_by(spc_df_melt, FieldID, SampleDate, PlotID, SampleID, wl)
#'   \item PlotID: group_by(spc_df_melt, FieldID, SampleDate, PlotID, wl)
#'   \item Treatment: group_by(spc_df_melt, FieldID, SampleDate, Treatment, wl)
#'   \item SampleDate: group_by(spc_df_melt, FieldID, SampleDate, wl)
#' }
#'
#' @return spc obj
#' @export
#'
spc_ave <- function(spc, by = "SampleID") {
  # melt
  spc_df_melt <- spc_melt(spc)
  out <- NULL

  if (by == "SampleID") {
    out <- group_by(spc_df_melt, FieldID, SampleDate, PlotID, SampleID, wl)
  } else if (by == "PlotID") {
    out <- group_by(spc_df_melt, FieldID, SampleDate, PlotID, wl)
  } else if (by == "Treatment") {
    out <- group_by(spc_df_melt, FieldID, SampleDate, Treatment, wl)
  } else if (by == "SampleDate") {
    out <- group_by(spc_df_melt, FieldID, SampleDate, wl)
  } else {
    stop("Error in parameter value!!!")
  }

  # handle the longitude and latitude
  if (all(c("longitude", "latitude") %in% names(df))) {
    out <- out %>%
      summarise(
        longitude = mean(longitude, na.rm = TRUE),
        latitude = mean(latitude, na.rm = TRUE),
        reflect = mean(reflect, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      spread(wl, reflect)
  } else {
    out <- out %>%
      summarise(reflect = mean(reflect, na.rm = TRUE)) %>%
      ungroup() %>%
      spread(wl, reflect)
  }

  out <- spc_fromDf(out)
  return(out)
}


#' Bandwise cor with responser
#'
#' @param spc  spc
#' @param biochemphy name of responsor in SI
#'
#' @return df with col ('wl', 'estimate', 'p.value')
#' @export
spc_cor <- function(spc, biochemphy) {
  ref <- spectra(spc)
  wl <- wavelength(spc)
  colnames(ref) <- wl
  reps <- SI(spc)[[biochemphy]]

  apply(ref, 2, cor.test, reps, na.rm = TRUE) %>%
    map_df(function(fit){fit[c('estimate', 'p.value')]}, .id = 'wl') %>%
    mutate(wl = parse_double(wl))
}


#' wrapper of \code{\link{spc_cor}} doing bandwise band cor with responser
#' grouped by stage
#'
#' @param stageValue levels of stage in SI
#' @param spc spc
#' @param biochemphy name of reponser in SI
#'
#' @return df
#' @export
spc_cor_stage <- function(stageValue, spc, biochemphy) {
  if (stageValue == "full") {
    spc_cor(spc, biochemphy)
  } else {
    # incase no records for specific growth stage
    df <- spc_2df(spc) %>% dplyr::filter(stage %in% stageValue)

    # if so return NULL
    if(nrow(df) == 0) return(NULL)

    # else
    spc_cor(spc_fromDf(df), biochemphy)
  }
}


#'  change the reflectance of masked bands to NA
#'
#' @param input spc or melt spc_df
#' @param mask  NA or vector
#'
#' @return df
#' @export
spcdf_mask2NA <- function(input, mask = NA){
  if(inherits(input, 'Speclib')) input <- spc_melt(input)

  # band within mask change to NA
  # thus can be easily handled by ggplot2
  while(length(mask) >= 2){
    lower <- mask[1]
    upper <- mask[2]
    mask <- mask[c(-1, -2)]

    # help function to use mutate
    foo <- function(wl, reflect, lower, upper){
      flag <- wl >= lower & wl <= upper
      reflect[flag] <- NA
      return(reflect)
    }

    input <- mutate(input, reflect = foo(wl, reflect, lower, upper))
  }

  return(input)

}
