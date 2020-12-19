# spc transform -----------------------------------------------------------

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



# spc filter --------------------------------------------------------------




# spc calc ----------------------------------------------------------------



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
    options(stringsAsFactors = FALSE)
    spc_df <- do.call(rbind, spc_df_list)
    out <- spc_fromDf(spc_df)
    return(out)
  } else {
    stop("wavelength not match, stop!!")
  }
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
        reflect = mean(reflect, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      ungroup() %>%
      spread(wl, reflect)
  } else {
    out <- out %>%
      summarise(reflect = mean(reflect, na.rm = TRUE), .groups = 'drop') %>%
      spread(wl, reflect)
  }

  out <- spc_fromDf(out)
  return(out)
}



# spc from/to df ----------------------------------------------------------

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
    msg <- sprintf('Input error: %s !!!', input_type)
    stop(msg)
  }
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
  if('wl' %in% names(df)) {
    df <- pivot_wider(df, names_from = 'wl',
                                            values_from = 'reflect')
  }


  # do select
  df_spectra <- dplyr::select(df, matches(bands_reg))
  df_meta <- dplyr::select(df, -matches(bands_reg))
  wl <- parse_double(names(df_spectra))

  # handle spc
  spc <- speclib(as.matrix(df_spectra), wl)
  SI(spc) <- as.data.frame(df_meta)

  return(spc)
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
  if(inherits(spc, 'data.frame')) {
    df <- mutate(spc, .ID = nrow(spc))
    return(df)

  }
  if(!inherits(spc, 'Speclib')) stop('A Speclib is needed!!!')

  df_meta <- SI(spc)
  df_spectra <- spectra(spc)

  colnames(df_spectra) <- wavelength(spc)
  df_spectra <- as_tibble(df_spectra)

  bind_cols(df_meta, df_spectra) %>%
    as_tibble() %>%
  mutate(.ID = 1:nrow(.))
}


#' Title
#'
#' @param spc
#'
#' @return
#' @export
#'
#' @examples
spc_2df_withMask <- function(spc){
  if (inherits(spc, 'data.frame')) return(spc)
  df <- spc_2df(spc)

  masks <- mask(spc)
  if(!is.null(masks)){
    wl_missing <- c()
    # geneate the wavelength to fill in by step of 5nm
    for(i in 1:nrow(masks)){
      wl_missing <- c(wl_missing, seq(masks$lb[i]+1,masks$ub[i]-1, by = 5))
    }
    # add NA to the df
    for(band in wl_missing){
      df[[as.character(band)]] <- NA
    }
  }

  return(df)
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

  bind_cols(df_meta, df_spectra) %>% as_tibble()

}



#' spc or spc_df to long df
#'
#' @param input speclib obj or df
#'
#' @return df (long df with "wl" and "reflect" column)
#' @export
#'
#' @examples
spc_melt <- function(input){
  if (inherits(input, 'Speclib')) input <- spc_2df(input)
  input %>%
    mutate(.ID = 1:nrow(.)) %>%
    pivot_longer(matches('^[[:digit:]]'), names_to = 'wl', values_to = 'reflect') %>%
    mutate(wl = parse_number(wl))
}



# spc cor -----------------------------------------------------------------

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
    summarise(!!(colLabel) := cor(.data[[yname]], .data[['reflect']], use = 'na.or.complete'),
              .groups = 'drop')
}



# vegindex addon ----------------------------------------------------------


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ssdxj_rededge <- function(x){
  n_sgolay <- floor((25/mean(x@fwhm))/2)*2+1
  if (n_sgolay < 5)
    n_sgolay <- 5

  D1 <- derivative.speclib(x,method="sgolay",m=1, n=n_sgolay)
  D2 <- derivative.speclib(D1,method="sgolay",m=1, n=n_sgolay)

  RedEdge_data <- as.data.frame(t(as.matrix(sapply(c(1:nspectra(x)),
                                                   FUN = ssdxj_rededge_apply, spectra(x), D1, D2), ncol = 6)))


  row.names(RedEdge_data) <- idSpeclib(x)
  names(RedEdge_data) <- c("R0","l0","Rp","lp","Rs","ls")

  #  if (round)
  #  {
  #    RedEdge_data[,1] <- round(RedEdge_data[,1], 0)
  #    RedEdge_data[,2] <- round(RedEdge_data[,2], 0)
  #    RedEdge_data[,3] <- round(RedEdge_data[,3], 0)
  #  }

  return(as_tibble(RedEdge_data))
}

#' Title
#'
#' @param i
#' @param x
#' @param D1
#' @param D2
#'
#' @return
#' @export
#'
#' @examples
ssdxj_rededge_apply <- function(i, x, D1, D2) {
  i <- i[1]

  # l0: wavelength of the minimum reflectance in the red spectrum
  # R0: reflectance at l0
  tmp <- wavelength(D2) >= 660 & wavelength(D2) <= 700
  R0 <- min(x[i,tmp],na.rm=TRUE)
  l0 <- wavelength(D2)[tmp]
  l0 <- l0[which.min(abs(R0 - x[i,tmp]))]


  # lp: wavelength of the inflection point
  # tmp <- wavelength(D2) >= 700 & wavelength(D3) <= 750
  tmp <- wavelength(D2) >= 680 & wavelength(D2) <= 750
  tmp2 <- spectra(D1)[i,]
  tmp2[!tmp] <- -99999.9
  lp <- which.max(tmp2) # index of max D1 in c(700, 750)
  Rp <- x[i,lp]
  lp <- wavelength(D2)[lp]

  # ls: wavelength of the reflectance shoulder
  tmp <- wavelength(D2) > lp & wavelength(D2) < 900 # wl length flag vector
  tmp2 <- sign(spectra(D2)[i,tmp]) # d2 sign (1, 0, -1)
  tmp3 <- tmp2[-c(1,2)]*tmp2[-c(length(tmp2)-1,length(tmp2))] # a(i) * a(i+2)
  tmp3 <- c(FALSE,tmp3==-1,FALSE) # D2 sign of a(i) and a(i+2) are not equal IS TRUE
  tmp4 <- wavelength(D2)[tmp]
  tmp3 <- tmp4[tmp3]
  ls <- tmp3[1]
  if (is.finite(ls))
  {
    Rs <- x[i,wavelength(D2)==ls]
  } else {
    Rs <- NA
  }
  return(c(R0,l0,Rp,lp,Rs,ls))

}


# depresed ----------------------------------------------------------------

# spc_2df4plot <- function(spc) {
#   out <- NULL
#   attri <- SI(spc)
#   ref <- spectra(spc)
#
#   # incase no attri
#   if (ncol(attri) == 0) {
#     out <- as_tibble(ref)
#   } else {
#     out <- as_tibble(cbind(attri, ref))
#   }
#
#   # handle colnames
#   names(out) <- c(names(attri), hsdar::wavelength(spc))
#
#   masked <- mask(spc)
#   if(!is.null(maskVector)){
#     for(rows in 1:nrow(masked)){
#       lb <- masked$lb[rows]
#       ub <- masked$ub[rows]
#
#       NA_wl <- seq(lb+1, ub-1, by = spc@fwhm)
#       names(NA_wl) <- NA_wl
#       NA_ref <- map_df(NA_wl, ~rep(NA, times = nspectra(spc)))
#       out <- c(out, NA_ref)
#     }
#   }
#
#   return(out)
# }







