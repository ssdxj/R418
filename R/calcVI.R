#' vegetion index calculation handler
#'
#' @param index: name of index, used to match function
#' @param spc: the speclib obj
#' @return the index value (vector)
#' @export
ssdxj_vegindex <- function(index, spc, weighted = FALSE, ...) {
  param <- list(...)

  nir <- get_reflectance(spc, wavelength = 830, weighted = weighted)
  red <- get_reflectance(spc, wavelength = 670, weighted = weighted)
  green <- get_reflectance(spc, wavelength = 550, weighted = weighted)


  # RED/NIR -----------------------------------------------------------------

  if(index == 'NDVI'){
    out <- vegindex(spc, index)

  }else if(index == 'MSR'){
    # Chen1996
    out <- (nir/red-2)/(sqrt(nir/red) + 1)

  } else if(index == 'DVI') {
    out <- nir -red

  } else if(index == 'RDVI'){
    # Roujean and Breon (1995)
    out <- vegindex(spc, index)

  } else if(index == 'IDVI'){
    # He, Li, Craig A. Coburn, Zhi-Jie Wang, Wei Feng, and Tian-Cai Guo. 2018. “Reduced Prediction Saturation and View Effects for Estimating the Leaf Area Index of Winter Wheat.” IEEE Transactions on Geoscience and Remote Sensing, 1–16. https://doi.org/10.1109/TGRS.2018.2868138.
    out <- (1+(nir-red))/(1-(nir-red))

  } else if(index %in% c('Green_NDVI', 'Green NDVI')){
    # Gitelson, Anatoly A., Yoram J. Kaufman, and Mark N. Merzlyak. 1996. “Use of a Green Channel in Remote Sensing of Global Vegetation from EOS-MODIS.” Remote Sensing of Environment 58 (3): 289–98. https://doi.org/10.1016/S0034-4257(96)00072-7.
    out <- vegindex(spc, 'Green NDVI')


    # soil adjust -------------------------------------------------------------
  } else if(index == 'MSAVI'){
    # Qi1994
    out <- vegindex(spc, index)

  } else if(index == 'OSAVI'){
    #  Rondeaux et al. 1996
    out <- vegindex(spc, index)

    # CI ----------------------------------------------------------------------
  } else if(index == 'CI_rededge'){

    out <- vegindex(spc, 'CI2')

  } else if(index == 'CI_green'){
    # LAI Gitelson et al,
    out <- nir / green - 1
  } else if(index == 'CI_rededge743'){
    # Gitelson et al, (2003)
    re <- get_reflectance(spc, wavelength = 743, weighted = FALSE)
    out <- nir / re - 1
  } else if(index == 'CI_rededge705'){
    # Gitelson et al, (2003)
    re <- get_reflectance(spc, wavelength = 705, weighted = FALSE)
    out <- nir / re - 1


    # Sims and Gamon (2002)
  } else if(index %in% c('mSR', 'mSR705','mNDVI', 'mND705')){
    # Sims, Daniel A., and John A. Gamon. 2002. “Relationships between Leaf Pigment Content and Spectral Reflectance across a Wide Range of Species, Leaf Structures and Developmental Stages.” Remote Sensing of Environment 81 (2–3): 337–54. https://doi.org/10.1016/S0034-4257(02)00010-X.
    out <- vegindex(spc, index)

    # Datt 1999 ---------------------------------------------------------------
  } else if(index %in% c('Datt', 'Datt2', 'Datt3')){
    # Datt, B. 1999. “Visible/near Infrared Reflectance and Chlorophyll Content in Eucalyptus Leaves.” International Journal of Remote Sensing 20 (14): 2741–59. https://doi.org/10.1080/014311699211778.
    out <- vegindex(spc, index)

    #   CARI varities -----------------------------------------------------------

  } else if(index == 'CARI'){
    # Kim, M. S. 1994. “The Use of Narrow Spectral Bands for Improving Remote Sensing Estimations of Fractionally Absorbed Photosynthetically Active Radiation (Fapar).” Master Thesis, University of Maryland.
    out <- vegindex(spc, index)

  } else if(index %in% c('MCARI', 'MCARI/OSAVI')){
    # Daughtry, C.S.T, C.L Walthall, M.S Kim, E.Brown de Colstoun, and J.E McMurtrey. 2000. “Estimating Corn Leaf Chlorophyll Concentration from Leaf and Canopy Reflectance.” Remote Sensing of Environment 74 (2): 229–39. https://doi.org/10.1016/S0034-4257(00)00113-9.
    out <- vegindex(spc, index)

  } else if(index %in% c('TCARI', 'TCARI/OSAVI')){
    # Haboudane, Driss, John R Miller, Nicolas Tremblay, Pablo J Zarco-Tejada, and Louise Dextraze. 2002. “Integrated Narrow-Band Vegetation Indices for Prediction of Crop Chlorophyll Content for Application to Precision Agriculture.” Remote Sensing of Environment 81 (2–3): 416–26.
    # Haboudane et al. (2002) Chl
    out <- vegindex(spc, index)


    # Wu2008 ------------------------------------------------------------------

  } else if(index %in% c('MCARI2', 'TCARI2', 'MCARI2/OSAVI2', 'TCARI2/OSAVI2')){
    # Wu, Chaoyang, Zheng Niu, Quan Tang, and Wenjiang Huang. 2008. “Estimating Chlorophyll Content from Hyperspectral Vegetation Indices: Modeling and Validation.” Agricultural and Forest Meteorology 148 (8–9): 1230–41.

    out <- vegindex(spc, index)


  } else if(index == 'MCARI2_LAI'){

    # Haboudane, Driss, John R. Miller, Elizabeth Pattey, Pablo J. Zarco-Tejada, and Ian B. Strachan. 2004. “Hyperspectral Vegetation Indices and Novel Algorithms for Predicting Green LAI of Crop Canopies: Modeling and Validation in the Context of Precision Agriculture.” Remote Sensing of Environment 90 (3): 337–52. https://doi.org/10.1016/j.rse.2003.12.013.

    r800 <- get_reflectance(spc, wavelength = 800, weighted = weighted)
    r550 <- get_reflectance(spc, wavelength = 550, weighted = weighted)
    r670 <- get_reflectance(spc, wavelength = 670, weighted = weighted)

    num <- 1.5 * (2.5 * (r800- r670) - 1.3 * (r800- r550))
    deo <- sqrt(2 * (r800 + 1)^2 - (6 * r800 - 5 * sqrt(r670)) - 0.5)
    out <- num / deo


    # } else if(index == 'MCARI2_wu2010'){
    #   r750 <- get_reflectance(spc, wavelength =  750, weighted = weighted)
    #   r705 <- get_reflectance(spc, wavelength = 705, weighted = weighted)
    #   num <- 1.5 * (2.5 * (r750 - r705) - 1.3 * (r750 - green))
    #   deo <- sqrt(2 * (r750 + 1)^2 + 6 * r750 + 5 * sqrt(r705) - 0.5)
    #   out <- num / deo




    # TVI varities ------------------------------------------------------------
  } else if(index == 'TVI'){
    # Broge and Leblanc 2001 LAI
    out <- vegindex(spc, index)

  } else if(index == 'TGI'){
    # Hunt, E. Raymond, Paul C. Doraiswamy, James E. McMurtrey, Craig S T Daughtry, Eileen M. Perry, and Bakhyt Akhmedov. 2012. “A Visible Band Index for Remote Sensing Leaf Chlorophyll Content at the Canopy Scale.” International Journal of Applied Earth Observation and Geoinformation 21 (1): 103–12. https://doi.org/10.1016/j.jag.2012.07.020.
    out <- vegindex(spc, index)

  } else if(index == 'MTVI2_LAI'){
    # Haboudane, Driss, John R. Miller, Elizabeth Pattey, Pablo J. Zarco-Tejada, and Ian B. Strachan. 2004. “Hyperspectral Vegetation Indices and Novel Algorithms for Predicting Green LAI of Crop Canopies: Modeling and Validation in the Context of Precision Agriculture.” Remote Sensing of Environment 90 (3): 337–52. https://doi.org/10.1016/j.rse.2003.12.013.
    r800 <- get_reflectance(spc, wavelength = 800, weighted = weighted)
    r550 <- get_reflectance(spc, wavelength = 550, weighted = weighted)
    r670 <- get_reflectance(spc, wavelength = 670, weighted = weighted)
    num <- 1.5*(1.2*(r800-r550)-2.5*(r670-r550))
    deo <- sqrt((2*r800+1)^2-(6*r800-5*sqrt(r670))-0.5)
    out <- num/deo

  } else if(index == 'TTVI'){
    # LAI, Xing2019 10.3390/rs12010016
    out <- .calc_TTVI(spc)


    # rededge -----------------------------------------------------------------
  } else if(index %in% c('l0', 'lp', 'ls', 'R0', 'Rp', 'Rs')){
    out <- rededge(spc)[[index]]
  } else if(index  ==  'REP_lp'){
    out <- .calc_REP_sg(spc)$lp
  } else if(index == 'REP_Rp'){
    out <- .calc_REP_sg(spc)$Rp
  } else if(index == 'REP_Dp'){
    out <- .calc_REP_sg(spc)$Dp
  } else if(index == 'REP_multi'){
    out <-.calc_REP_multi(spc)
  } else if(index %in% c('REP_gaussian', 'mREIP')){
    # MILLER, J. R., E. W. HARE, and J. WU. 1990. “Quantitative Characterization of the Vegetation Red Edge Reflectance 1. An Inverted-Gaussian Reflectance Model.” International Journal of Remote Sensing 11 (10): 1755–73. https://doi.org/10.1080/01431169008955128.
    out <- .calc_mREIP(spc)
  } else if(index %in% c('REP_Li' ,'REP_LE')){
    out <- vegindex(spc, index)
  } else if(index == 'MTCI'){
    out <- vegindex(spc, index)

  } else if(index == 'IRECI'){
    # Frampton2013
    r783 <- get_reflectance(spc, 783)
    r665 <- get_reflectance(spc, 665)
    r705 <- get_reflectance(spc, 705)
    r740 <- get_reflectance(spc, 740)

    out <- (r783-r665)/(r705/r740)

  } else if(index == 'S2REP'){
    # Frampton, William James, Jadunandan Dash, Gary Watmough, and Edward James Milton. 2013. “Evaluating the Capabilities of Sentinel-2 for Quantitative Estimation of Biophysical Variables in Vegetation.” ISPRS Journal of Photogrammetry and Remote Sensing 82 (August): 83–92. https://doi.org/10.1016/j.isprsjprs.2013.04.007.
    r664 <- get_reflectance(spc, 664)
    r705 <- get_reflectance(spc, 705)
    r740 <- get_reflectance(spc, 740)
    r783 <- get_reflectance(spc, 783)

    out <-705 + 35*((r783-r664)/2-r705)/(r740-r705)

  } else if(index == 'IRECI'){
    # Frampton, William James, Jadunandan Dash, Gary Watmough, and Edward James Milton. 2013. “Evaluating the Capabilities of Sentinel-2 for Quantitative Estimation of Biophysical Variables in Vegetation.” ISPRS Journal of Photogrammetry and Remote Sensing 82 (August): 83–92. https://doi.org/10.1016/j.isprsjprs.2013.04.007.
    r664 <- get_reflectance(spc, 664)
    r705 <- get_reflectance(spc, 705)
    r740 <- get_reflectance(spc, 740)
    r783 <- get_reflectance(spc, 783)

    out <- (r783-r664)/(r705/r740)


    # LbyL --------------------------------------------------------------------
  } else if(index == 'NDVI_Delegido2013'){
    # Delegido, J., Verrelst, J., Meza, C. M. M., Rivera, J. P. P., Alonso, L., & Moreno, J. (2013). A red-edge spectral index for remote sensing estimation of green LAI over agroecosystems. European Journal of Agronomy, 46, 42–52. https://doi.org/10.1016/j.eja.2012.12.001
    out <- ssdxj_vegindex('NDVI_712_674', spc)


  }  else if (str_detect(index, "^NDVI_\\d+")) { # incase NRI calc by hsdar::nri
    # wls <- str_extract_all(index, '\\d+(.\\d+)?') %>% unlist() %>% as.numeric()
    wls <- str_split(index, "_")[[1]][c(2, 3)] %>% parse_double()
    if (length(wls) != 2) {
      msg <- glue::glue('error: {index} parsed wl is not with length 2({wls})')
      stop(msg)
    }
    b1 <- get_reflectance(spc, wavelength = wls[1], weighted = FALSE)
    b2 <- get_reflectance(spc, wavelength = wls[2], weighted = FALSE)
    out <- (b2 - b1) / (b2 + b1)


  } else if (str_detect(index, "^DVI_\\d+")) { # incase NRI calc by hsdar::nri
    # wls <- str_extract_all(index, '\\d+(.\\d+)?') %>% unlist() %>% as.numeric()
    wls <- str_split(index, "_")[[1]][c(2, 3)] %>% parse_double()
    if (length(wls) != 2){
      msg <- glue::glue('error: {index} parsed wl is not with length 2({wls})')
      stop(msg)
    }
    b1 <- get_reflectance(spc, wavelength = wls[1], weighted = FALSE)
    b2 <- get_reflectance(spc, wavelength = wls[2], weighted = FALSE)
    out <- b2 - b1
  } else if(str_detect(index, "^SR_\\d+")){
    wls <- str_split(index, "_")[[1]][c(2, 3)] %>% parse_double()
    if (length(wls) != 2){
      msg <- glue::glue('error: {index} parsed wl is not with length 2({wls})')
      stop(msg)
    }
    b1 <- get_reflectance(spc, wavelength = wls[1], weighted = FALSE)
    b2 <- get_reflectance(spc, wavelength = wls[2], weighted = FALSE)
    out <- b1/b2

    # others ------------------------------------------------------------------
  } else if(index == 'VNAI'){
    # Yue, J., Feng, H., Tian, Q., & Zhou, C. (2020). A robust spectral angle index for remotely assessing soybean canopy chlorophyll content in different growing stages. Plant Methods, 16(1), 104. https://doi.org/10.1186/s13007-020-00643-z
    out <- .calc_VNAI(spc)

  } else if(index == 'DANIR'){
    # Dutta, Dibyendu, Prabir Kumar Das, Kazi Arif Alam, P Safwan, Soubhik Paul, Manoj Kumar Nanda, and Vinay Kumar Dadhwal. 2016. “Delta Area at Near Infrared Region (DANIR)—A Novel Approach for Green Vegetation Fraction Estimation Using Field Hyperspectral Data.” IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing 9 (9): 3970–81. https://doi.org/10.1109/JSTARS.2016.2539359.
    out <- .calc_DANIR(spc)

  } else if(index == 'NAOC'){
    # Delegido, J., Alonso, L., González, G., & Moreno, J. (2010). Estimating chlorophyll content of crops from hyperspectral data using a normalized area over reflectance curve (NAOC). International Journal of Applied Earth Observation and Geoinformation, 12(3), 165–174. https://doi.org/10.1016/j.jag.2010.02.003

    out <- .calc_NAOC(spc)


    # finally -----------------------------------------------------------------

  } else {
    msg <- sprintf('Not defined index: %s!!!', index)
    stop(msg)
  }


  #check
  if (is.matrix(out)) {
    n <- dim(out)
    msg <- glue::glue('multicol of VI output.\n index:{index}\n dim:{n}')
    stop(n)
  }

  if (all(is.na(out))){
    msg <- glue::glue('NA output of index: {index}!!!')
    stop(msg)
  }

  return(out)
}


.calc_VNAI <- function(spc){
  .radians2degree <- function(x) x/pi*180.0

  wl_B <- 492.4
  wl_G <- 559.8
  wl_R <- 664.6
  wl_NIR <- 832.8
  ref_B <- get_reflectance(spc, wl_B, weighted = FALSE)
  ref_G <- get_reflectance(spc, wl_G, weighted = FALSE)
  ref_R <- get_reflectance(spc, wl_R, weighted = FALSE)
  ref_NIR <- get_reflectance(spc, wl_NIR)


  alpha <- pi - atan((ref_G-ref_B)/(wl_G-wl_B)*2500) - atan((ref_G-ref_R)/(wl_G-wl_R)*2500)
  beta <- pi - atan((ref_G-ref_B)/(wl_G-wl_B)*2500) + atan((ref_G-ref_NIR)/(wl_G-wl_NIR)*2500)

  .radians2degree(alpha + beta)

}


.calc_TTVI <- function(spc){
  r865 <- get_reflectance(spc, wavelength = 865, weighted = FALSE)
  r783 <- get_reflectance(spc, wavelength = 783, weighted = FALSE)
  r740 <- get_reflectance(spc, wavelength = 740, weighted = FALSE)
  out <- nno0.5*((783-740)*(r865-r740)-(865-740)*(r783-r740))

}

.calc_SF <- function(spc){
  r800 <- get_reflectance(spc, wavelength = 800, weighted = FALSE)
  r560<- get_reflectance(spc, wavelength = 560, weighted = FALSE)

  (r800^2-r560)/r800

}


.calc_DANIR <- function(spc){

  n <- nspectra(spc)
  wl <- wavelength(spc)
  ref <- spectra(spc)
  re_info <- rededge(spc)
  lp <- re_info[, 'lp']
  ls <- re_info[, 'ls']

  out <- rep(0, times = n)
  for(i in 1:n){
    # bands between lp and ls
    wl_flag <- wl>=lp[i] & wl <=ls[i]
    ref_in <- ref[i,wl_flag]
    R_ls <- ref_in[length(ref_in)]
    out[i] <- sum(R_ls-ref_in)
  }

  out

}


.calc_NAOC <- function(spc, a = 643, b = 795){
  wl <- wavelength(spc)
  ref <- spectra(spc)
  wl_in_flag <- wl >= a & wl <= b
  wl_in <- wl[wl_in_flag]
  n_wl_in <- length(wl_in)
  ref_in <- ref[,wl_in_flag]

  ref_int <- apply(ref_in, 1, sum)
  ref_max <- apply(ref_in, 1, max)
  area_sqare <- n_wl_in * ref_max

  NAOC <- 1 - ref_int/area_sqare

  return(NAOC)

}


.calc_mREIP <- function(spc){

# the hsdar::vegindex(spc, 'mREIP') DO NOT add sigma

  spec <- spectra(spc)
  wl <- wavelength(spc)
  wl_R0_flag <- wl >= 670 & wl <= 685
  wl_Rs_flag <- wl >= 780 & wl <= 795
  wl_Rl_flag <- wl >= 685 & wl <= 780

  if (wl[1] > 670)  return(rep.int(NA, nrow(y)))
  if (wl[length(wl)] < 795) return(rep.int(NA, nrow(y)))

  R0 <- matrix(data = apply(spec[, wl_R0_flag], 1, mean), ncol = 1)
  Rs <- matrix(data = apply(spec[, wl_Rs_flag], 1, mean), ncol = 1)
  Rl <- as.matrix(spec[, wl_Rl_flag])
  dat <- cbind(Rl, Rs, R0)

  mREIP_fun <- function(x, wl_Rl) {
    Rs_local <- x[length(x) - 1]
    R0_local <- x[length(x)]
    Rl_local <- x[1:(length(x) - 2)]

    # incase of na
    flag <- (Rs_local - Rl_local)  > 0
    Rl_local <- Rl_local[flag]
    wl_Rl <- wl_Rl[flag]

    Bl <- -1 * log(sqrt((Rs_local - Rl_local)/(Rs_local - R0_local)))
    coef <- summary(lm(Bl ~ wl_Rl))$coefficients
    a0 <- coef[1,1]
    a1 <- coef[2,1]
    lambda_0 <- -a0/a1
    sigma <- 1/(sqrt(2)*a1)
    return(lambda_0 + sigma)

  }

  apply(dat, 1, mREIP_fun, wl[wl_Rl_flag])
}



.calc_REP_sg <- function(spc){
  n_sgolay <- floor((25/mean(spc@fwhm))/2)*2+1
  if (n_sgolay < 5) n_sgolay <- 5
  D1 <- derivative.speclib(spc, m = 1, method = 'sgolay', n = n_sgolay)

  flag_region <- wavelength(spc) >= 680 & wavelength(spc) <= 780
  spec <- spectra(D1)
  spec[,!flag_region] <- -99999.9
  lp <- apply(spec, 1, which.max)
  Rp <- map2_dbl(as.data.frame(t(spectra(spc))), lp, ~.x[.y])
  Dp <- map2_dbl(as.data.frame(t(spectra(D1))), lp, ~.x[.y])
  lp <- wavelength(spc)[lp]

  names(Rp) <- NULL
  names(Dp) <- NULL

  return(list(Rp = Rp, lp = lp, Dp = Dp))
}



.calc_REP_multi <- function(spc){
  wl <- wavelength(spc)
  wl_sub_flag <- wl >= 680 & wl <= 750
  wl_sub <- wl[wl_sub_flag]
  n <- length(wl_sub)
  Raw_spec <- spectra(spc)[,wl_sub_flag]
  D1_spec <- spectra(derivative.speclib(spc, m = 1,
                                        method = 'sgolay', n = 21))[,wl_sub_flag]
  D2_spec <- spectra(derivative.speclib(spc, m = 2,
                                        method = 'sgolay', n = 21))[,wl_sub_flag]

  .REP_multi_apply <- function(i, n, wl_sub, Raw_spec, D1_spec, D2_spec){
    # incase
    i <- i[1]
    lp <- c()
    Rp <- c()
    Dp <- c()

    for(j in 3:(n-2)){
      if(D2_spec[i, j-2] > 0 &
         D2_spec[i, j-1] > 0 &
         D2_spec[i, j+1] <= 0 &
         D2_spec[i, j+2] <= 0){

        lp <- c(lp, wl_sub[j])
        Rp <- c(Rp, Raw_spec[i,j])
        Dp <- c(Dp, D1_spec[i,j])

        D2_spec[i, 1:j] <- -9999 # to exclude the j+1 channel
      }
    }

    if(length(lp) < 1){
      lp <- NA
      Rp <- NA
      Dp <- NA
    }

    # if(length(lp == 1)){
    #   lp <- c(lp, NA)
    #   Rp <- c(Rp, NA)
    #   Dp <- c(Dp, NA)
    # }

    if(length(lp) > 2) {
      print(i)
      print(lp)
      print('more than 2 lp are founded!!!')
    }

    c(lp1 = lp[1], lp2 = lp[2], Rp1 = Rp[1], Rp2 = Rp[2], Dp1 = Dp[1], Dp2 = Dp[2])

  }

  ids <- 1:nspectra(spc1nm)
  names(ids) <- ids

  map_df(ids, .REP_multi_apply, n, wl_sub, Raw_spec, D1_spec, D2_spec)
}


