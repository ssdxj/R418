#' wet Chl calc
#'
#' @param df data.frame (with A665, A649, A470, vol, area and mass column)
#'
#' @return original df with additional columns: LCC and LCD
#' @export
#'
calc_Chl <- function(df){

  # Dong, T., Shang, J., Chen, J. M., Liu, J., Qian, B., Ma, B., Morrison, M. J., Zhang, C., Liu, Y., Shi, Y., Pan, H., & Zhou, G. (2019). Assessment of Portable Chlorophyll Meters for Measuring Crop Leaf Chlorophyll Concentration. Remote Sensing, 11(22), 2706. https://doi.org/10.3390/rs11222706

  A664 <- df$A665
  A649 <- df$A649
  # A470 <- df$A470
  vol <- df$vol
  area <- df$area
  mass <- df$mass


  # Chla_solu <- (13.36*A664)-(5.19*A649)
  # Chlb_solu <- (27.43*A649) - (8.12*A664)
  # Chlxc_solu <- ((1000*A470) - (2.13*Chla_solu) - (97.63*Chlb_solu))/209

  # content in solution (mg/L)
  Chla_solu <- (13.95*A664) - (6.88*A649)
  Chlb_solu <- (24.96*A649) - (7.32*A664)
  # Chlxc_solu <- ((1000*A470) - (2.05*Chla_solu) - (114.8*Chlb_solu))/245


  # ug/cm^2
  Chla_area <- Chla_solu * vol / area * 1000
  Chlb_area <- Chlb_solu * vol / area * 1000
  # Chlxc_area <- Chlxc_solu * vol / area * 1000

  # mg/g
  Chla_mass <- Chla_solu * vol / mass
  Chlb_mass <- Chlb_solu * vol / mass
  # Chlxc_mass <- Chlxc_solu * vol / mass

  df$LCC <- Chla_area + Chlb_area
  df$LCD <- Chla_mass + Chlb_mass

  return(df)
}
