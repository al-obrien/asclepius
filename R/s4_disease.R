#' S4 class for disease
#'
#' @description Class for containing information on disease to simulation across a population.
#'
#' @name S4Disease
#'
#' @slot inf_prob Numeric value for probability of infection; used in calculation of \code{beta (new infected) = inf_prob * # of successful contacts}.
#' @slot infectious_period Numeric value for infectious period.
#' @slot latent_period Numeric value for length of latent period (useful for 'exposure period').
#' @slot incubation_period Numeric value for length of incubation period (NOT IN USE).
#' @slot convalescence_period Numeric value for length of convalescence period; defines minimum time until can be infected again, otherwise known as the immunity period (NOT IN USE).
#' @slot severity_coef Numeric value for severity factor of disease (NOT IN USE).
#'
#' @seealso \code{\link{S4Population}}
setClass("Disease",
         slots = c(
           inf_prob = 'numeric',
           infectious_period = 'numeric',
           latent_period = 'numeric',
           incubation_period = 'numeric',
           convalescence_period = 'numeric',
           severity_coef = 'numeric' # ?? Use as multiplier that builds towards hosp, icu, or death? How easy it is to progress into severe outcomes
         ))

#' Constructor for S4 Disease class.
#' @describeIn S4Disease Constructor for Disease class
#' @param inf_prob Numeric value for probability of infection.
#' @param inf_p Numeric value for infectious period.
#' @param latent_p Numeric value for latent period.
#' @param conv_p Numeric value for convalescence period.
#' @export
create_disease <- function(inf_prob, inf_p, latent_p = NULL,  conv_p = NULL) {
  inf_prob <- as.numeric(inf_prob)
  inf_p <- as.numeric(inf_p)
  conv_p <- as.numeric(conv_p)
  latent_p <- as.numeric(latent_p)

  new("Disease",
      inf_prob = inf_prob,
      infectious_period = inf_p,
      latent_period = latent_p,
      convalescence_period = conv_p)
}


# ------------------------------------- #
# Print method for Disease
# ------------------------------------- #

#' Show method for S4 Disease class.
#' @param object S4Disease class.
#' @describeIn S4Disease Show method for Disease class
setMethod("show", "Disease", function(object) {
  cat(is(object)[[1]], "\n",
      "  Infection probability:  ", object@inf_prob, "\n",
      "  Infectious period:      ", object@infectious_period, "\n",
      "  Latent period:          ", object@latent_period, "\n",
      "  Incubation period:      ", object@incubation_period, "\n",
      "  Convalescence period:   ", object@convalescence_period, "\n",
      sep = " "
  )
})

#TODO add base values on recovery rate, incubation period, infectious period, hosprate, icurate, asymp, fatality rate
#TODO each of the base disease values could be impacted by event range or population class...
#TODO change beta to tau? beta is inf prob and contact rate...
