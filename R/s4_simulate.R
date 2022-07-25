#' S4 class for disease simulation
#' @description Class for simulating disease process using population and disease classes
#' @name S4Simulate
#' @slot timesteps Numeric value.
#' @slot Population S4Population class.
#' @slot Disease S4Disease class.
#'
setClass("Simulate",
         slots = c(
           timesteps = 'numeric',
           init_inf = 'numeric',
           states = 'list',
           totals = 'list',
           Population = 'Population',
           Disease = 'Disease'
         ))

#' Constructor for S4 Simulation class for disease propagation
#'
#' @description Creates the simulation class based upon provided input classes.
#'
#' @describeIn S4Simulate Constructor for disease simulation class
#' @param timesteps Numeric value.
#' @param init_inf Numeric value.
#' @param replications Integer value (NOT INTEGERATED YET...PERHAPS METHOD INSTEAD?).
#' @param Population S4Population class.
#' @param Disease S4Disease class.
#'
#' @export
run_simulation <- function(timesteps, init_inf, replications, Population, Disease) {

  timesteps <- as.numeric(timesteps)

  #first possible generic... (SIS)
  inf <- rep(0,Population@n)
  inf[sample(1:Population@n, init_inf, replace = FALSE)] <- 1
  stepdata <- vector(mode = "list", length = timesteps)
  stepdata[[1]] <- inf
  recov <- rep(0,Population@n)
  for(i in 2:timesteps) {
    inf <- calculate_prob_inf(Population@contact_structure@adj_matrix, inf, Disease@inf_prob)
    #inf <- ifelse(inf-recov<=0, 0, 1); # if SIR instead of SIS, drop the inf before new resolution calc
    recov <- calculate_resolution(inf, inf_period = Disease@infectious_period) #TODO dont include infected that just joined on this cycle
    inf <- inf-recov
    stepdata[[i]] <- inf # Track inf at each step to return
  }
  names(stepdata) <- paste0('t',as.character(1:timesteps))

  status_i <- as.data.frame(stepdata)
  status_s <- 1-status_i

  inf_t <- colSums(status_i) #TODO add via methods?
  sus_t <- colSums(status_s)

  # Second possible generic/method (SIR)

  # Third possible generic/method (SEIR)

  # Fourth.. possible generic/method (SEIRF)

  #TODO BETTER WAY TO ALLOW ANY COMBINATION!

  new("Simulate",
      timesteps = timesteps,
      init_inf = init_inf,
      states = list(susceptible = status_s, infected = status_i),
      totals = list(susceptible = sus_t, infected = inf_t),
      Population = Population,
      Disease = Disease)
}


setMethod("plot", signature =  "Simulate", function(x) {
  plot.new()
  plot.window(xlim = c(0, x@timesteps), ylim = c(1, x@Population@n))
  axis(1)
  axis(2)
  title(xlab = 'Time step', ylab = 'Count', main = 'SIS')
  lines(x = 1:x@timesteps, x@totals$infected, col = 'red')
  lines(x = 1:x@timesteps, x@totals$susceptible, col = 'green')
  box()
})


#TODO ADD METHOD TO EXTRACT TOTALS
#TODO ADD METHOD TO EXTRACT STATUS

# Add option for...how to sample, not exact sample?
# setGeneric("run_simulation_sis", function( Population, Disease, timesteps, init_inf) standardGeneric("run_simulation_sis"))
# setMethod("run_simulation_sis", signature("Population", "Disease"), function(Population, Disease, timesteps, init_inf) {
#
#   inf <- rep(0,Population@n)
#   inf[sample(1:Population@n, init_inf, replace = FALSE)] <- 1
#   stepdata <- vector(mode = "list", length = timesteps)
#   stepdata[[1]] <- inf
#
#   recov <- rep(0,Population@n) # May not be necessary if recovered == S in terms of reinfection ability
#
#   for(i in 2:timesteps) {
#     inf <- calculate_prob_inf(Population@contact_structure@adj_matrix, inf, Disease@inf_prob)
#     #inf <- ifelse(inf-recov<=0, 0, 1); # if SIR instead of SIS, drop the inf before new resolution calc
#     recov <- calculate_resolution(inf, inf_period = Disease@infectious_period) #TODO dont include infected that just joined on this cycle
#     inf <- inf-recov
#     stepdata[[i]] <- inf # Track inf at each step to return
#   }
#
#   names(stepdata) <- paste0('t',as.character(1:timesteps))
#   as.data.frame(stepdata)
#   #TODO Add specific object contains each slice details
# })

#TODO, for perhaps when adding events... may use signature for another METHOD
#TODO... option here is constructor of SIMULATION can have various methods inside that change based on what signature is provied (no events? then much simpler)
#TODO final decision on if want to returns sim object... or takes in sim object?
