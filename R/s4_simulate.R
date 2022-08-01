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
           states = 'matrix',
           compartment_model = 'character',
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
#' @param random_init Boolean to control if random set of infected in population selected for each simulation.
#'
#' @export
run_simulation <- function(timesteps, init_inf, replications = 1, Population, Disease, random_init = TRUE) {

  timesteps <- as.numeric(timesteps)

  # Perform replications
  if(length(Disease@convalescence_period) > 0) {
    state_list <- replicate(replications, sir(timesteps, init_inf, replications, Population, Disease, random_init))
    compartment_model <- 'sir(s)'
  } else {
    state_list <- replicate(replications, sis(timesteps, init_inf, replications, Population, Disease, random_init))
    compartment_model <- 'si(s)'
  }

  colnames(state_list) <- paste0('replicant', 1:replications)
  rownames(state_list) <- names(state_list[,1])

  new("Simulate",
      timesteps = timesteps,
      init_inf = init_inf,
      states = state_list,
      compartment_model = compartment_model,
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


#TODO BETTER WAY TO ALLOW ANY COMBINATION FOR STATE SELECTION, AND CUSTOM ONES?
#TODO ADD METHOD TO EXTRACT TOTALS
#TODO ADD METHOD TO EXTRACT STATUS
#TODO network at each step!
#TODO standardize the y axis, so its 0-1 instead of counts...

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
