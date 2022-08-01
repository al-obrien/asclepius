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
#' @param transition_model Character or function to describe how transition occurs between disease states. Defaults include: SIS, SIR, and SEIR.
#'
#' @export
run_simulation <- function(timesteps, init_inf, replications = 1, Population, Disease, random_init = TRUE, transition_model) {

  timesteps <- as.numeric(timesteps)

  # Perform replications on provided model
  if(is.function(transition_model)) {
    state_list <- replicate(replications, transition_model(timesteps, init_inf, replications, Population, Disease, random_init)) #TODO provide better way to validate user provided function (class?)
    compartment_model <- 'Custom state model'
  } else {
    switch(transition_model,
      'SIS' = {
        state_list <- replicate(replications, sis(timesteps, init_inf, replications, Population, Disease, random_init))
        compartment_model <- 'si(s)'
        },
      'SIR' = {
        stopifnot(length(Disease@convalescence_period) > 0) # Must have this to be valid
        state_list <- replicate(replications, sir(timesteps, init_inf, replications, Population, Disease, random_init))
        compartment_model <- 'sir(s)'
      },
      'SEIR' = {
        state_list <- replicate(replications, seir(timesteps, init_inf, replications, Population, Disease, random_init))
        compartment_model <- 'seir(s)'
      })
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

#TODO allow to work with replicants (SE and all plots)
setMethod("plot", signature =  "Simulate", function(x, network = FALSE, replication = 1, timestep = 1, layout_f = igraph::layout_in_circle) {
  if(!network) {

    switch(x@compartment_model,
           'si(s)' = {
             plot.new()
             plot.window(xlim = c(0, x@timesteps), ylim = c(0, x@Population@n))
             axis(1)
             axis(2)
             if(ncol(x@states)>1) {
               title(xlab = 'Time step', ylab = paste0('Average count (', ncol(x@states), ' replications)'), main = 'SI(S)')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['infected_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'red')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['susceptible_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'green')
             } else {
               title(xlab = 'Time step', ylab = 'Count', main = 'SI(S)')
               lines(x = 1:x@timesteps,
                     x@states['infected_total',1][[1]],
                     col = 'red')
               lines(x = 1:x@timesteps,
                     x@states['susceptible_total',1][[1]],
                     col = 'green')
             }
             legend(x = "topright", legend = c("S", "I"), col = c('green', 'red'), lty = 1)
             box()
           },
           'sir(s)' = {
             plot.new()
             plot.window(xlim = c(0, x@timesteps), ylim = c(0, x@Population@n))
             axis(1)
             axis(2)
             if(ncol(x@states)>1) {
               title(xlab = 'Time step', ylab = paste0('Average count (', ncol(x@states), ' replications)'), main = 'SIR(S)')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['infected_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'red')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['susceptible_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'green')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['recovered_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'blue')
             } else {
               title(xlab = 'Time step', ylab = 'Count', main = 'SIR(S)')
               lines(x = 1:x@timesteps,
                     x@states['infected_total',1][[1]],
                     col = 'red')
               lines(x = 1:x@timesteps,
                     x@states['susceptible_total',1][[1]],
                     col = 'green')
               lines(x = 1:x@timesteps,
                     x@states['recovered_total',1][[1]],
                     col = 'blue')
             }
             legend(x = "topright", legend = c("S", "I", "R"), col = c('green', 'red', 'blue'), lty = 1)
             box()
           },
           'seir(s)' = {
             plot.new()
             plot.window(xlim = c(0, x@timesteps), ylim = c(0, x@Population@n))
             axis(1)
             axis(2)
             if(ncol(x@states)>1) {
               title(xlab = 'Time step', ylab = paste0('Average count (', ncol(x@states), ' replications)'), main = 'SEIR(S)')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['infected_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'red')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['susceptible_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'green')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['exposed_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'purple')
               lines(x = 1:x@timesteps,
                     Reduce(x@states['recovered_total',1:ncol(x@states)], f = `+`)/ncol(x@states),
                     col = 'blue')
             } else {
               title(xlab = 'Time step', ylab = 'Count', main = 'SEIR(S)')
               lines(x = 1:x@timesteps,
                     x@states['infected_total',1][[1]],
                     col = 'red')
               lines(x = 1:x@timesteps,
                     x@states['exposed_total',1][[1]],
                     col = 'purple')
               lines(x = 1:x@timesteps,
                     x@states['susceptible_total',1][[1]],
                     col = 'green')
               lines(x = 1:x@timesteps,
                     x@states['recovered_total',1][[1]],
                     col = 'blue')
             }
             legend(x = "topright", legend = c("S","E", "I", "R"), col = c('green', 'purple', 'red', 'blue'), lty = 1)
             box()
           })


  }
  if(network) {
    network_plot <- igraph::graph.adjacency(x@Population@contact_structure@adj_matrix, mode = 'undirected')
    igraph::V(network_plot)$color <- ifelse(x@states['infected', replication][[1]][[timestep]] == 1, 'red', 'lightgreen')
    igraph::E(network_plot)$color <- ifelse(igraph::tail_of(network_plot, igraph::E(network_plot))$color == 'red', 'tomato1', 'lightgrey')
    if (x@compartment_model %in% c('sir(s)', 'seir(s)')) {
      igraph::V(network_plot)$color <- ifelse(x@states['recovered', replication][[1]][[timestep]] == 1, 'lightblue', igraph::V(network_plot)$color)
      igraph::E(network_plot)$color <- ifelse(igraph::head_of(network_plot, igraph::E(network_plot))$color == 'lightblue' & igraph::E(network_plot)$color == 'tomato1',
                                              'steelblue',
                                              igraph::E(network_plot)$color)
      }
    coord <- layout_f(network_plot)
    plot(network_plot, layout = coord)
  }
})

#TODO BETTER WAY TO ALLOW ANY COMBINATION FOR STATE SELECTION, AND CUSTOM ONES? INSTEAD OF AUTOSLECT, HAVE IT PREDEFINED OR A CUSTOM FUNCTION  THAT THEN PERFORMS SOME BAIS CHECKS!
#TODO possible convert compartment transition models to a method class? and user provided ones generated on the fly?
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
