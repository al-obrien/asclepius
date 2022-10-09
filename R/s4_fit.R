#TODO SEPARATE FROM SIMULATION CLASS, BUT WILL CREATE WITHIN FOR LOSS FUNCTION?
# PROB CAN SEPARATE SINCE NEED TO REPEAT SIM FOR LOSS FUNCTION?

#' S4 class for fitting of disease simulation
#' @description Class for fitting the simulations of disease processes
#' @name S4Fit
#' @slot timesteps Numeric value.
#'
setClass("DiseaseFit",
         slots = c(
           data = 'data.frame', # Data to fit against
           loss_function = 'function' # Which target to minimize difference?
         ))



#' Loss function for parameter estimation
#'
#' Helper function to take a set of parameters, constants, and incorporate them dynamically into
#' disease simulation based upon the transition model. Estimation is based upon a particular target state and
#' the data vector to match the trajectory to. This can be useful in conjunction with optimization algorithms suitable
#' for bounded integers such as select routines in `optim()`.
#'
#' Constant values provided will be filled into the params vector. This will determine which values one wants to optimize.
#'
#' @param params Parameters to optimize over.
#' @param constant Values that should remain constant during loss calculation and following optimization.
#' @param data Vector of values to calculate loss based on output of the disease simulation of the `target` state (e.g. infected_new).
#' @examples
#' \dontrun{
#' # Create a sample population
#' sim_pop <- create_population(1000)
#' sim_pop <- set_age(sim_pop, range = c(1:100))
#' sim_pop <- set_gender(sim_pop, range = c('M','F'))
#' sim_pop <- set_contacts(pop_obj = sim_pop, range = c(1,10), vars = c('age_structure', 'gender_structure'), mu = 0.33, variance = .001)
#'
#' data(ab_covid)
#' loss_func(c(3, 0.25, 10, 3, 15),
#'           Population = sim_pop,
#'           data = ab_covid$Number.of.cases,
#'           times = 50,
#'           transition_model = 'SEIR')
#'
#' # Terrible grid search
#' samplegrid <- expand.grid(seq(1:5), seq(.05, 0.5,.05), 3:15, 2:5, 5:15)
#' samplegrid <- samplegrid[sample(1:nrow(samplegrid), 500, replace = FALSE),]
#' model_fits <- apply(samplegrid, MARGIN = 1, simplify = T,
#'                     FUN = function(x) {loss_func(x,
#'                                                  Population = sim_pop,
#'                                                  data = ab_covid$Number.of.cases,
#'                                                  times = 50,
#'                                                  transition_model = 'SEIR')})
#'
#' # Run terrible grid search in parallel
#' library(parallel)
#' core_n <- detectCores()
#' cl <- makeCluster(core_n-1)
#' clusterExport(cl, 'samplegrid')
#' clusterExport(cl, 'sim_pop')
#' clusterEvalQ(cl, {library(asclepius)})
#' model_fits <- parApply(cl, samplegrid, MARGIN = 1,
#'                        FUN = function(x) {loss_func(x,
#'                                                     Population = sim_pop,
#'                                                     data = ab_covid$Number.of.cases,
#'                                                     times = 50, transition_model = 'SEIR')})
#' #samplegrid[which(rownames(samplegrid) == names(model_fits[which.min(model_fits )])),]
#'
#' }
#' @export
loss_func <- function(params,
                      constants = NULL,
                      data,
                      Disease = NULL,
                      Population = NULL,
                      times,
                      transition_model,
                      loss_function = 'rmse',
                      target = 'infected_new',
                      rescale = 'maxmin',
                      ...) {

  # Which are null (to be selective of which params I care about)
  #TODO checks on valid transiton based on params + constants
  #TODO outline list where param1 = init, param2 = inf_prob...
  null_list <- which(is.null(params))
  if(!is.null(constants) && length(null_list) == 0) warning('Provide constants for those not given in params')
  stopifnot(length(constants)==length(null_list))
  params[null_list] <- constants

  # Dis
  if(is.null(Disease)){
    Disease <- create_disease(inf_prob = params[2],
                              inf_p = params[3],
                              latent_p = params[4],
                              conv_p = params[5])
  }

  # Pop
  # if(is.null(Population)) {
  #   Population <- create_population()
  # }


  # Sim with parameters
  out <- run_simulation(Population = Population,
                        Disease = Disease,
                        init_inf = params[1],
                        replications = 5,
                        timesteps = times,
                        transition_model = transition_model)

  # Combine replications... and make DF, standardize
  states <- out@states[target,,drop = F]
  y_hat <- apply(states, 1, function(x) Reduce(`+`, x)/length(x))

  # Normalize 0-1
  y_hat_s <- calculate_rescale(y_hat, rescale)

  # Assign comparisons
  y <- calculate_rescale(data[1:times], rescale)

  # Cost to minimize
  #TODO allow more params to be passed into here? Custom function?
  # Other methods to add....
  # wrss = {loss <- sum(w * (y_hat-y)^2)}
  # mle_p = {loss <- sum(dpois(y, RATEATTIMEX, log = TRUE))})

  switch(loss_function,
         rmse = {loss <- sqrt(mean((y_hat-y)^2))},
         rss = {loss <- sum((y_hat-y)^2)},
         mse = {loss <- mean((y_hat-y)^2)},
         mae = {loss <- mean(abs(y_hat-y))})



  loss
}



# NB!!!!! NEED TO DETERMINE UTILITY OF THIS OVER JUST OPTIM OR OTHER GA TYPE METHODS????

#' Constructor for S4 Fitting of simulations to data
#'
#' @description Creates the fitting class for disease simulations.
#'
#' @describeIn S4DiseaseFit Constructor for fitting disease simulations
#' @param data Provided data for fitting.
#' @param Population S4Class for Population.
#' @param Disease S4Class for disease.
#' @param loss_function Function or name of how loss is calculated for fitting.
#' @param optim_method
#' @noRd
run_fitting <- function(data, Population, Disease, loss_function, optim_method) {

  # #TODO COMPLETE A LOSS FUNCTION
  # loss_f <- function(loss_function = loss_function) {
  #   #
  #   y_hat <- 1 # replace with simulated values?
  #   y <- data$count # Provided timeseries count...
  #
  #   # Cost to minimize
  #   switch(loss_function,
  #          'rmse' = {loss <- sqrt(mean((y_hat-y)^2))},
  #          'mse' = {loss <- mean((y_hat-y)^2)},
  #          'mae' = {loss <- mean(abs(y_hat-y))})
  #
  #   return(loss)
  # }

  result <- optim(par = init_par,
                  fn = loss_f,
                  method = optim_method,
                  times = times,
                  data = sim_data,
                  ...)

  new("DiseaseFit")
}

















