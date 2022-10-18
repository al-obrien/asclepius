# ------------------------------------- #
# Utility functions
# ------------------------------------- #

#' Calculate parameters of beta distribution
#'
#' Calculate alpha beta values of beta distribution from mean and variance.
#'
#' @param mu Numeric value for mean.
#' @param var Numeric value for variance
#'
calculate_alphabeta <- function(mu, var){

  if(mu < 0 | mu > 1) stop('mu must be between 0 and 1')
  if(!(var < mu*(1-mu))) stop('invalid var entered for provided mu')

  alpha <-  mu^2 * ( (1-mu)/var - 1/mu )
  beta <- mu * (1-mu) * ( ((1-mu)/var) - (1/mu))

  if(alpha <= 0 | beta <= 0) warning('alpha or beta values are <= 0 based upon these inputs...')

  return(list(alpha = alpha, beta = beta))

}

#' Calculate variance from alpha/beta parameters of beta distribution
#'
#' @param alpha Numeric value.
#' @param beta Numeric value.
#'
var_beta <- function(alpha, beta){
  (alpha * beta) / ( (alpha + beta)^2 * (alpha + beta + 1) )
}

#' Calculate mean from alpha/beta parameters of beta distribution
#'
#' @inheritParams var_beta
#'
mean_beta <- function(alpha, beta){
  alpha / (alpha + beta)
}

#' Base R merge with order maintained
#'
#' Can be used with R code to ensure adj matrix order is maintained.
#'
#' @param data1 Input data 1.
#' @param data2 Input data 2.
#' @param idrow Column that contains rows of IDs
#' @param sort Boolean to determine sorting.
#'
merge_keep_order <- function(data1, data2, idrow = NULL, sort = FALSE, ...) {
  data_comb <- merge(data1, data2, sort = sort, ...)
  row.names(data_comb) <- data_comb[[idrow]]
  data_comb <- data_comb[order(as.numeric(row.names(data_comb))),]
  data_comb
}


#TODO make it so if combos with similar category values can still be searched.... i.e. subset columns first?
#' Create validity matrix for contact patterns
#'
#' Used in creating contact patterns with restrictions.
#'
#' @param data Dataset containing the combinations of contact patterns.
#' @param rule_set List of rules.
#'
ruletest <- function(data, rule_set){

  # Check that sum of those that match is same as total entries that should match;
  sum_check <- colSums(apply(data, 1, function(x) {tabulate(match(x,rule_set), length(rule_set))} )) == length(rule_set)

  # Ensure that unique length matches the unique values provided
  length_check <- apply(data, 1, function(x) match(x,rule_set))
  length_check <- apply(length_check, 2, function(x) {length(unique(na.omit(x))) == length(unique(rule_set))} )

  comb_check <- length_check & sum_check
  return(comb_check)
}

#' Check if contact probabilities have violations to rule set
#'
#' Compare list of rules againt contact matrix probabilities to help assign violations to 0.
#'
#' @param var_list Columns containing values of contact patterns combinations (e.g. gender and age).
#' @param rule_list List of rules.
#' @param as_matrix Boolean.
#'
rule_check <- function(var_list, rule_list, as_matrix = TRUE){
  #var_list <- lapply(var_list, as.character) # do.call expand.grid makes all factrors...
  combos <- do.call(expand.grid, args = var_list)
  ncombo <- nrow(combos)
  combos <- do.call(expand.grid, args = c(var_list,var_list))

  temp <- Reduce(`+`, lapply(rule_list, ruletest, data = combos)) # ? .... + for all or just check for 1s?

  if(as_matrix) {
    return(matrix(temp, nrow = ncombo, ncol = ncombo))
  } else {
    return(temp)
  }
}


# TODO, combine w/ infectious period?
calculate_prob_inf <- function(contact_m, infected_v, beta, only_new = FALSE, only_all = FALSE) {
  if(beta > 1 || beta < 0) warning('Non-sense beta provided, outside of [0,1] bounds.')
  tau <- beta # beta * c(.1,.2,.3) ... Can add adjusted later... dep on contact type and other variables and time, offset by infec and by contact or events
  num_contacts <- contact_m %*% infected_v
  #TODO Add block vector to account for if state of exposed is dead or immune... maybe cut from network?
  #TODO if not blocked then calculate adjusted inf_prob by each contact...
  #TODO Add another function for exposure time instead of inf (so new 'infected' are swapped to E and then each start determines if goes to Inf)
  prob <- 1 - (1 - tau)^num_contacts # 1 - (Reduce(`*`, (1-tau))) where tau length is number of contacts  of various prob
  prob[infected_v==1] <- 0 # force existing infected to 0
  new_infected <- rbinom(length(prob), size = 1, as.vector(prob, mode = 'numeric'))
  all_inf <- infected_v + new_infected
  all_inf[all_inf > 1] <- 1 # Just in case...

  if(only_new) {
    return(new_infected)
    } else if (only_all) {
      return(all_inf)
      } else {
        return(list(new_inf = new_infected,
                    all_inf = all_inf))
      }
  }

calculate_resolution <- function(infected_v, inf_period) {
  if(inf_period <= 0) warning('Non-sense inf_period provided (<=0).')
  #TODO add gamma distirubiotn or custom params
  #TODO do so by having the inf_period as a vector assigned sampled from dis of choice...
  #TODO ...#comb_recov <- rate # Consider other factors to incr/decr rates?
  rate <- 1/inf_period # Rate at each time step
  new_resol <- rbinom(length(infected_v), infected_v, rate) # vector of inf b/c only trial is if infected already... and rate is constant
  new_resol
}

calculate_convalescence <- function(recov_v, conv_period) {

  rate <- 1/conv_period # Rate at each time step
  new_susc <- rbinom(length(recov_v), recov_v, rate)
  new_susc
}

#TODO combine with similar utility calcs based on rbinom method?
calculate_infectiousness <- function(expo_v, latent_p) {

  rate <- 1/latent_p # Rate at each time step
  new_expo <- rbinom(length(expo_v), expo_v, rate)
  new_expo
}

# Helper function to do various rescaling from one place
calculate_rescale <- function(x, type = 'maxmin', ...) {
  switch(type,
         'maxmin' = {scales::rescale(x, ...)}, # Between 0,1 or other arbitrary range
         'none' = {x},
         'simple' = {x / max(x, ...)}, # Similar to maxmin if min is 0
         'standardize' = {scale(x, ...)}, # mean 0, with sd of x
         'log' = {log(x, ...)},
         'robust' = { (x - median(x, ...)) / (quantile(x, probs = 0.75) - quantile(x, probs = 0.25))}, # Useful when many outliers
         'unitlength' = { x / sqrt(sum(x^2, ...))}, # Euclidean length (x / ||x||)
         'mean' = {(x - mean(x, ...)) / (max(x, ...) - min(x, ...))} # -1 to 1, mean 0
  )
}

