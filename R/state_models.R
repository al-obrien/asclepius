# SIS MODEL
# (SIS), will be SI if inf period is INF
sis <- function(timesteps, init_inf, replications, Population, Disease, random_init) {
  inf <- rep(0,Population@n)
  if(random_init) {
    inf[sample(1:Population@n, init_inf, replace = FALSE)] <- 1
  } else {inf[1:init_inf] <- 1}
  stepdata <- vector(mode = "list", length = timesteps)
  stepdata[[1]] <- inf
  recov <- rep(0,Population@n)
  for(i in 2:timesteps) {
    inf <- calculate_prob_inf(Population@contact_structure@adj_matrix, inf, Disease@inf_prob)
    recov <- calculate_resolution(inf, inf_period = Disease@infectious_period) #TODO dont include infected that just joined on this cycle
    inf <- inf-recov
    stepdata[[i]] <- inf # Track inf at each step to return
  }
  names(stepdata) <- paste0('t',as.character(1:timesteps))

  status_i <- as.data.frame(stepdata)
  status_s <- 1-status_i

  inf_t <- colSums(status_i) #TODO add via methods?
  sus_t <- colSums(status_s)

  list(susceptible = status_s, infected = status_i, susceptible_total = sus_t, infected_total = inf_t)
}

# SIR MODEL
# Second possible generic/method (SIR and SIRS)
sir <- function(timesteps, init_inf, replications, Population, Disease, random_init) {
  inf <- rep(0,Population@n)
  if(random_init) {
    inf[sample(1:Population@n, init_inf, replace = FALSE)] <- 1
  } else {inf[1:init_inf] <- 1}
  stepdata <- vector(mode = "list", length = timesteps)
  stepdata[[1]] <- inf

  recov <- rep(0,Population@n)
  stepdata_r <- vector(mode = "list", length = timesteps)
  stepdata_r[[1]] <- recov

  for(i in 2:timesteps) {
    inf <- calculate_prob_inf(Population@contact_structure@adj_matrix, inf, Disease@inf_prob)
    inf <- ifelse(inf-recov<=0, 0, 1) # If was in recovered state at beginning revert infection back to 0

    conv <- calculate_convalescence(recov, Disease@convalescence_period) # Determine who will be susc again..
    recov <- ifelse(recov-conv<=0, 0, 1) # Flip recov to 0 if made past conval period

    new_recov <- calculate_resolution(inf, inf_period = Disease@infectious_period) #TODO dont include infected that just joined on this cycle
    inf <- inf-new_recov # Remove recovered for end of this cycle from inf

    recov <- recov + new_recov

    stepdata[[i]] <- inf # Track inf at each step to return
    stepdata_r[[i]] <- recov
  }
  names(stepdata) <- paste0('t',as.character(1:timesteps))

  status_i <- as.data.frame(stepdata)
  status_r <- as.data.frame(stepdata_r)
  status_s <- 1-status_i-status_r

  inf_t <- colSums(status_i) #TODO add via methods?
  sus_t <- colSums(status_s)
  rec_t <- colSums(status_r)

  list(susceptible = status_s, infected = status_i,recovered = status_r,
       susceptible_total = sus_t, infected_total = inf_t, recovered_total = rec_t)
}

# Third possible generic/method (SEIR)

# Fourth.. possible generic/method (SEIRF)
