
# -------------------------------- #
# SIS MODEL, will be SI if inf period is INF
# -------------------------------- #
sis <- function(timesteps, init_inf, replications, Population, Disease, assign_init) {

  # Initialize infected
  inf <- rep(0,Population@n)
  if(assign_init == 'random') {
    inf[sample(1:Population@n, init_inf, replace = FALSE)] <- 1
  } else if (assign_init == 'index') { # Not random but want to assign specific Infected
    inf[init_inf] <- 1
  } else { # If not a vector, assign from top down by number
    inf[1:init_inf] <- 1
  }
  stepdata <- vector(mode = "list", length = timesteps)
  stepdata[[1]] <- inf

  new_i <- vector(mode = "numeric", length = timesteps)
  new_i[[1]] <- sum(inf) # Sum b/c could be using index method

  # Initialize recovered (instantly susceptible again)
  recov <- rep(0,Population@n)

  for(i in 2:timesteps) {
    inf_l <- calculate_prob_inf(Population@contact_structure@adj_matrix, inf, Disease@inf_prob, )
    inf <- inf_l$all_inf
    new_inf <- inf_l$new_inf

    recov <- calculate_resolution(inf, inf_period = Disease@infectious_period) #TODO dont include infected that just joined on this cycle
    inf <- inf-recov
    stepdata[[i]] <- inf # Track inf at each step to return
    new_i[[i]] <- sum(new_inf)
  }

  names(stepdata) <- paste0('t',as.character(1:timesteps))
  names(new_i) <- paste0('t',as.character(1:timesteps))

  status_i <- as.data.frame(stepdata)
  status_s <- 1-status_i

  inf_t <- colSums(status_i) #TODO add via methods?
  sus_t <- colSums(status_s)

  list(susceptible = status_s,
       infected = status_i,
       infected_new = new_i,
       susceptible_total = sus_t,
       infected_total = inf_t)
}


# -------------------------------- #
# SIR MODEL, generic/method (SIR and SIRS)
# -------------------------------- #

sir <- function(timesteps, init_inf, replications, Population, Disease, assign_init) {

  # Initialize infected
  inf <- rep(0,Population@n)
  if(assign_init == 'random') {
    inf[sample(1:Population@n, init_inf, replace = FALSE)] <- 1
  } else if (assign_init == 'index') { # Not random but want to assign specific Infected
    inf[init_inf] <- 1
  } else { # If not a vector, assign from top down by number
    inf[1:init_inf] <- 1
  }
  stepdata <- vector(mode = "list", length = timesteps)
  stepdata[[1]] <- inf
  new_i <- vector(mode = "numeric", length = timesteps)
  new_i[[1]] <- sum(inf) # Sum b/c could be using index method

  # Initialize recovered (different from new_recov)
  recov <- rep(0,Population@n)
  stepdata_r <- vector(mode = "list", length = timesteps)
  stepdata_r[[1]] <- recov

  for(i in 2:timesteps) {
    inf_l <- calculate_prob_inf(Population@contact_structure@adj_matrix, inf, Disease@inf_prob)
    inf <- inf_l$all_inf
    new_inf <- inf_l$new_inf
    inf[inf-recov<=0] <- 0 # If was in recovered state at beginning revert infection back to 0

    conv <- calculate_convalescence(recov, Disease@convalescence_period) # Determine who will be susc again..
    recov[recov-conv<=0] <- 0 # Flip recov to 0 if made past conval period

    new_recov <- calculate_resolution(inf, inf_period = Disease@infectious_period) #TODO dont include infected that just joined on this cycle
    inf <- inf-new_recov # Remove recovered for end of this cycle from inf

    recov <- recov + new_recov

    stepdata[[i]] <- inf # Track inf at each step to return
    stepdata_r[[i]] <- recov
    new_i[[i]] <- sum(new_inf)
  }

  names(stepdata) <- paste0('t',as.character(1:timesteps))
  names(stepdata_r) <- paste0('t',as.character(1:timesteps))
  names(new_i) <- paste0('t',as.character(1:timesteps))

  status_i <- as.data.frame(stepdata)
  status_r <- as.data.frame(stepdata_r)
  status_s <- 1-status_i-status_r

  inf_t <- colSums(status_i) #TODO add via methods?
  sus_t <- colSums(status_s)
  rec_t <- colSums(status_r)

  list(susceptible = status_s,
       infected = status_i,
       recovered = status_r,
       infected_new = new_i,
       susceptible_total = sus_t,
       infected_total = inf_t,
       recovered_total = rec_t)
}


# -------------------------------- #
# SEIR
# -------------------------------- #

#NOTE YET FUNCTIONING WITH EXPOSED?

seir <- function(timesteps, init_inf, replications, Population, Disease, assign_init) {

  # Create empty states
  recov <- expo <- inf <- rep(0,Population@n)
  stepdata_r <- stepdata_e <- stepdata <- vector(mode = "list", length = timesteps)
  new_i <- vector(mode = "numeric", length = timesteps)

  # Add initial states infected to inf vector
  if(assign_init == 'random') {
    inf[sample(1:Population@n, init_inf, replace = FALSE)] <- 1
  } else if (assign_init == 'index') { # Not random but want to assign specific Infected
    inf[init_inf] <- 1
  } else { # If not a vector, assign from top down by number
    inf[1:init_inf] <- 1
  }
  stepdata[[1]] <- inf
  stepdata_r[[1]] <- recov
  stepdata_e[[1]] <- expo
  new_i[[1]] <- sum(inf) # Sum b/c could be using index method

  # Update states
  for(i in 2:timesteps) {

    # Newly exposed, S->E
    new_expo <- calculate_prob_inf(Population@contact_structure@adj_matrix, inf, Disease@inf_prob, only_new = TRUE)
    new_expo[new_expo-recov <= 0] <- 0 # If was in recovered state at beginning revert infection back to 0
    new_expo[new_expo-expo <= 0] <- 0 # Necessary to ensure exposures that transition this cyle arent recycled in expo

    # Newly susceptible (R->S)
    conv <- calculate_convalescence(recov, Disease@convalescence_period) # Determine who will be susc again..
    recov[recov-conv<=0] <- 0 # Flip recov to 0 if made past conval period

    # Newly recovered (I->R)
    new_recov <- calculate_resolution(inf, inf_period = Disease@infectious_period) #TODO dont include infected that just joined on this cycle
    inf <- inf-new_recov # Remove recovered for end of this cycle from inf

    # Newly infected (E->I)
    new_inf <- calculate_infectiousness(expo, Disease@latent_period)
    expo[expo-new_inf<=0] <- 0 # Flip to 0 due to transition out of expo and into infected

    # Final state updates
    inf <- inf + new_inf; inf[inf>1] <- 1; # Ensure state is still 1 for overlapping infected?
    recov <- recov + new_recov
    expo <- expo + new_expo; expo[expo>1] <- 1;

    stepdata[[i]] <- inf # Track inf at each step to return
    stepdata_r[[i]] <- recov
    stepdata_e[[i]] <- expo
    new_i[[i]] <- sum(new_inf)
  }
  names(stepdata_e) <- names(stepdata_r) <- names(stepdata) <- paste0('t',as.character(1:timesteps))
  names(new_i) <- paste0('t',as.character(1:timesteps))

  status_i <- as.data.frame(stepdata)
  status_r <- as.data.frame(stepdata_r)
  status_e <- as.data.frame(stepdata_e)
  status_s <- 1-status_i-status_r-status_e

  inf_t <- colSums(status_i) #TODO add via methods?
  exp_t <- colSums(status_e)
  sus_t <- colSums(status_s)
  rec_t <- colSums(status_r)

  list(susceptible = status_s,
       exposed = status_e,
       infected = status_i,
       recovered = status_r,
       infected_new = new_i,
       susceptible_total = sus_t,
       exposed_total = exp_t,
       infected_total = inf_t,
       recovered_total = rec_t)
}

#TODO add fatal status...exit contact network
# Fourth.. possible generic/method (SEIRF)
seirf <- function(timesteps, init_inf, replications, Population, Disease, assign_init) { }
