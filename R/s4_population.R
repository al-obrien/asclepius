#' S4 class for population
#'
#' @description An S4 class to create a population for disease simulation.
#'
#' @name S4Population
#'
#' @slot n Integer, size of population.
#' @slot age_structure S4 class of 'BasicTrait' depicting age details of population.
#' @slot gender_structure S4 class of 'BasicTrait' depicting gender details of population.
#' @slot contact_structure S4 class of 'ContactTrait' depicting contact structure of population.
#' @include s4_traits.R
setClass("Population",
         slots = c(
           n = 'integer',
           age_structure = 'BasicTrait',
           gender_structure = 'BasicTrait',
           sex_orienation_structure = 'BasicTrait',
           contact_structure = 'ContactTrait'
         ),
         prototype = list(
           n = 0L,
           age_structure = new('BasicTrait', value = rep(NA_integer_, 1), range = NA_integer_),
           gender_structure = new('BasicTrait', value = rep(NA_character_, 1), range = NA_character_),
           sex_orienation_structure = new('BasicTrait', value = rep(NA_character_, 1), range = NA_character_),
           contact_structure = new('ContactTrait', value = rep(NA_integer_, 1), range = NA_integer_)
         ))


# Assign validity checks for S4 class Population
# TODO expand checks...
setValidity("Population", function(object) {
  if (object@n != length(object@age_structure@value)) {
    "@n and @age_structure must be same length"
  } else {
    TRUE
  }
})


#' Constructor for S4 Population class.
#' @describeIn S4Population Constructor for Population class
#' @param n Integer, size of population.
#' @export
create_population <- function(n) {
  n <- as.integer(n)
  new("Population", n = n, age_structure = new('BasicTrait', value = rep(NA_integer_, n), range = NA_integer_))
}


# ------------------------------------- #
# Print method for Population
# ------------------------------------- #

#' Show method for S4 Population class.
#' @param object S4Population class.
#' @describeIn S4Population Show method for Population class
setMethod("show", "Population", function(object) {
  cat(is(object)[[1]], "\n",
      "  Population size:  ", object@n, "\n",
      "  Age range:        ", range(object@age_structure@range)[[1]], ":", range(object@age_structure@range)[[2]], "\n",
      "  Gender values:    ", factor(object@gender_structure@range), "\n",
      "  Average contacts: ", mean(object@contact_structure@value, na.rm = TRUE), "\n",
      sep = " "
  )
})


# ------------------------------------- #
# Age generics/methods
# ------------------------------------- #

# S4Population Generic for retrieving age structure of Population class
setGeneric("get_age", function(pop_obj, value) standardGeneric("get_age"))

#' Retrieve age structure of Population (S4 Class).
#' @param pop_obj S4Population class.
#' @describeIn S4Population Method for retrieving age structure of Population class
#' @export
setMethod("get_age", "Population", function(pop_obj, value = NULL) {pop_obj@age_structure@value})


# S4Population Generic for setting age structure of Population class
# TODO Add option for...how to sample, not exact sample?
setGeneric("set_age", function(pop_obj, value = NULL, range) standardGeneric("set_age"))

#' Set age structure of Population (S4 Class).
#' @param value Vector of user-provided values for set operation (default: NULL, assigns a sample based on \code{range} parameter).
#' @param range Vector of min/max values for valid age.
#' @describeIn S4Population Method for setting age structure of Population class
#' @export
setMethod("set_age", "Population", function(pop_obj, value = NULL, range) {
  if(is.null(value)) value <- sample(range, pop_obj@n, replace = TRUE)
  pop_obj@age_structure <- new('BasicTrait', value = value, range = range)
  validObject(pop_obj)
  pop_obj
})

# ------------------------------------- #
# Gender generics/methods
#
# TODO add heterosexual... have sanity check against gender...
# ------------------------------------- #

# S4Population Generic for retrieving gender structure of Population class
setGeneric("get_gender", function(pop_obj, value) standardGeneric("get_gender"))

#' Retrieve gender structure of Population (S4 Class).
#' @describeIn S4Population Method for retrieving gender structure of Population class
#' @export
setMethod("get_gender", "Population", function(pop_obj, value = NULL) {pop_obj@gender_structure@value})


# S4Population Generic for setting gender structure of Population class
setGeneric("set_gender", function(pop_obj, value = NULL, range) standardGeneric("set_gender"))

#' Set gender structure of Population (S4 Class).
#' @describeIn S4Population Method for setting gender structure of Population class
#' @export
setMethod("set_gender", "Population", function(pop_obj, value = NULL, range) {
  if(is.null(value)) value <- sample(range, pop_obj@n, replace = TRUE)
  value <- as.character(value)
  pop_obj@gender_structure <- new('BasicTrait', value = value, range = range)
  validObject(pop_obj)
  pop_obj
})


# ------------------------------------- #
# Contact structure generics/methods

# Assign the number of contacts, determine distance
#TODO add complete homogeneous (i.e. ignore contact), or create at random from igraph options
#TODO ASSIGN RANDOM GRAPH OF THEIR CHOICE OF SET ALGORITHM
# ------------------------------------- #

# S4Population Generic for setting contact structure of Population class
setGeneric("set_contacts", function(pop_obj, value = NULL, range, vars = c(), rule_list = NULL, mu, variance, random_fill = TRUE, progress = TRUE) standardGeneric("set_contacts"))

#' Set contact structure of Population (S4 Class).
#' @param vars Vector of variable names of interest to create contact patterns (e.g. gender and age).
#' @param rule_list List of values defining valid contact patterns (default: NULL, all combinations of \code{vars} considered probable).
#' @param mu Numeric value for \eqn{mu} for beta distribution (\eqn{mu = a / (a + b)}); determines heterogeneity of contact patterns.
#' @param variance Numeric value for beta distribution; used with \eqn{mu} to determine beta distribution shape parameters for contact pattern sampling.
#' @param random_fill Boolean value to determine if population is filled in random ordering (default: \code{TRUE}); otherwise will go in provided row order.
#' @param progress Boolean value to determine whether or not a progress bar is shown for creating the network; helfpul if you need networks for population of 10,000 and over.
#' @describeIn S4Population Method for setting contact structure of Population class
#' @export
setMethod("set_contacts", "Population", function(pop_obj,
                                                 value = NULL,
                                                 range,
                                                 vars = c(),
                                                 rule_list = NULL,
                                                 mu, variance,
                                                 random_fill = TRUE,
                                                 progress = TRUE) {

  # Create blank matrix
  adj_mat_sp <- Matrix::Matrix(data = 0, nrow = pop_obj@n, ncol = pop_obj@n, sparse = TRUE)

  # Fill order by contact info
  if(random_fill) fillO <- sample(seq_along(1:pop_obj@n), size = pop_obj@n, replace = F) else fillO <- seq_along(1:pop_obj@n)

  # Assign contacts based on provided value/range or just sample from a range
  if(length(range) != 2) stop('range paramter must be length of 2')
  if(is.null(value)) value <- sample(seq(range[1], range[2]), pop_obj@n, replace = TRUE)
  pop_obj@contact_structure@value <- as.integer(value)
  pop_obj@contact_structure@range <- range

  # Create combinations of variables of interest for contact
  if(length(vars > 0)) {

    # Vars provided must be one of the available traits in the class
    stopifnot(!names(vars) %in% slotNames(pop_obj))

    var_list <- lapply(vars, \(x) slot(pop_obj, x)@range)
    dat_list <- lapply(vars, \(x) slot(pop_obj, x)@value)

    stopifnot("At least one of provided `vars` is all NA." = vapply(dat_list, \(x) !all(is.na(x)), logical(1))) # If all are NA for a var then dont run...

    names(var_list) <- vars
    names(dat_list) <- vars
    dat_list$rowid <- seq_along(1:pop_obj@n)
    dat_list$contacts <- pop_obj@contact_structure@value
    dat_list <- as.data.frame(dat_list, stringsAsFactors = FALSE)

    combined_attr <- do.call(expand.grid, args = var_list)
    combined_attr$ID <- 1:nrow(combined_attr)
    num_combs <- nrow(combined_attr)

    dis_matrix <- as.matrix(cluster::daisy(metric = "gower", x = combined_attr))
    diag(dis_matrix) <- .Machine$double.eps

    # Sort distance for each combination and draw 0,1 from beta for probabilities based on that distance (0 being similar, 1 being far)
    #TODO allow mixed? .5*dbeta(seq(0,1,0.01) ,1,3) + .5*dbeta(seq(0,1,0.01) ,5,1.75)... provide custom function...
    tmp_prob <- apply(dis_matrix[,,drop=FALSE], 2,
                      FUN = function(X) stats::dbeta(sort(X),
                                                     shape1 = calculate_alphabeta(mu, variance)[[1]],
                                                     shape2 = calculate_alphabeta(mu, variance)[[2]]) / length(X),
                      simplify = FALSE)

    # Convert to matrix of prob with ordering restored
    cnct_matrix <- matrix(unlist(lapply(tmp_prob,
                                        function(x) x[order(as.numeric(names(x)))])),
                          nrow = num_combs,
                          ncol = num_combs,
                          dimnames = list(as.character(1:num_combs),
                                          as.character(1:num_combs)))

    if(!is.null(rule_list)) {
      #TODO need to improve to allow different patterns (e.g. cant figure out if two combos both use 0,1 patterns!)
      # Flip those that violate rule to 0 probability
      cnct_matrix <- cnct_matrix * rule_check(var_list, rule_list, as_matrix = TRUE)
    }

    cmb_dat_wAttr <- merge_keep_order(dat_list, combined_attr, 'rowid')



  # ------------------- #
  # Fill adj matrix
  # ------------------- #
  adj_mat_sp <- init_adj_matrix_cpp(adj_mat_sp, cnct_matrix, cmb_dat_wAttr$ID, cmb_dat_wAttr$contacts, fillO, display_progress = progress)


  } else {

    # ------------------- #
    # Fill adj matrix
    # ------------------- #

    # Matrix of 1 cell for all same contact pattern
    cnct_matrix <- matrix(mu, 1, dimnames = list('1', '1'))

    adj_mat_sp <- init_adj_matrix_cpp(adj_mat_sp,
                                      cnct_matrix,
                                      rep(1, pop_obj@n), # All have same contact pattern
                                      pop_obj@contact_structure@value,
                                      fillO, display_progress = progress)

    # If matrix has no contacts, throw a warning
    if(sum(adj_mat_sp, na.rm = TRUE) < 1) warning('Adjacency matrix is all 0 (i.e. has no contacts filled). Double check that this is intentional!');
  }

  # Populate slot for final outputs
  pop_obj@contact_structure@contact_matrix <- cnct_matrix
    pop_obj@contact_structure@adj_matrix <- adj_mat_sp
    validObject(pop_obj)
  pop_obj

})


