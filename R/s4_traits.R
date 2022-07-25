#' S4 class for basic population traits
#'
#' @description Class for containing basic demographic traits.
#'
#' @name S4BasicTrait
#'
#' @slot value Input value for a particular trait of the population.
#' @slot range Range of possible values for a particular trait of the population.
#'
#' @seealso \code{\link{S4Population}}
setClass("BasicTrait",
         slots = c(
           value = 'ANY',
           range = 'ANY'
         ))

#' S4 class for population contact traits
#'
#' @description Class for containing contact trait details
#'
#' @name S4ContactTrait
#'
#' @slot value Input value for a particular trait of the population.
#' @slot range Range of possible values for a particular trait of the population.
#' @slot contact_matrix Contact probabilities between defined types of individuals in population.
#' @slot adj_matrix Sparse adjacency matrix for population.
#'
#' @seealso \code{\link{S4Population}}
#'
#' @importClassesFrom Matrix dgCMatrix
setClass("ContactTrait",
         contains = 'BasicTrait',
         slots = c(
           contact_matrix = 'matrix',
           adj_matrix = 'dgCMatrix'
         ))
