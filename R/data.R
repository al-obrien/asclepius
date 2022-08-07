#' COVID-19 case data for Alberta.
#'
#' A dataset containing case counts in Alberta, Canada from March 2020 to August 2022.
#'
#' @format A data frame with 879 rows and 8 columns.
#' \describe{
#'   \item{Date.reported.to.Alberta.Health}{Date information reported to Alberta Health}
#'   \item{Number.of.lab.tests}{Total number of lab tests performed on reporting day}
#'   \item{Cumulative.number.of.lab.tests}{Rolling total for lab tests for COVID-19}
#'   \item{Number.of.cases}{Number of new cases for COVID-19 reported to Alberta Health}
#'   \item{Cumulative.number.of.cases}{Rolling total of newly reported COVID-19 cases}
#'   \item{Cumulative.number.of.deaths}{Rolling total of deaths related to COVID-19 cases}
#'   \item{Number.of.deaths}{New deaths related to COVID-19 cases}
#'   \item{Percent.positivity}{Percent of lab tests positive for COVID-19}
#' }
#' @source \url{https://www.alberta.ca/stats/covid-19-alberta-statistics.htm#data-export}
"ab_covid"
