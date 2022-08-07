ab_covid <- read.csv('data-raw/covid-19-alberta-statistics-summary-data.csv')

usethis::use_data(ab_covid, overwrite = TRUE)
