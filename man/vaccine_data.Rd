% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vaccine_data.R
\docType{data}
\name{vaccine_data}
\alias{vaccine_data}
\title{Austrian vaccination rates dataset}
\format{
A  \code{data.frame}; each row represents one municipality;
\code{ID} Is an identifier for each municipality.
\code{agexy} represents the number of people that was of age xy at the beginning of 2021, for xy between 0 and 99.
\code{partly_vaccinated_p100} contains the percentage of the population that had received at least one vaccine against COVID-19,
\code{fully_vaccinated_p100} the percentage of the population that was considered fully vaccinated in October 2021.
}
\source{
\href{https://statcube.at}{STATcube - Statistische Datenbank von Statistik Austria}. Bevölkerung zu Jahresbeginn 2021, 2021. Accessed: 23.10.2021.

\href{https://www.data.gv.at/katalog/dataset/covid-19-schutzimpfungen-impfungen-in-gemeinden}{Bundesministerium für Soziales, Gesundheit, Pflege und Konsumentenschutz (BMSGPK)}. Covid-19 Schutzimpfungen – Impfungen in Gemeinden, 2021. Licensed under \href{https://creativecommons.org/licenses/by/4.0/}{CC-BY-4.0}. Accessed: 23.10.2021.
}
\usage{
data(vaccine_data)
}
\description{
Rates of vaccination against COVID-19 and age structure
in Austrian municipalities 2021.
}
\examples{
data(vaccine_data)
hist(vaccine_data$fully_vaccinated_p100, freq = FALSE,
xlab="\% fully vaccinated", main="Vaccination rates")
age.curves=vaccine_data[,3:102]
matplot(t(age.curves), type="l", xlab="Age", ylab="Number of People", main="Age structure")
}
\keyword{datasets}
