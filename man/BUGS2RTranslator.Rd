% Generated by roxygen2 (4.0.1): do not edit by hand
\name{BUGS2RTranslator}
\alias{BUGS2RTranslator}
\title{translates from BUGS to R}
\usage{
BUGS2RTranslator(bugsname = "norm", bugsparameters = list(c(1, 1)))
}
\arguments{
\item{bugsnames}{a list/vector of model names (as characters)}

\item{bugsparameters}{a list of model parameters
which correspond to the models in bugsnames
and are in the correct and standard order (e.g. c(mu, sigma)).}
}
\description{
BUGS2Rtranslator translates probabilty density functions
from BUGS to R.
}
\details{
Given a list of functions and optionally a list of parameters,
BUGS2Rtranslator will return the corresponding
functions in R and translate the parameters from their
BUGS formulation into a R compatible form.
}
\examples{
bugs2rTranslator("norm",list(c(10,100)))
}
\author{
Marco D. Visser
}
\concept{
MCMC diagnostics
}

