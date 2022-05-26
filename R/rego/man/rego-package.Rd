\name{rego-package}
\alias{rego-package}
\alias{rego}
\docType{package}

\title{
Automatic Time Series Forecasting and Missing Value Imputation.
}

\description{
rego is a machine learning algorithm for predicting and imputing time series. It can automatically set all the parameters needed, thus in the minimal configuration it only requires the target variable and the regressors if present. It can address large problems with hundreds or thousands of dependent variables and problems in which the number of dependent variables is greater than the number of observations. Moreover it can be used not only with time series but also with any other real valued target variable. The algorithm implemented includes a Bayesian stochastic search methodology for model selection and a robust estimation based on bootstrapping. rego is fast because all the code is C++.
}

\details{
\tabular{ll}{
Package: \tab rego\cr
Type: \tab Package\cr
Version: \tab 1.5.1 \cr
Date: \tab 2022-05-26 \cr
License: GPL3 | see file LICENCE
}
 Package contains a function for automatic time series forecasting and missing values imputation. 
}

\author{
  Davide Altomare (\email{info@channelattribution.io}).
}

\references{
\href{https://channelattribution.io/docs/rego}{rego Documentation}
}

\keyword{time series forecasting}
\keyword{time series prediction}
\keyword{time series}
\keyword{sarimax}
\keyword{regression}
\keyword{linear regression}
\keyword{missing imputation}
\keyword{prediction intervals}

