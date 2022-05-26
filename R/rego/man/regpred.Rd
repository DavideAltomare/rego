\name{regpred}


\alias{regpred}


\title{Automatic Time Series forecasting and Missing Value Imputation.}

\description{Automatic time series prediction and missing value imputation.}



\usage{
  regpred(Data, from_lag=1, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, 
          direction="->", loss_function="MSE", flg_const=TRUE, flg_diff=FALSE, model=NULL)
}


\arguments{
  \item{Data}{data.frame containing target variable at first column and regressors if present from second to last column.}
  \item{from_lag}{minimum time lag to be considered in the autoregressive moving average part of the algorithm.}
  \item{max_lag}{maximum time lag to be considered in the autoregressive moving average part of the algorithm. If "auto" then the algorithm will set a suitable value. Set to 0 or NULL if you want to remove the autoregressive moving average part as in case of non time series data.}
  \item{alpha}{significance level for the confidence interval produced around predictions. If 0.05 then the algorithm will calculate a 95\% confidence interval around predictions.}
  \item{nsim}{number of bootstrap replications used for producing confidence interval around predictions.}
  \item{flg_print}{if 1 some information during the evaluation will be printed.}
  \item{direction}{if "->" then only a forward prediction will be executed, if "<-" then only a backward prediction will be executed, if "<->" then both a forward than a backward prediction will be executed if possible. For imputing missing values is convenient to leave default "<->".}
  \item{loss_function}{if "MAE" then mean absolute error will be used as loss function in parameters estimation, if "MSE" then mean squared error will be set.}
  \item{flg_const}{if 1 then a constant is included into the model.}
  \item{flg_diff}{if 1 and no regressor is present then if the target variable exhibits a trend, it is one-step differentiated up to two times.}
  \item{model}{estimated models from a previous train to be used in new data prediction without retraining.}
}

\value{An object of \code{class} \code{list} with predictions and models.}

\author{
  Davide Altomare (\email{info@channelattribution.io}).
}


\references{
\href{https://channelattribution.io/docs/rego}{rego Documentation}
}

\examples{

\dontrun{

#example 1: seasonal time series

library(rego)

data(Data)

res=regpred(Data)

#print final prediction 
print(res$prediction)

#example 2: high dimensional problem

Data=read.csv(url("https://channelattribution.io/csv/Data_sim_1000.csv"),header=FALSE)

res=regpred(Data, max_lag=NULL)

#print final prediction 
print(res$prediction)


}

}

\keyword{time series forecasting}
\keyword{time series prediction}
\keyword{time series}
\keyword{sarimax}
\keyword{regression}
\keyword{linear regression}
\keyword{missing imputation}
\keyword{prediction intervals}
