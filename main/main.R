# rego: Automatic time series forecasting and missing values imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

library(rego)
library(data.table)

#usecase 1

Data=fread("https://app.channelattribution.net/data/Data_air.csv",sep=";",header=FALSE)

res=regpred(Data, max_lag="auto", alpha=0.05,nsim=1000, flg_print=1, direction="<->")

x=res$prediction$final
plot(x$real,type="l",ylim=c(min(x,na.rm=TRUE)*0.95,max(x,na.rm=TRUE)*1.05))
lines(x$predicted,col="blue")
lines(x$lower_bound,col="red")
lines(x$upper_bound,col="red")

#use the saved model to get predictions
res1=regpred(Data, max_lag="auto", alpha=0.05,nsim=1000, flg_print=1, direction="<->", model=res[["model"]])

#usecase 2

Data=fread("https://app.channelattribution.net/data/Data_sim_1000.csv",sep=",",header=FALSE)

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->")

x=res$prediction$final
plot(x$real,type="l",ylim=c(min(x,na.rm=TRUE)*0.95,max(x,na.rm=TRUE)*1.05))
lines(x$predicted,col="blue")
lines(x$lower_bound,col="red")
lines(x$upper_bound,col="red")
