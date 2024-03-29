# rego: Automatic time series forecasting and missing value imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

library(rego)

r=c()

#USECASE 1

Data=read.csv("data/Data_air.csv",sep=";",header=FALSE)

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->")
print(res)

res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", model=res$model)
print(res1)


#USECASE 2

Data=read.csv("data/Data_sim_1000.csv",sep=",",header=FALSE)

res=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="->")
print(res)

res1=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="->", model=res$model)
print(res1)

