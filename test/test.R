# rego: Automatic time series forecasting and missing value imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

library(rego)

dir_princ="/home/Projects/GIT/personal/rego/repo/data/"

r=c()

#USECASE 1

Data=read.csv(paste0(dir_princ,"Data_air.csv"),sep=";",header=FALSE)


res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE", model=res$model)

r=c(r,all.equal(res$prediction,res1$prediction))

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE", model=res$model)

r=c(r,all.equal(res$prediction,res1$prediction))

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE", model=res$model)

r=c(r,all.equal(res$prediction$final,res1$prediction$final))

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE", model=res$model)

r=c(r,all.equal(res$prediction$final,res1$prediction$final))


#USECASE 2

Data=read.csv(paste0(dir_princ,"Data_sim_1000.csv"),sep=",",header=FALSE)

res=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE")
res1=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE", model=res$model)

r=c(r,all.equal(res$prediction,res1$prediction))

res=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE")
res1=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE",model=res$model)

r=c(r,all.equal(res$prediction,res1$prediction))

res=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE")
res1=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE",model=res$model)

r=c(r,all.equal(res$prediction$final,res1$prediction$final))

res=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE")
res1=regpred(Data, max_lag=NULL, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE",model=res$model)

r=c(r,all.equal(res$prediction$final,res1$prediction$final))


#USECASE 3

Data=read.csv(paste0(dir_princ,"Data_trading.csv"),sep=";")

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE", model=res$model)

r=c(r,all.equal(res$prediction,res1$prediction))

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE", model=res$model)

r=c(r,all.equal(res$prediction,res1$prediction))

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE", model=res$model)

r=c(r,all.equal(res$prediction$final,res1$prediction$final))

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE", model=res$model)

r=c(r,all.equal(res$prediction$final,res1$prediction$final))


print(r)