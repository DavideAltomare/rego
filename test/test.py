# rego: Automatic time series forecasting and missing value imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 


import pandas as pd
from rego import *

r=[]

#USECASE 1

Data=pd.read_csv("data/Data_air.csv",sep=";",header=None)
Data.columns=["target"]

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE", model=res['model'])

r=r + [res['prediction'].equals(res1['prediction'])]

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE", model=res['model'])

r=r + [res['prediction'].equals(res1['prediction'])]

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE", model=res['model'])

r=r + [res['prediction']['final'].equals(res1['prediction']['final'])]

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE", model=res['model'])

r=r + [res['prediction']['final'].equals(res1['prediction']['final'])]

#USECASE 2

Data=pd.read_csv("data/Data_sim_1000.csv",sep=",",header=None)
Data.columns=["target"]+["X"+str(x) for x in range(1,len(Data.columns))]

res=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE")
res1=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE", model=res['model'])

r=r + [res['prediction'].equals(res1['prediction'])]

res=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE")
res1=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE",model=res['model'])

r=r + [res['prediction'].equals(res1['prediction'])]

res=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE")
res1=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE",model=res['model'])

r=r + [res['prediction']['final'].equals(res1['prediction']['final'])]

res=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE")
res1=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE",model=res['model'])

r=r + [res['prediction']['final'].equals(res1['prediction']['final'])]

#USECASE 3

Data=pd.read_csv("data/Data_trading.csv",sep=";")

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE", model=res['model'])

r=r + [res['prediction'].equals(res1['prediction'])]

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<-", loss_function="MSE", model=res['model'])

r=r + [res['prediction'].equals(res1['prediction'])]

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MSE", model=res['model'])

r=r + [res['prediction']['final'].equals(res1['prediction']['final'])]

res=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE")
res1=regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE", model=res['model'])

r=r + [res['prediction']['final'].equals(res1['prediction']['final'])]


print(r)

