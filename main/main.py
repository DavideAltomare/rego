# rego: Automatic time series forecasting and missing values imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 


import pandas as pd
from rego import *
import plotly.express as px


#usecase 1

Data=pd.read_csv("https://app.channelattribution.net/data/Data_air.csv",sep=";",header=None)

res=regpred(Data, max_lag="auto", alpha=0.05,nsim=1000, flg_print=1, direction="<->")

pred=res["prediction"]["final"]
fig = px.line(pred)
fig.show()

#use the saved model to get predictions
res1=regpred(Data, max_lag="auto", alpha=0.05,nsim=1000, flg_print=1, direction="<->", model=res['model'])

#usecase 2

Data=pd.read_csv("https://app.channelattribution.net/data/Data_sim_1000.csv",sep=",",header=None)
Data.columns=["target"]+["X"+str(x) for x in range(1,len(Data.columns))]

res=regpred(Data, max_lag=None, alpha=0.05, nsim=1000, flg_print=1, direction="<->")

pred=res["prediction"]["final"]
fig = px.line(pred)
fig.show()