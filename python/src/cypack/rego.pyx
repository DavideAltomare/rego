# rego: Automatic Time Series Forecasting and Missing Value Imputation
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.list cimport list
from libcpp.pair cimport pair
import pandas as pd


__version="1.3.2"

print("Visit https://www.channelattribution.net/docs/rego for more information about rego")
print("Version: " + str(__version))
    
cdef extern from "functions.h":

    pair[pair[list[vector[unsigned long int]],list[vector[string]]],pair[list[vector[vector[double]]],list[double]]]  regpred_py(vector[vector[double]]& Y, double max_lag, double alpha, unsigned long int nsim, int flg_print, string direction, string loss_function);


def __regpred_py(vector[vector[double]] Y, double max_lag, double alpha, unsigned long int nsim, int flg_print, string direction, string loss_function):
    return(regpred_py(Y,max_lag,alpha,nsim,flg_print,direction,loss_function))


#start documentation

"""

**Automatic time series forecasting and missing values imputation.**
rego is a machine learning algorithm for predicting and imputing time series. It can automatically set all the parameters needed, thus in the minimal configuration it only requires the target variable and the dependent variables if present. It can address large problems with hundreds or thousands of dependent variables and problems in which the number of dependent variables is greater than the number of observations. Moreover it can be used not only for time series but also for any other real valued target variable. The algorithm implemented includes a Bayesian stochastic search methodology for model selection and a robust estimation based on bootstrapping. rego is fast because all the code is C++.

"""
        
def regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->", loss_function="MAE"):

    '''
    
    Parameters
    ----------
    Data : DataFrame
        data.frame containing target variable at first column and regressors if present from second to last column.
    max_lag: string
        maximum time lag to be considered in the autoregressive moving average part of the algorithm. If "auto" then the algorithm will set a suitable value. Set to 0 or None if you want to remove the autoregressive moving average part as in case of non time series data.
    alpha : string
        significance level for the confidence interval produced around predictions. If 0.05 then the algorithm will calculate a 95\% confidence interval around predictions.
    nsim : string
        number of bootstrap replications used for producing confidence interval around predictions.
    flg_print : string, optional, default None
        if 1 some information during the evaluation will be printed.
    direction : string, default "<->"
        if "->" then only a forward prediction will be executed, if "<-" then only a backward prediction will be executed, if "<->" then both a forward than a backward prediction will be executed if possible. For imputing missing values is convenient to leave default "<->".        

    loss_function : string, default "MAE"
        if "MAE" then mean absolute error is used as penalty function in regressions, if "MSE" then mean squared error is used as penalty function in regressions        


    Returns
    -------
    dictionary
        final : final predictions
            (DataFrame) predictions : predictions and confidence interval
            (float) L : mean absolute error of the model selected divided by the mean absolute error of the trivial predictor (average of the observations) 
            (float) L_adj : L penalized by the number of coefficients 
        forward : forward predictions
            (DataFrame) predictions : predictions and confidence interval
            (list) var_x_names : names of the regressors selected by the algorithm
            (list) var_ar_idx : AR retards selected by the algorithm
            (list) var_ma_idx : MA retards selected by the algorithm
            (float) L : mean absolute error of the model selected divided by the mean absolute error of the trivial predictor (average of the observations) 
            (float) L_adj : L penalized by the number of coefficients
        backward: backward predictions
            (DataFrame) predictions : predictions and confidence interval
            (list) var_x_names : names of the regressors selected by the algorithm
            (list) var_ar_idx : AR retards selected by the algorithm
            (list) var_ma_idx : MA retards selected by the algorithm
            (float) L : mean absolute error of the model selected divided by the mean absolute error of the trivial predictor (average of the observations) 
            (float) L_adj : L penalized by the number of coefficients
                
                        
    Examples
    --------
    
    >>> import pandas as pd    
    >>> from rego import *    

    Example 1: seasonal time series

    >>> Data=pd.read_csv("https://app.channelattribution.net/data/Data_air.csv",header=None)

    >>> res=regpred(Data)

    print final prediction 
    
    >>> print(res['final']['predictions'])

    print final L_adj

    >>> print(res['final']['L_adj'])                                                

    Example 2: high dimensional problem

    >>> Data=pd.read_csv("https://app.channelattribution.net/data/Data_sim_1000.csv",header=None)

    >>> res=regpred(Data, max_lag=None)

    print final prediction 
    
    >>> print(res['final']['predictions'])

    print final L_adj

    >>> print(res['final']['L_adj'])                                                

    '''


    if max_lag != None:
        if max_lag!="auto":
            max_lag=float(max_lag)
            if (max_lag < 0):
                raise NameError("max_lag must be >= 0")
        
    if ((alpha < 0) | (alpha > 1)):
       raise NameError("alpha must be in [0,1]")

    if (nsim < 0):
       raise NameError("nsim must be > 0")

    if (flg_print not in [0,1]):
       raise NameError("flg_print must be 0 or 1")

    if (direction not in ["->","<-","<->"]):
       raise NameError("direction must be '->', '<-' or '<->'")

    if (loss_function not in ["MAE","MSE"]):
       raise NameError("loss_function must be 'MAE' or 'MSE'")
    
    if (max_lag == None):
        max_lag=0

    if (max_lag == "auto"):
        max_lag=-1    
        
    if(str(type(Data))!="<class 'pandas.core.frame.DataFrame'>"):
        raise NameError("Data must be a Dataframe'")
    
    cols_Y=[str(x) for x in Data.columns.tolist()]
    Y=Data.to_numpy()
    del Data
                    
    res0=__regpred_py(Y, max_lag, alpha, nsim, flg_print, direction.encode('utf-8'), loss_function.encode('utf-8'))
    
    res={'final':{},'forward':{},'backward':{}}

    res['final']['predictions']=pd.DataFrame(res0[1][0][0]) 
    res['final']['predictions'].columns=['real','fitted','lower_bound','predicted','upper_bound'] 
    del res['final']['predictions']['fitted']
    res['final']['L']=res0[1][1][0]
    res['final']['L_adj']=res0[1][1][1]
    
    res['forward']['predictions']=pd.DataFrame(res0[1][0][1]) 
    if(len(res['forward']['predictions'])>0):
        res['forward']['predictions'].columns=['real','fitted','lower_bound','predicted','upper_bound'] 
        del res['forward']['predictions']['fitted']
    res['forward']['var_x_names']=res0[0][0][0] 
    if len(res['forward']['var_x_names']):
        res['forward']['var_x_names']=(pd.Series(cols_Y)[res['forward']['var_x_names']]).tolist()
    res['forward']['var_ar_idx']=res0[0][0][1]
    res['forward']['var_ma_idx']=res0[0][0][2]
    res['forward']['L']=res0[1][1][2]
    res['forward']['L_adj']=res0[1][1][3]

    res['backward']['predictions']=pd.DataFrame(res0[1][0][2]) 
    if(len(res['backward']['predictions'])>0):
        res['backward']['predictions'].columns=['real','fitted','lower_bound','predicted','upper_bound'] 
        del res['backward']['predictions']['fitted']
    res['backward']['var_x_names']=res0[0][0][3] 
    if len(res['backward']['var_x_names']):
        res['backward']['var_x_names']=(pd.Series(cols_Y)[res['backward']['var_x_names']]).tolist()
    res['backward']['var_ar_idx']=res0[0][0][4]
    res['backward']['var_ma_idx']=res0[0][0][5]

    res['backward']['L']=res0[1][1][4]
    res['backward']['L_adj']=res0[1][1][5]
    
    return(res)

    
