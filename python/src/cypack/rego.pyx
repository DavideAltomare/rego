# rego: Automatic Time Series Forecasting and Missing Value Imputation
#
# Copyright (C) Davide Altomare and David Loris <https://channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.list cimport list
from libcpp.pair cimport pair
import pandas as pd


__version="1.6.1"

print("Visit https://channelattribution.io/docs/rego for more information about rego")
print("Version: " + str(__version))
    
#vector[vector[vector[vector[double]]]]
#pair < vec3, pair < vec2, vec4 > >


    
cdef extern from "functions.h":

    pair[vector[vector[vector[double]]],pair[vector[vector[double]],vector[vector[vector[vector[double]]]]]] regpred_py(vector[vector[double]]& Y, double from_lag, double max_lag, double alpha, unsigned long int nsim, int flg_print, string direction, string loss_function, int pred_only, int flg_const, int flg_diff, double h_c, vector[vector[vector[vector[double]]]]& vmodels);


def __regpred_py(vector[vector[double]] Y, double from_lag, double max_lag, double alpha, unsigned long int nsim, int flg_print, string direction, string loss_function, int pred_only, int flg_const, int flg_diff, double h_c, vector[vector[vector[vector[double]]]] vmodels):
    return(regpred_py(Y,from_lag,max_lag,alpha,nsim,flg_print,direction,loss_function,pred_only,flg_const,flg_diff,h_c,vmodels))


#start documentation

"""

**Automatic time series forecasting and missing values imputation.**
rego is a machine learning algorithm for predicting and imputing time series. It can automatically set all the parameters needed, thus in the minimal configuration it only requires the target variable and the dependent variables if present. It can address large problems with hundreds or thousands of dependent variables and problems in which the number of dependent variables is greater than the number of observations. Moreover it can be used not only for time series but also for any other real valued target variable. The algorithm implemented includes a Bayesian stochastic search methodology for model selection and a robust estimation based on bootstrapping. rego is fast because all the code is C++.

"""
        
def regpred(Data, from_lag=1, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", flg_const=True, flg_diff=False, model=None):

    '''
    
    Parameters
    ----------
    Data : DataFrame
        data.frame containing target variable at first column and regressors if present from second to last column.
    from_lag: int
        minimum time lag to be considered in the autoregressive moving average part of the algorithm.
    max_lag: string
        maximum time lag to be considered in the autoregressive moving average part of the algorithm. If "auto" then the algorithm will set a suitable value. Set to 0 or None if you want to remove the autoregressive moving average part as in case of non time series data.
    alpha : double
        significance level for the confidence interval produced around predictions. If 0.05 then the algorithm will calculate a 95\% confidence interval around predictions.
    nsim : int
        number of bootstrap replications used for producing confidence interval around predictions.
    flg_print : bool, optional, default None
        if 1 some information during the evaluation will be printed.
    direction : string, default "->"
        if "->" then only a forward prediction will be executed, if "<-" then only a backward prediction will be executed, if "<->" then both a forward than a backward prediction will be executed if possible. For imputing missing values is convenient to leave default "<->".        
    flg_const : bool, default True
        if True then a constant is included into the model
    flg_diff : bool, default False
        if True and no regressor is present then if the target variable exhibits a trend, it is one-step differentiated up to two times
    model: list
        estimated models from a previous train to be used in new data prediction without retraining


    Returns
    -------
    dictionary

        (DataFrame) predictions : predictions and confidence interval
        (DataFrame) models : estimated models
    
    Examples
    --------
    
    >>> import pandas as pd
    >>> from rego import *

    Example 1: seasonal time series

    >>> Data=pd.read_csv("https://channelattribution.io/csv/Data_air.csv",header=None)

    >>> res=regpred(Data)

    print final prediction 
    
    >>> print(res['predictions'])

    Example 2: high dimensional problem

    >>> Data=pd.read_csv("https://channelattribution.io/csv/Data_sim_1000.csv",header=None)

    >>> res=regpred(Data, max_lag=None)

    print final prediction 
    
    >>> print(res['predictions'])

    '''

    loss_function="MSE" 
    h=None 

    if (from_lag < 1):
       raise NameError("from_lag must be > 1")

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
       
    if (flg_const not in [0,1]):
       raise NameError("flg_const must be 0 or 1")
       
    if (flg_diff not in [0,1]):
       raise NameError("flg_diff must be 0 or 1")

    if h != None:
        if (type(h) != int):
            raise NameError("flg_diff must be 0 or 1")
    else:
        h=-1
    
    if (max_lag == None):
        max_lag=0

    if (max_lag == "auto"):
        max_lag=-1
    
    if (model == None):
        pred_only=0
        model=[]
    else:
        pred_only=1
        if direction=="<->":
            fw_model=model['forward']
            bw_model=model['backward']
        elif direction=="->":
            fw_model=model
            bw_model=[]
        elif direction=="<-":
            fw_model=[]
            bw_model=model
        model=[fw_model, bw_model]
            
    if(str(type(Data))!="<class 'pandas.core.frame.DataFrame'>"):
        raise NameError("Data must be a Dataframe'")
    
    cols_Y=[str(x) for x in Data.columns.tolist()]
    Y=Data.to_numpy()
    del Data
                    
    res0=__regpred_py(Y, from_lag, max_lag, alpha, nsim, flg_print, direction.encode('utf-8'), loss_function.encode('utf-8'), pred_only, flg_const, flg_diff, h, model)
    
    prediction=pd.DataFrame(res0[0][0])
    prediction.columns=['real','fitted','lower_bound','predicted','upper_bound']
    del prediction["fitted"]
    
    if direction=="<->":
        fw_prediction=pd.DataFrame(res0[0][1])
        if(len(fw_prediction)>0):
            fw_prediction.columns=['real','fitted','lower_bound','predicted','upper_bound']
            del fw_prediction["fitted"]
        
        bw_prediction=pd.DataFrame(res0[0][2])
        if(len(bw_prediction)>0):
            bw_prediction.columns=['real','fitted','lower_bound','predicted','upper_bound']
            del bw_prediction["fitted"]
            
    if pred_only==0:
        fw_model=res0[1][1][0]
        bw_model=res0[1][1][1]
    
    res={'prediction':{}, 'model':{}}
    
    if direction=="<->":
        res['prediction']['final']=prediction
        res['prediction']['forward']=fw_prediction
        res['prediction']['backward']=bw_prediction
        res['model']['forward']=fw_model
        res['model']['backward']=bw_model
    elif direction=="->":
        res['prediction']=prediction
        res['model']=fw_model
    elif direction=="<-":
        res['prediction']=prediction
        res['model']=bw_model
        
    # if 0!=0:
        
    #     performance=pd.DataFrame(res0[1][0][0])
    #     performance.columns=['performance']
    #     performance.index=["L","L_adj"]
    #     fw_performance=pd.DataFrame(res0[1][0][1])
    #     fw_performance.columns=['performance']
    #     fw_performance.index=["L","L_adj"]
    #     bw_performance=pd.DataFrame(res0[1][0][2])
    #     if(len(bw_performance)>0):
    #         bw_performance.columns=['performance']
    #         bw_performance.index=["L","L_adj"]
        
    #     res['performance']['final']=performance
    #     res['performance']['forward']=fw_performance
    #     res['performance']['backward']=bw_performance
        
    return(res)

    
