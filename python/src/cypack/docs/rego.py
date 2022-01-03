
"""

**Automatic Time Series Forecasting and Missing Value Imputation.**
rego is intended for predicting and imputing time series. Its algorithm can automatically set all the parameters needed, thus in the minimal configuration it only requires the target variable and the regressors if present. It can address large problems with hundreds or thousands of dependent variables and problems in which the number of dependent variables is greater than the number of observations. Moreover it can be used not only for time series but also for any other real valued target variable. The algorithm implemented includes a Bayesian stochastic search methodology for model selection and a robust estimation method based on bootstrapping. rego is fast because all the code is C++. 

"""
        
def regpred(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->"):

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
    direction : int, default 10
        if "->" then only a forward prediction will be executed, if "<-" then only a backward prediction will be executed, if "<->" then both a forward than a backward prediction will be executed if possible. For imputing missing values is convenient to leave default "<->".        
            
    Returns
    -------
    dictionary
        final : final predictions
            (DataFrame) predictions : predictions and confidence interval
            (float) R2 : coefficient of determination
            (float) R2_adj : adjusted coefficient of determination
            (float) abs_err_ratio : mean absolute error of the model selected divided by the mean absolute error of the trivial predictor (average of the observations)  
        forward : forward predictions
            (DataFrame) predictions : predictions and confidence interval
            (float) R2 : coefficient of determination
            (float) R2_adj : adjusted coefficient of determination
            (float) abs_err_ratio : mean absolute error of the model selected divided by the mean absolute error of the trivial predictor (average of the observations) 
        backward: backward predictions
            (DataFrame) predictions : predictions and confidence interval
            (float) R2 : coefficient of determination
            (float) R2_adj : adjusted coefficient of determination
            (float) abs_err_ratio : mean absolute error of the model selected divided by the mean absolute error of the trivial predictor (average of the observations) 
                
                        
    Examples
    --------
    
    >>> import pandas as pd    
    >>> from rego import *    

    Example 1: seasonal time series

    >>> Data=pd.read_csv("https://app.channelattribution.net/data/Data_air.csv",header=None)

    >>> res=regpred(Data)

    print final prediction 
    
    >>> print(res['final']['predictions'])

    print final R2_adj

    >>> print(res['final']['R2_adj'])                                                

    Example 2: high dimensional problem

    >>> Data=pd.read_csv("https://app.channelattribution.net/data/Data_sim_1000.csv",header=None)

    >>> res=regpred(Data, max_lag=None)

    print final prediction 
    
    >>> print(res['final']['predictions'])

    print final R2_adj

    >>> print(res['final']['R2_adj'])                                                

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
    
    if (max_lag == None):
        max_lag=0

    if (max_lag == "auto"):
        max_lag=-1    
        
    if(str(type(Data))!="<class 'pandas.core.frame.DataFrame'>"):
        raise NameError("Data must be a Dataframe'")
    
    cols_Y=[str(x) for x in Data.columns.tolist()]
    Y=Data.to_numpy()
    del Data
                    
    res0=__regpred_py(Y, [x.encode('utf-8') for x in cols_Y], max_lag, alpha, nsim, flg_print, direction.encode('utf-8'))
    
    res={'final':{},'forward':{},'backward':{}}

    res['final']['predictions']=pd.DataFrame(res0[1][0][0]) 
    res['final']['predictions'].columns=['real','lower_bound','predicted','upper_bound'] 
    res['final']['R2']=res0[1][1][0]
    res['final']['R2_adj']=res0[1][1][1]
    res['final']['abs_err_ratio']=res0[1][1][2]

    res['forward']['predictions']=pd.DataFrame(res0[1][0][1]) 
    if(len(res['forward']['predictions'])>0):
        res['forward']['predictions'].columns=['real','lower_bound','predicted','upper_bound'] 
    res['forward']['I_step_var_x_idx']=res0[0][0][0] 
    res['forward']['I_step_var_x_idx']=res0[0][0][0] 
    res['forward']['I_step_var_ar_idx']=res0[0][0][1]
    res['forward']['I_step_var_ma_idx']=res0[0][0][2]
    res['forward']['II_step_var_x_idx']=res0[0][0][3] 
    res['forward']['II_step_var_ar_idx']=res0[0][0][4]
    res['forward']['II_step_var_ma_idx']=res0[0][0][5]
    res['forward']['R2']=res0[1][1][3]
    res['forward']['R2_adj']=res0[1][1][4]
    res['forward']['abs_err_ratio']=res0[1][1][5]

    res['backward']['predictions']=pd.DataFrame(res0[1][0][2]) 
    if(len(res['backward']['predictions'])>0):
        res['backward']['predictions'].columns=['real','lower_bound','predicted','upper_bound'] 
    res['backward']['bw_I_step_var_x_idx']=res0[0][0][6] 
    res['backward']['bw_I_step_var_ar_idx']=res0[0][0][7]
    res['backward']['bw_I_step_var_ma_idx']=res0[0][0][8]
    res['backward']['bw_II_step_var_x_idx']=res0[0][0][9] 
    res['backward']['bw_II_step_var_ar_idx']=res0[0][0][10]
    res['backward']['bw_II_step_var_ma_idx']=res0[0][0][11]
    res['backward']['R2']=res0[1][1][6]
    res['backward']['R2_adj']=res0[1][1][7]
    res['backward']['abs_err_ratio']=res0[1][1][8]
    
    return(res)

    
