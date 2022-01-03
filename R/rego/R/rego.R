# rego: Automatic Time Series Forecasting and Missing Value Imputation
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

.v=as.character(packageVersion("rego"))
  
.onAttach = function(libname, pkgname) {
  
 packageStartupMessage("Visit https://www.channelattribution.net/docs/rego for more information about rego")
 packageStartupMessage(paste0("rego ",.v))

}

regpred=function(Data, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="<->"){

    if(!("data.frame"%in%class(Data)|"data.table"%in%class(Data))){
     stop("Data must be a data.frame or a data.table")
    } 
     
    if(length(max_lag)>0){ 
      if(max_lag!="auto"){
        max_lag=as.numeric(max_lag)
        if(max_lag < 0){
           stop("max_lag must be >= 0")
        }
	   }
    }
	
    if((alpha < 0) | (alpha > 1)){
       stop("alpha must be in [0,1]")
    }

    if(nsim < 0){
       stop("nsim must be > 0")
    }

    if(!flg_print %in% c(0,1)){
       stop("flg_print must be 0 or 1")
    }

    if(!direction %in% c("->","<-","<->")){ 
       stop("direction must be '->', '<-' or '<->'")
    }
	
    if(length(max_lag)==0){
        max_lag=0
    }
    if(max_lag=="auto"){
        max_lag=-1
    }	
		
    cols_Y=colnames(Data)
    Y=as.matrix(Data)
    rm(Data)
    
    res=.Call("regpred_R", Y , max_lag, alpha, nsim, flg_print, direction)
 
    res$final$predictions=as.data.frame(res$final$predictions)
    colnames(res$final$predictions)=c('real', 'fitted', 'upper_bound','predicted','lower_bound')
    res$final$predictions=res$final$predictions[,c('real', 'upper_bound','predicted','lower_bound')]

	 res$forward$predictions=as.data.frame(res$forward$predictions)
	 if(nrow(res$forward$predictions)>0){
	   colnames(res$forward$predictions)=c('real','fitted', 'upper_bound','predicted','lower_bound')
      res$forward$predictions=res$forward$predictions[,c('real', 'upper_bound','predicted','lower_bound')]
      if(length(res$forward$var_x_names)>0){
        res$forward$var_x_names=cols_Y[res$forward$var_x_names]  
      }
    }

    res$backward$predictions=as.data.frame(res$backward$predictions) 
	 if(nrow(res$backward$predictions)>0){
	   colnames(res$backward$predictions)=c('real', 'fitted', 'upper_bound','predicted','lower_bound')
      res$backward$predictions=res$backward$predictions[,c('real', 'upper_bound','predicted','lower_bound')]
      if(length(res$backward$var_x_names)>0){
       res$backward$var_x_names=cols_Y[res$backward$var_x_names]  
      }
    }
	  
    return(res)
	
}

