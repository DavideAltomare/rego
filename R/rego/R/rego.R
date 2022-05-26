# rego: Automatic Time Series Forecasting and Missing Value Imputation
#
# Copyright (C) Davide Altomare and David Loris <https://channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

.v=as.character(packageVersion("rego"))
  
.onAttach = function(libname, pkgname) {
  
 packageStartupMessage("Visit https://channelattribution.io/docs/rego for more information about rego")
 packageStartupMessage(paste0("rego ",.v))

}

regpred=function(Data, from_lag=1, max_lag="auto", alpha=0.05, nsim=1000, flg_print=1, direction="->", loss_function="MSE", flg_const=TRUE, flg_diff=FALSE, model=NULL){

    if(!("data.frame"%in%class(Data)|"data.table"%in%class(Data))){
     stop("Data must be a data.frame or a data.table")
    } 
     
    
    if(from_lag < 1){
       stop("from_lag must be > 1")
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

    if(!loss_function %in% c("MSE","MAE")){ 
       stop("loss_function must be 'MSE' or 'MAE'")
    }
    
    if(!flg_print %in% c(0,1)){
       stop("flg_print must be 0 or 1")
    }
	
    if(length(max_lag)==0){
        max_lag=0
    }
    if(max_lag=="auto"){
        max_lag=-1
    }
    
    if(!flg_const %in% c(0,1)){
       stop("flg_const must be 0 or 1")
    }
    
    if(!flg_diff %in% c(0,1)){
       stop("flg_diff must be 0 or 1")
    }

    if(length(model)==0){
       pred_only=0
       model=list()
    }else{
       pred_only=1
       if(direction=="<->"){
         fw_model=model[['forward']]
         bw_model=model[['backward']]
       }else if(direction=="->"){
         fw_model=model
         model=list()
         model[['forward']]=fw_model
         model[['backward']]=fw_model
       }else if(direction=="<-"){
         bw_model=model
         model=list()
         model[['forward']]=bw_model
         model[['backward']]=bw_model
       }
    }
      
    cols_Y=colnames(Data)
    Y=as.matrix(Data)
    rm(Data)
    
    res0=.Call("regpred_R", Y , from_lag, max_lag, alpha, nsim, flg_print, direction, loss_function, pred_only, flg_const, flg_diff, model)
 
    res=list()
 
    prediction=as.data.frame(res0$prediction$final)
    colnames(prediction)=c('real', 'fitted', 'lower_bound','predicted','upper_bound')
    prediction=prediction[,c('real', 'lower_bound','predicted','upper_bound')]
    
    if(direction=="<->"){
      fw_prediction=as.data.frame(res0$prediction$forward)
	   if(nrow(fw_prediction)>0){
	     colnames(fw_prediction)=c('real','fitted', 'lower_bound','predicted','upper_bound')
        fw_prediction=fw_prediction[,c('real', 'lower_bound','predicted','upper_bound')]
      }
    
      bw_prediction=as.data.frame(res0$prediction$backward)
	   if(nrow(bw_prediction)>0){
	     colnames(bw_prediction)=c('real','fitted', 'lower_bound','predicted','upper_bound')
        bw_prediction=bw_prediction[,c('real', 'lower_bound','predicted','upper_bound')]
      }
    }
    
    if(pred_only==0){
        fw_model=res0$model$forward
        bw_model=res0$model$backward
    }
    
    if(direction=="<->"){
      res$prediction$final=prediction
      res$prediction$forward=fw_prediction
      res$prediction$backward=bw_prediction
      res$model$forward=fw_model
      res$model$backward=bw_model
    }else if(direction=="->"){
       res$prediction=prediction
       res$model=fw_model
    }else if(direction=="<-"){
       res$prediction=prediction
       res$model=bw_model
    }
   
    return(res)
	
}

