//rego: Automatic time series forecasting and missing value imputation.
//
//Copyright (C) Davide Altomare and David Loris <channelattribution.io>
//
//This source code is licensed under the MIT license found in the
//LICENSE file in the root directory of this source tree. 


#include <iostream>
#include <vector>
#include <set>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <sstream>
#include <list>
#include <string>
#include <random>
#include <numeric>
#include <time.h> 
#include <thread>
#include <map>
#include <algorithm>
#include <list>
#include <limits> 

#include <armadillo>

#define uli unsigned long int

#include "functions.h"

using namespace std;
using namespace arma;


#define language_cpp 

int main(void) 
{

   
 mat Y;
 vec vnan;
 double max_lag,alpha;
 string direction;
 uli nsim;
 bool flg_print;
 string loss_function;
 bool pred_only=0;
 vector < field<vec> > vmodels;
 str_output out;
 

 //USECASE 1
  
 Y.load("data/Data_air.csv", csv_ascii);
 Y.replace(0,datum::nan); 

 max_lag=-1; //-1 means auto
 direction="<->";
 alpha=0.05;
 nsim=1000;
 flg_print=1;
 loss_function="MSE";
 
 out=regpred_cpp(&Y, max_lag, alpha, nsim, flg_print, direction, loss_function, pred_only,vmodels);

 out.predictions.print("Usecase 1 predictions");

 //USECASE 2

 Y.load("data/Data_sim_1000.csv", csv_ascii);
 Y.replace(0,datum::nan); 

 max_lag=-1; //-1 means auto
 direction="<->";
 alpha=0.05;
 nsim=1000;
 flg_print=1;
 loss_function="MSE";
 

 out=regpred_cpp(&Y, max_lag, alpha, nsim, flg_print, direction, loss_function,pred_only,vmodels);

 out.predictions.print("Usecase 2 predictions");


 

 return 0;

}