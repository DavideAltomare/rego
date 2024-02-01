//rego: Automatic time series forecasting and missing values imputation.
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
#include <ctime>

// #define OPTIM_ENABLE_ARMA_WRAPPERS
// #include "optim.hpp"

#include <armadillo>

#define uli unsigned long int

#include "functions.h"

using namespace std;
using namespace arma;

#define language_cpp 

int main(void) 
{
  
///AIR
 
 mat Y,Y0; 
 Y0.load("data/Data_air.csv", csv_ascii);
 mat Mnan(12,Y0.n_cols);
 Mnan.fill(datum::nan);
 Y0=join_cols(Y0,Mnan);
 vec vnan=zeros<vec>(12)+datum::nan;
 Y=Y0;
 Y.submat(59,0,70,0)=vnan;
 double from_lag=1;
 double max_lag=-1; //-1 means auto
 string direction="->";
 double alpha=0.05;
 uli nsim=1000;
 bool flg_print=1;
 string loss_function="MSE";
 bool pred_only=0;
 bool flg_const=1;
 bool flg_diff=1;
 vector < field<vec> > vmodels;
 double h_c=1;
 
  //Y.print("Y_0");
 
 str_output out=regpred_cpp(&Y, from_lag, max_lag, alpha, nsim, flg_print, direction, loss_function, pred_only, flg_const, h_c, flg_diff, vmodels);
 
 out.predictions.print("predictions");
 
 Y.print("Y_1");
 
 out=regpred_cpp(&Y, from_lag, max_lag, alpha, nsim, flg_print, direction, loss_function, pred_only, flg_const, h_c, flg_diff, vmodels);
 
 out.predictions.print("predictions");
 
 return 0;

}