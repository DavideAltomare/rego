//rego: Automatic time series forecasting and missing value imputation.
//
//Copyright (C) Davide Altomare and David Loris <channelattribution.io>
//
//This source code is licensed under the MIT license found in the
//LICENSE file in the root directory of this source tree. 

// #ifndef CALIB_H
// #define CALIB_H

#define language_cpp
//#define language_python
//#define language_Rgo
#include <armadillo>

#define uli unsigned long int

using namespace std;
using namespace arma;

struct str_output
{

 mat predictions;
 double L=datum::nan;   
 double L_adj=datum::nan;
 
 mat fw_predictions;
 vector<uli> fw_var_x_idx;
 vector<uli> fw_var_ar_idx;
 vector<uli> fw_var_ma_idx;
 double fw_L=datum::nan;   
 double fw_L_adj=datum::nan;

 mat bw_predictions;
 vector<uli> bw_var_x_idx;
 vector<uli> bw_var_ar_idx;
 vector<uli> bw_var_ma_idx;
 double bw_L=datum::nan;
 double bw_L_adj=datum::nan;

};


str_output regpred_cpp(Mat<double>* Y, double max_lag, double alpha, uli nsim, bool flg_print, string direction);

#ifdef language_py 
pair < pair < list< vector<uli> >, list< vector<string> > > , pair < list< vector< vector<double> > >, list< double> > >
regpred_py(vector< vector<double> >& Y, double max_lag, double alpha, uli nsim, bool flg_print, string direction);
#endif

#ifdef language_R
 RcppExport SEXP regpred_R(SEXP Y_p, SEXP max_lag_p, SEXP alpha_p, SEXP nsim_p, SEXP flg_print_p, SEXP direction_p);
#endif
