//rego: Automatic time series forecasting and missing value imputation.
//
//Copyright (C) Davide Altomare and David Loris <https://channelattribution.io>
//
//This source code is licensed under the MIT license found in the
//LICENSE file in the root directory of this source tree. 

// #ifndef CALIB_H
// #define CALIB_H

#define language_cpp
//#define language_python
//#define language_R
#include <armadillo>

using namespace std;
using namespace arma;

#define uli unsigned long int
using svec1=vector<double>;
using svec2=vector< vector <double> >;
using svec3=vector < vector< vector <double> > >;
using svec4=vector < vector < vector< vector <double> > > >;
//https://stackoverflow.com/questions/44663890/c-vectorvectordouble-using-typename-alias

struct str_model_selection
{
 
 field<vec> models;
 mat predictions;

};


struct str_output
{

 mat predictions;
 vec performances;
 
 mat fw_predictions;
 field<vec> fw_models;
 vec fw_performances;
 
 mat bw_predictions;
 field<vec> bw_models;
 vec bw_performances;

};


str_output regpred_cpp(mat* Y, double from_lag, double max_lag, double alpha, uli nsim, bool flg_print, string direction, string loss_function, bool pred_only, bool flg_const, bool flg_diff, double h_c, vector<field<vec>> models);

#ifdef language_py 
pair < svec3, pair < svec2, svec4 > >
regpred_py(svec2& Y, double from_lag, double max_lag, double alpha, uli nsim, bool flg_print, string direction, string loss_function, bool pred_only, bool flg_const, bool flg_diff, double h_c, svec4& models);
#endif

#ifdef language_R
 RcppExport SEXP regpred_R(SEXP Y_p, SEXP from_lag_p, SEXP max_lag_p, SEXP alpha_p, SEXP nsim_p, SEXP flg_print_p, SEXP direction_p, SEXP loss_function_p, SEXP pred_only_p, SEXP flg_const_p, SEXP flg_diff_p, SEXP h_c_p, SEXP vmodels_p);
#endif
