# rego: Automatic time series forecasting and missing values imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

#!/bin/bash

main=$1

rm -rf /tmp/build
mkdir -p /tmp/build
cp -r c++ /tmp/build 
cp -r python /tmp/build 
cp ${main} /tmp/build/c++/test.cpp

dir_princ=$PWD

cd /tmp/build/c++
g++ -I"${dir_princ}/python/src/cypack/armadillo/include" -std=c++11 functions.cpp test.cpp -o test.o -pthread  -llapack -lblas

