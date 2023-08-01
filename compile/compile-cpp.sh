# rego: Automatic time series forecasting and missing value imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

#!/bin/bash

dir_rego=/home/Projects/GIT/personal/rego/repo

file_main=$1

rm -rf /tmp/build
mkdir -p /tmp/build
cp -r c++ /tmp/build  
cp $file_main /tmp/build/c++/main.cpp

dir_princ=$PWD

cd /tmp/build/c++
dir_optim=${dir_rego}/libs/optim-master/header_only_version
dir_armadillo=${dir_rego}/libs/armadillo/include
# g++ -I${dir_armadillo} -I${dir_optim}"/header_only_version/" -std=c++11 functions.cpp main.cpp -o main.o -pthread  -llapack -lblas
g++ -I${dir_armadillo} -I${dir_optim} -std=c++11 functions.cpp main.cpp -o main.o -pthread


