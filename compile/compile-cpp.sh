# rego: Automatic time series forecasting and missing value imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

#!/bin/bash

rm -rf /tmp/build
mkdir -p /tmp/build
cp -r c++ /tmp/build 
cp -r python /tmp/build 
cp /home/Projects/GIT/personal/rego/scripts/main.cpp /tmp/build/c++/test.cpp

dir_princ=$PWD

cd /tmp/build/c++
dir_optim=/home/Projects/GIT/personal/channelattributionpro/src
g++ -I"${dir_princ}/python/src/cypack/armadillo/include" -I${dir_optim}"/c++/optim-master/header_only_version/" -std=c++11 functions.cpp test.cpp -o test.o -pthread  -llapack -lblas

