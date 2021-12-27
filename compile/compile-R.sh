# rego: Automatic time series forecasting and missing values imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

#!/bin/bash

rm -rf /tmp/build/R
mkdir -p /tmp/build/R
cp -r R /tmp/build

cp c++/functions.cpp  /tmp/build/R/rego/src/rego.cpp
cp c++/functions.h  /tmp/build/R/rego/src/rego.h

sed -i 's/\#define language_cpp/\/\/\#define language_cpp/' /tmp/build/R/rego/src/rego.cpp
sed -i 's/\#define language_cpp/\/\/\#define language_cpp/' /tmp/build/R/rego/src/rego.h

sed -i 's/\/\/\#define language_R/\#define language_R/' /tmp/build/R/rego/src/rego.cpp
sed -i 's/\/\/\#define language_R/\#define language_R/' /tmp/build/R/rego/src/rego.h

cd /tmp/build/R

rm -f rego/src/rego.so rego/src/rego.o rego/src/init.o 
find /tmp/build/R/  -name "rego*.tar.gz" -delete

R CMD REMOVE rego
R CMD INSTALL --build rego



