# rego: Automatic time series forecasting and missing values imputation.
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

#!/bin/bash

rm -rf /tmp/build/python
mkdir -p /tmp/build/python
cp -r python /tmp/build 

cp c++/functions.cpp  /tmp/build/python/src/cypack/
cp c++/functions.h  /tmp/build/python/src/cypack/

sed -i 's/\#define language_cpp/\/\/\#define language_cpp/' /tmp/build/python/src/cypack/functions.cpp
sed -i 's/\#define language_cpp/\/\/\#define language_cpp/' /tmp/build/python/src/cypack/functions.h

sed -i 's/\/\/\#define language_python/\#define language_py/' /tmp/build/python/src/cypack/functions.cpp
sed -i 's/\/\/\#define language_python/\#define language_py/' /tmp/build/python/src/cypack/functions.h

#cd /tmp/build/python
#python3 setup.py build_ext --inplace
cd /tmp/build/
pip3 install python/



