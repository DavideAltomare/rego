# rego: Automatic time series forecasting and missing value imputation.
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

if [ "$1" = "test" ]; then
    now=$(date +"%Y-%m-%d %H:%M:%S")
    sed -i "1s/^/print('${now}')\n/" /tmp/build/python/src/cypack/rego.pyx
fi

# pip3 uninstall -y rego
# cd /tmp/build/python
# python3 setup.py build_ext --inplace
pip3 uninstall -y rego
cd /tmp/build
python3 -m pip uninstall rego -y
python3 -m pip --no-cache-dir install python/




