# rego: Automatic Time Series Forecasting and Missing Value Imputation
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 


import os

dir=os.path.dirname(__file__)

TAG = '#start documentation'

tag_found = False
with open(os.path.join(dir, 'rego.pyx')) as in_file:
    with open(os.path.join(dir,'docs/rego.py'), 'w') as out_file:
        for line in in_file:
            if not tag_found:
                if line.strip() == TAG:
                    tag_found = True
            else:
                out_file.write(line)
				

if not os.path.exists(os.path.join(dir,"docs/")):
    os.makedirs(os.path.join(dir,"docs/"))				

os. chdir(os.path.join(dir,"docs/"))
os.system("make clean && make html")
os.system("sphinx-build -b rinoh . _build/rinoh")

