# rego: Automatic Time Series Forecasting and Missing Value Imputation
#
# Copyright (C) Davide Altomare and David Loris <channelattribution.io>
# 
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree. 

from setuptools import setup
from Cython.Distutils import build_ext, Extension

# from distutils.core import setup
# from distutils.extension import Extension
#from Cython.Distutils import build_ext 

extensions = [Extension(name="rego", 
                        sources=["src/cypack/rego.pyx", "src/cypack/functions.cpp"],
						include_dirs=["src/cypack/armadillo/include"], 
						language='c++', 
						extra_compile_args=['-std=c++11'],
						extra_link_args=['-pthread','-llapack', '-lblas'],
						cython_directives={"language_level":'3',"embedsignature": True})
]
 
setup(
	ext_modules = extensions,    
	cmdclass={'build_ext': build_ext},
	install_requires=['numpy', 'pandas']
)


