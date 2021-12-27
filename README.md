<div align="center">
    <img width="50%" src="img/rego.png?">
</div>

<br />
*rego* is a machine learning algorithm for predicting and imputing time series. It can automatically set all the parameters needed, thus in the minimal configuration it only requires the target variable and the regressors if present. It can address large problems with hundreds or thousands of regressors and problems in which the number of regressors is greater than the number of observations. Moreover it can be used not only with time series but also with any other real valued target variable. The algorithm implemented includes a bayesian stochastic search methodology for model selection and a robust estimation based on boostrapping. *rego* is fast because all the code is C++. 

Installation notes
-------------------

Compilation requires a C++11 compiler, lapack and blas installed. For a debian based SO, they can be easily installed running:

```bash
apt-get install build-essential liblapack-dev libblas-dev
```

Installation on Windows requires Microsoft Visual C++ 14.0 or greater. (https://visualstudio.microsoft.com/it/downloads/)


Compiling and testing C++ code
-------------------------

```bash
sh compile/compile-cpp.sh test/test.cpp
/tmp/build/c++/test.o
```

Compiling and testing Python code
----------------------------

Note! Only Python3 is supported! 

```bash
pip install -r python/requirements.txt


sh compile/compile-py.sh
python test/test.py 
```

Compiling and testing R code
-----------------------

```bash
sh compile/compile-R.sh
Rscript BATCH test/test.R
```

Installation from PyPi
----------------------

```bash
pip install --upgrade setuptools
pip install Cython
pip install pandas
pip install rego
```

Installation from CRAN
----------------------

```bash
R --vanilla -e 'install.packages(c("Rcpp", "RcppArmadillo"), repos="http://cran.us.r-project.org")'
R --vanilla -e 'install.packages(c("rego"), repos="http://cran.us.r-project.org")'
```

Generating Python documentation
-------------------------------

```bash
pip install sphinx
pip install numpydoc
pip install pygments --upgrade
pip install rinotype

cd python/src/cypack

python generate_doc.py
```

The following .pdf will be generated:

```bash
python/src/cypack/docs/_build/rinoh/rego.pdf
```
