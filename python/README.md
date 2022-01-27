*rego*
======

*Automatic Time Series Forecasting and Missing Value Imputation*

*rego* is a machine learning algorithm for predicting and imputing time series. It can automatically set all the parameters needed, thus in the minimal configuration it only requires the target variable and the dependent variables if present. It can address large problems with hundreds or thousands of dependent variables and problems in which the number of dependent variables is greater than the number of observations. Moreover it can be used not only for time series but also for any other real valued target variable. The algorithm implemented includes a Bayesian stochastic search methodology for model selection and a robust estimation based on bootstrapping. *rego* is fast because all the code is C++.

PyPi installation
-----------------

Note! Only Python3 is supported!

```bash
pip install --upgrade setuptools
pip install wheel
pip install Cython
pip install pandas
pip install rego
```


### Compile from source

```bash
cd /python

python -m venv env
source env/bin/activate

pip install -r requirements.txt
python setup.py build_ext --inplace
