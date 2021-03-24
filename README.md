Repository for the Python package for Computation of Lagrangian Descriptors of Dynamical Systems.

  -   [Summary](#summary)
  -   [Installation](#installation)
  -   [Usage](#usage)
  -   [Tests](#tests)
  -   [Contributing](#contributing)
  -   [Acknowledgements](#acknowledgements)
  -   [Copyright and License](#copyright-and-license)
  -   [References](#references)

LDDS
====

[![Documentation Status](https://readthedocs.org/projects/ldds/badge/?version=latest)](https://ldds.readthedocs.io/en/latest/?badge=latest)

[![DOI](linktoZenodorepo.svg)](linktoZenodo)


Summary
-------


Installation
------------

# Install minimal package dependencies before running Notebooks

First, make sure to use `pip` within the `conda` environment using

```bash
# FIRST
conda install pip
# SECOND
pip install -r requirements.txt
```

Clone the git repository and install `ldds` as a module using

``` bash
$ git clone git@github.com:champsproject/ldds.git
$ cd ldds
$ python setup.py install
```

The `setup.py` should also install the modules listed in
[requirements.txt](https://github.com/champsproject/ldds/blob/develop/requirements.txt)
and which are typically installed using

``` bash
$ pip install -r requirements.txt (or pip3 install -r requirements.txt)
```

Now `ldds` module is available for import along with the methods and
example systems as submodules, for example (SHOW HOW TO IMPORT HERE)

``` python


```
 

Contributing (THIS SECTION NEEDS COMPLETING)
------------

Guidelines on how to contribute to this package can be found

along with the code of conduct

for engaging with the fellow contributors.



Acknowledgements
----------------

We acknowledge the support of EPSRC Grant No. EP/P021123/1 (CHAMPS project). 
