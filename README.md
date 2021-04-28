**LDDS**: Python package for Computation of Lagrangian Descriptors of Dynamical Systems.

## Table of contents
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


## Description



## Dependencies and installation

The `setup.py` should install the dependencies listed in
[requirements.txt](https://github.com/champsproject/ldds/blob/develop/requirements.txt) using

``` bash
> pip install -r requirements.txt (or pip3 install -r requirements.txt)
```


### Installing from source
Clone the git repository and install `ldds` as a module using

``` bash
> git clone git@github.com:champsproject/ldds.git
> cd ldds
> python setup.py install
```

Now `ldds` module is available for import along with the functions to run examples and tests. The following lines of code imports the `ldds` module, imports the function to compute Lagrangian descriptor, imports the vector field and Hamiltonian for the two degrees of freedom quadratic normal form Hamiltonian with index one saddle.

``` python
import ldds
from ldds.base import compute_lagrangian_descriptor
from ldds.vector_fields import quadratic_normalform_saddlecenter
from ldds.hamiltonians import quadratic_normal_form_saddlecenter_ha
```

## Documentation

**LDDS** uses [Sphinx](http://www.sphinx-doc.org/en/stable/) for documentation and is made available online [here](https://ldds.readthedocs.io/en/latest/?badge=latest#) using [Read the Docs](https://readthedocs.org/). To build the html version of the docs locally simply:

```bash
> cd docs
> make html
```

The generated html can be viewed by opening `docs/_build/html/index.html`.


## Contributing 

Guidelines on how to contribute to this package can be found [here](https://github.com/champsproject/ldds/blob/develop/contributing.md) along with the code of conduct [here](https://github.com/champsproject/ldds/blob/develop/code_of_conduct.md) for engaging with the fellow contributors. As and when we receive improvements to the package, we will acknowledge the pull request and the contributor in this section.


## Acknowledgements

We acknowledge the support of EPSRC Grant No. EP/P021123/1 [CHAMPS project](https://champsproject.com). 
