Repository hosting work in progress towards the construction of a Python package for Computation of Lagrangian Descriptors for Nonlinear Dynamical Systems.


# Install minimal package dependencies before running Notebooks


First, make sure to use `pip` within the `conda` environment.

So, do as shown below.


```bash
# FIRST
conda install pip
# SECOND
pip install -r requeriments.txt
```


<span style='color:red'><b>NOTE</b></span> It is important to keep in mind that IN GENERAL portability across different _Operating Systems_ is not guaranteed. This means that exact package versions for some libraries may not be able to be installed in both Windows or Linux, for instance. This has some technical reasons. Also, installing identical libraries may find troubles for different CPU architectures. 


# CHANGE LOG

* Incorporated suggested changes for code transparency (Thanks v!) 
* Added `requirements.txt`


# TO-DOs

* Add `docstrings` for description within functions
* Continue commenting functions

# NEXT STAGES
* Break down functions to achieve a more modular structure
* Maybe define classes of objects
