Repository hosting work in progress towards the construction of a Python package for Computation of Lagrangian Descriptors for Nonlinear Dynamical Systems.


# CHANGE LOG

* Redefined inputs/outputs as `(t, x)`, where `x` has to be a 1D-array of length `(2*Nx*Ny)` for every vector field function. This is needed to be able to use `scipy.integrate.solv_ivp`. Internally every vector field function, need to reshape the array to `(Nx*Ny, 2)` for matrix operations though.

* Temporarily, `forcing` was removed from `HamiltonCentre` and `HamiltonSaddle` vector field functions, but internally incorporated within `HamiltonDuffing`.

* Removed redundant call to perturbation parameters `pert_params` by LD calculator.

* Swapped variable and function names for more self-descriptive ones.



# TO-DOs

* Continue renaming variables and functions to be self-descriptive
* Remove redundant variable calls
* Reduce number of input arguments for all functions, use lists of variables or dictionaries to pass parameters to functions, for example. `PARAMETERS`
* Add description within functions
* Continue commenting functions


## NEXT STAGES
* Break down functions to achieve a more modular structure
* Maybe define classes of objects
