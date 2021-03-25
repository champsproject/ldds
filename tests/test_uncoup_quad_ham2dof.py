# Test suite for obtaining unstable periodic orbits in the uncoupled 
# quartic Hamiltonian

import numpy as np
from scipy.integrate import solve_ivp

import unittest

import ldds
from ldds.base import compute_lagrangian_descriptor
from ldds.vector_fields import quadratic_normalform_saddlecenter
from ldds.hamiltonians import quadratic_normal_form_saddlecenter_ham

# import uposham.uncoupled_quartic_hamiltonian as uncoupled
# import uposham.differential_correction_uncoupled as diff_corr_unc
# import uposham.turning_pt_uncoupled as turning_pt_uncoupled 
# import uposham.turning_pt_coord_diff_uncoupled as tpcd_uncoupled 

# import uposham.differential_correction as diffcorr

import os
# path_to_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), \
                            # 'data/')

# from scipy.spatial.distance import directed_hausdorff

# def hausd_dist_numeric_analytic(orbit, t, deltaE_val, parameters):
#     """
#     Obtain Hausdorff distance between the unstable periodic orbits obtained
#     using a numerical method and analytical solution for the uncoupled 
#     quartic Hamiltonian
#     """

#     total_energy = deltaE_val + parameters[2]
#     # y, py = uncoupled.upo_analytical(total_energy, t, parameters)

#     numerical_orbit = np.array([orbit[:,1], orbit[:,3]])
#     analytical_orbit = np.array([y, py])

#     hausd_dist = directed_hausdorff(numerical_orbit, analytical_orbit)[0]

#     return hausd_dist

class TestContourMap(unittest.TestCase):
    """ 
    Test Lagrangian descriptor contour map for uncoupled quadratic normal form Hamiltonian with an index one saddle 
    """

    def test_uncoupled_ham2dof(self):
        """ Compare with provided data for the system """
    
        Hamiltonian = quadratic_normal_form_saddlecenter_ham
        vector_field = quadratic_normalform_saddlecenter

        H0 = 0.1

        # Integration parameters
        tau = 10

        # LDp, p-value
        p_value = 0.5

        # Mesh parameters
        x1_min, x1_max = [-1.6, 1.6]
        x2_min, x2_max = [-1.6, 1.6]
        x1_res, x2_res = [60, 60]
        slice_parameters = [[x1_min, x1_max, x1_res],[x2_min, x2_max, x2_res]]

        dims_fixed = [0,1,0,0] # Variable ordering (x1 x2 y1 y2)
        dims_fixed_values = [0] # This can also be an array of values
        dims_slice = [1,0,1,0] # Visualisation slice
        momentum_sign = 1 # Direction of momentum that defines the slice - (1) positive / (-1) negative
        grid_parameters = {
            'slice_parameters' : slice_parameters,
            'dims_slice' : dims_slice,
            'dims_fixed' : dims_fixed,
            'dims_fixed_values' : dims_fixed_values,
            'momentum_sign' : momentum_sign,
            'Hamiltonian': Hamiltonian,
            'energy_level': H0
        }
        

        # Obtain LD from the package
        forward_ld = compute_lagrangian_descriptor(grid_parameters, quadratic_normalform_saddlecenter, tau)

        # Load benchmark data
        forward_ld_benchmark = np.loadtxt('./benchmark_data/test.hdf5')


        np.testing.assert_array_almost_equal(forward_ld, 
                                            forward_ld_benchmark, decimal = 3)



if __name__ == "__main__":
    unittest.main()








