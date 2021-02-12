#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:36:43 2020

@author: vk17590
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from pylds.base import fit_pes

def discretise_potential(coords, potential):
    """
    Returns a 1- or 2-dimensional array of function (potential energy surface) values on a grid of points.

    Parameters
    ----------
    coords : list of ndarrays
        [x] or [x,y] contain coordinates.

    potential : function
        Potential to be discretised.

    Returns
    -------
    pes_data : ndarray, shape(len(x)) or shape(len(y),len(x))
        pes_data is an array of potential values.
    """

    if len(coords) == 1:
        pes_data = potential(coords)
        return pes_data

    elif len(coords) == 2:
        x, y = coords
        X,Y = np.meshgrid(x,y,indexing='xy')
        points = np.column_stack((X.flatten(),Y.flatten()))
        pes_data = potential(points)
        return pes_data.reshape(len(x),len(y))

    else:
        print('splines in +3D are not implemented yet')

def generate_pes_data(coords, potential, filename):
    """
    Saves 1- or 2-dimensional array of function (potential energy surface) values on a grid of points to pylds/pes_files/filename.hdf5. File format fixed to HDF5 by default.

    Parameters
    ----------
    coords : list (or ndarray) of ndarrays
        [x] or [x,y] contain coordinates.

    potential : function
        Potential to be discretised.

	filename : string
		Name of hdf5 file to be saved.
    """

    pes_data = discretise_potential(coords, potential)
    
    dirname = 'pylds/pes_files'
    filepath = os.path.join(dirname, filename+'.hdf5')
    hf = h5py.File(filepath,'w')
    hf.create_dataset('coords', data=np.array(coords).astype('float64'))
    hf.create_dataset('pes_data', data=pes_data.astype('float64'))
    hf.close()
