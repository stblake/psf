

import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer


lib_fill_array_2d = ctypes.cdll.LoadLibrary('/home/sblake/pypsf/_lib_fill_array_2d.so')
_fill_array_2d = lib_fill_array_2d.fill_array_2d


# void fill_array_2d(float *array_re, float *array_im, int d0, int d1, int *cs, int n_cs, 
#						int *ns, int n_ns, int *polyxy, int degx, int degy, int n_sums, int phase);

_fill_array_2d.argtypes = ( 
	ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), 
	ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), 
	ctypes.c_int, ctypes.c_int, 
	ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int, 
	ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
	ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_int,
	ctypes.c_int, ctypes.c_int)

_fill_array_2d.restype = None


def fill_array_2d(array, n_sums, Cs, Ns, pxy, phase):

	nrows, ncols = array.shape
	degx, degy = pxy.shape[1] - 1, pxy.shape[2] - 1

	array_re = np.zeros(array.shape, dtype = np.float32)
	array_im = np.zeros(array.shape, dtype = np.float32)

# void fill_array_2d(float *array_re, float *array_im, int d0, int d1, int *cs, int n_cs, 
#	int *ns, int n_ns, int *polyxy, int degx, int degy, int n_sums, int phase);

	_fill_array_2d(array_re, array_im, nrows, ncols, Cs, Cs.shape[0], Ns, Ns.shape[0], pxy, degx, degy, n_sums, phase)	

	array[:,:] = array_re[:,:] + 1j*array_im[:,:]

	return 
