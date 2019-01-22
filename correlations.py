#
#    PERIODIC CORRELATIONS VIA FFT 
#

# Written by Sam Blake.

# TODO: use FFTW via pyFFTW

import math
import numpy as np
from scipy.fftpack import fftn, ifftn
import pyfftw 

# See https://hgomersall.github.io/pyFFTW/pyfftw/interfaces/scipy_fftpack.html#module-pyfftw.interfaces.scipy_fftpack

def crosscorrelate_fftw(a,b):
	"As with crosscorrelate_fftpack, but much much faster due to fftw."
	# TODO: threads should not be manually set here. 
	# TODO: use real ffts.
	a_fft = pyfftw.interfaces.scipy_fftpack.fftn(a, threads = 1, planner_effort = 'FFTW_ESTIMATE')
	b_fft = np.conj(pyfftw.interfaces.scipy_fftpack.fftn(b, threads = 1, planner_effort = 'FFTW_ESTIMATE')) 
	ab_fft = np.multiply(a_fft, b_fft)
	ab = pyfftw.interfaces.scipy_fftpack.ifftn(ab_fft, threads = 1, planner_effort = 'FFTW_ESTIMATE') 

	#if np.allclose(np.imag(ab), 0.0, 1.0e-5):
	#	return np.real(ab)
	#else:
	#	return ab

	return ab


def crosscorrelate_fftpack(a,b):
	"Computes the periodic cross-correlation of the N-dimensional arrays /a/ and /b/ using the \
FFT and IFFT. This is the slow version which uses scipy.ffpack"

	a_fft = fftn(a)
	b_fft = np.conj(fftn(b))
	ab_fft = np.multiply(a_fft, b_fft)
	ab = ifftn(ab_fft)

	if np.allclose(np.imag(ab), 0.0, 1.0e-5):
		return np.real(ab)
	else:
		return ab

def autocorrelate_fftw(a): 
	"Computes the periodic autocorrelation from the cross-correlation."
	
	auto = crosscorrelate_fftw(a,a)

	return auto

def autocorrelate_fftpack(a): 
	"Computes the periodic autocorrelation from the cross-correlation."
	return crosscorrelate_fftpack(a,a)


def max_off_peak(a):
	"Computes the absolute maximum off-peak correlation."
	na = np.abs(a)
	na.flat[0] = 0.0
	return np.max(na)





