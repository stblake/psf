# psf
Perfect Sequence Finder (psf) - finds constructions for perfect sequences and arrays. 

This program searches for n-phase two dimensional arrays with perfect periodic autocorrelation. The index function of the array is of the form 

S(i,j) = \omega^{\sum_{k=1}^{N}floor(p_k(x,y)/d_k)}, 

where \omega = e^(2 \pi \sqrt{-1}/n), p_k(x,y) are bivariate polyomials and d_k \elem Z_n. Per Ma, Ng "On non-existence of perfect and nearly perfect sequences", Int. J. Information and Coding Theory, Vol. 1, No. 1, 2009, n should not be prime. The polynomials, p_k(x,y), and denominators are randomly generated within user-specified constraints. 

The program requires you to specify the following: 

* The minimum and maximum number of phases. For example: phases = [6,8,12] 

* Bounds on the number of elements in the arrays. For example: size_ranges = [[36,1296]] will generate arrays with all sizes between 36 and 1296. 

* The modulus of the d_k's and p_k(x,y)'s. By default the modulus is the phase of the array. For example: modulus = 'Automatic' will use the array phase. 

* The denominators, by default are generated randomly in Z_n. However they can be user-specified. For example: denominators = [2,3,5]. 

* The number of trials for each configuration of phase, array size, and denominator. For example, n_trials = 100. 

* The number of sums, N. For example n_sums = 2.

* The maximum degree of the polynomial p_k(x,y) in x. For example degX = 2. 

* The maximum degree of the polynomial p_k(x,y) in y. For example degY = 2. 

* The tolerance for a perfect array. For example: tol_perfect = 0.01 specifies that the off-peak shifts be 1% or less of the peak. 

* A flag to specify that the array is square. square_only = True. By default square_only = False. 

* A flag to specify that the dimensions of the arrays are coprime only. diagonal_only = True. By default diagonal_only = False. 

If a perfect array is found the array is then checked for the Array Orthogonality Property. If so, a perfect sequence has been discovered. If the array has coprime dimensions then a perfect sequence has been discovered. 

Numerous examples are given in searches.py. 

Running this program with efficient correlations uses FFTW (http://www.fftw.org/). Instructions for installing FFTW using brew is given at https://brewformulas.org/Fftw and instructions for the installation of pyfftw is given at https://github.com/pyFFTW/pyFFTW. 

Alternatively you can use FFTpack (https://docs.scipy.org/doc/scipy/reference/tutorial/fftpack.html) by replacing crosscorrelate_fftw with crosscorrelate_fftpack in max_abs_off_peak(). 

In the near future I'll test out replacing the C code for fill_array_2d with a numba compiled python equivalent. Ideally this should not slow down the program too much, if at all, while making everything significantly simpler. 




