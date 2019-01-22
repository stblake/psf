

#include <math.h>
#include <stdio.h>

void fill_array_2d(float *array_re, float *array_im, 
		int d0, int d1, 
		int *cs, int n_cs, 
		int *ns, int n_ns, 
		int *polyxy, int degx, int degy, 
		int n_sums, int phase);
long int abstract_index_fn(int i, int j, int n_sums, int *cs, int n_cs, 
	int *ns, int n_ns, int *polyxy, int degx, int degy, int phase);
long int poly_eval_2d(int *polyxy, int degx, int degy, int y, int x, int k, int phase);
long int power(int base, int exponent);