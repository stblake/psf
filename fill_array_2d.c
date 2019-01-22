




#include "fill_array_2d.h"

#define PI 3.14159265358979

void fill_array_2d(float *array_re, float *array_im, 
		int d0, int d1, 
		int *cs, int n_cs, 
		int *ns, int n_ns, 
		int *polyxy, int degx, int degy, 
		int n_sums, int phase) {

	int i, j, k;
	long int indx;
	float c;

	// printf("\ndegx,degy = %d,%d\n", degx,degy);

	c = 2.0*PI/((float) phase);

	k = 0;

	for (i = 0; i < d0; i++) {
		for (j = 0; j < d1; j++) {
			
			indx = abstract_index_fn(j, i, n_sums, cs, n_cs, ns, n_ns, polyxy, degx, degy, phase);
			indx = indx%phase;

			// e^(i t) = cos(t) + I sin(t), t = 2 pi indx / phase
			array_re[k] = cos(c*((float) indx));
			array_im[k] = sin(c*((float) indx));
			k++;
		}
	}
}



long int abstract_index_fn(int i, int j, int n_sums, int *cs, int n_cs, int *ns, int n_ns, int *polyxy, int degx, int degy, int phase) {

	int k;
	long int indx, flr;

	indx = 0;

	for (k = 0; k < n_sums; k++) {
		flr = poly_eval_2d(polyxy, degx, degy, i, j, k, phase)/ns[k];
		indx += cs[k]*flr;
		// printf("\ncs,ns,indx,flr = %d,%d,%d,%d\n", cs[k], ns[k], indx, flr);
	}

	return indx;
}



long int poly_eval_2d(int *polyxy, int degx, int degy, int x, int y, int k, int phase) {

	int i, j, n;

	long int pev = 0;
	n = 0;

	for (i = 0; i <= degx; i++) {
		for (j = 0; j <= degy; j++) {
			pev += polyxy[k*(degx + 1)*(degy + 1) + n++]*power(x,i)*power(y,j);
			// printf("\nx,i,y,j,coeff,pev = %d, %d, %d, %d, %d, %d", x,i, y,j, polyxy[k*(degx + 1)*(degy + 1) + n - 1], pev);
		}
	}

	return pev;
}


long int power(int base, int exponent) {

  if (exponent == 0) {
    return 1;
  } else if (exponent % 2) {
    return base * power(base, exponent - 1);
  } else {
    long int temp = power(base, exponent / 2);
    return temp * temp;
  }
}




