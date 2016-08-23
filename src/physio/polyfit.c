#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polyfit.h"

void solve(double *, float *, int);

/* polyfitlin:							*/
/* Arguments -							*/
/*   data   : an array of floats containing data points         */
/*   n      : # of data points					*/
/*   beta   : array of floats for holding the coefficients 	*/
/*   order  : order of the polynomial to fit to 		*/
/* Returns -							*/
/*   mean squared error of the data with the fitted polynomial  */

float polyfitlin(float *data, int n, float *beta, int order) {
	int *x, i;
	float result;

	x = (int *) malloc(n * sizeof(int));
	for(i = 0; i < n; i++) 
		x[i] = i;

	result = polyfit(data, x, n, beta, order);

	free(x);

	return result;
}

/* polyfit:							*/
/* Arguments - 							*/
/*   data   : an array of float containing data points 		*/
/*   x      : an array of x values for the data points          */
/*   n      : # of data points					*/
/*   beta   : array of floats for holding the cofficients	*/
/*   order  : order of the polynomial to fit to			*/
/* Returns -							*/
/*   mean squared error of the data with the fitted polynomial	*/

float polyfit(float *data, int *x, int n, float *beta, int order) {
	int i, j, k, k2, k3;
	double sum, xi, val, error;

	double *basis = (double *) malloc(order * n * sizeof(double));
	double *alpha = (double *) malloc(order * order * sizeof(double));

	/* Calculate basis */

	for(i = 0; i < order; i++) 
		for(j = 0; j < n; j++) {
			k = i + j * order;
			if(i == 0) basis[k] = 1;
			else basis[k] = x[j] * basis[k - 1];
		}

	/* Calculate alpha */

	for(i = 0; i < order; i++) {
		for(j = 0; j <= i; j++) {
			sum = 0;
			for(k = 0; k < n; k++) {
				k2 = i + k * order;
				k3 = j + k * order;

				sum += basis[k2] * basis[k3];
			}

			k2 = i + j * order;
			k3 = j + i * order;
			alpha[k2] = sum;
			alpha[k3] = sum;

		}
	}
	

	/* Calculate beta */

	for(i = 0; i < order; i++) {
		sum = 0;
		for(j = 0; j < n; j++) {
			k2 = i + j * order;
			sum += data[j] * basis[k2];
		}
		beta[i] = sum;
	}

	/* Solve for everything */

	solve(alpha, beta, order);

	/* Now calculate and return the mean squared error */

	sum = 0;
	for(i = 0; i < n; i++) {
		for(j = 0, xi = 1, val = 0; j < order; j++, xi *= x[i])
			val += beta[j] * xi;
		error = val - data[i];
		sum += error * error;
	}

	free(basis);
	free(alpha);

	return sum/n;
	

}


/* solve: 							*/
/* Uses Gaussian elimination to finish up the job for polyfit 	*/

void solve(double a[], float b[], int n) {
	double mag, mag2, temp;
	int pivot;
	int i, j, k;

	for(i = 0; i < n; i++) {
		mag = 0;
		pivot = -1;

		for(j = i; j < n; j++) {
			mag2 = fabs(a[i + j * n]);
			if(mag2 > mag) {
				mag = mag2;
				pivot = j;
			}
		}

		if(pivot == -1 || mag == 0) return;

		if(pivot != i) {
			for(j = i; j < n; j++) {
				temp = a[j + i * n];
				a[j + i * n] = a[j + pivot * n];
				a[j + pivot * n] = temp;
			}

			temp = b[i];
			b[i] = b[pivot];
			b[pivot] = temp;
		}

		mag = a[i + i * n];
		for(j = i; j < n; j++) 
			a[j + i * n] /= mag;
		b[i] /= mag;

		for(j = 0; j < n; j++) {
			if(j == i) continue;
		
			mag2 = a[i + j * n];

			for(k = i; k < n; k++)
				a[k + j * n] -= mag2 * a[k + i * n];
		
			b[j] -= mag2 * b[i];
		}

	}
}
