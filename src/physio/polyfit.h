/* Polyfit:							*/
/* Arguments - 							*/
/*   data   : an array of doubles containing data points	*/
/*   n      : # of data points					*/
/*   beta   : array of doubles for holding the cofficients	*/
/*   order  : order of the polynomial to fit to			*/
/* Returns -							*/
/*   mean squared error of the data with the fitted polynomial	*/

float polyfit(float *, int *, int, float *, int);
float polyfitlin(float *, int, float *, int);

