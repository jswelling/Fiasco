extern void cdf_bet(int* which ,double* p, double* q, double* x, double* y, 
		    double* a, double* b, int* status, double* bound);
extern void cdf_bin(int* which, double* p, double* q, double* s, double* xn, 
		    double* pr, double* ompr, int* status, double* bound);
extern void cdf_chi(int* which, double* p, double* q, double* x, double* df, 
		    int* status, double* bound);
extern void cdf_chn(int* which, double* p, double* q, double* x, double* df,
		    double* pnonc, int* status, double* bound);
extern void cdf_f(int* which, double* p, double* q, double* f, double* dfn,
		  double* dfd, int* status, double* bound);
extern void cdf_fnc(int* which, double* p, double* q, double* f, double* dfn,
		    double* dfd, double* phonc, int* status, double* bound);
extern void cdf_gam(int* which, double* p, double* q, double* x, 
		    double* shape, double* scale, int* status, double* bound);
extern void cdf_nbn(int* which, double* p, double* q, double* s, double* xn,
		    double* pr, double* ompr, int* status, double* bound);
extern void cdf_nor(int* which, double* p, double* q, double* x, 
		    double* mean, double* sd, int* status, double* bound);
extern void cdf_poi(int* which, double* p, double* q, double* s, double* xlam,
		    int* status, double* bound);
extern void cdf_t(int* which, double* p, double* q, double* t, double* df,
		  int* status, double* bound);
extern int cdf_ipmpar(int* which);

