#ifndef _ALP_H_
#define _ALP_H_

#define ALP_PERF

#include "vector.h"
#include "eval.h"
#include "approx.h"
#include <math.h>

#define ALP_SOLVE_TOL     1e-10
#define ALP_SOLVE_GUESS_N 5
#define ALP_SOLVE_NMAX    50
#define ALP_SOLVE_NR_N    10

#define ALP_CHK          30
#define ALP_CHECKRES     150
#define ALP_GUESS_FAC    4.0
#define ALP_MIN_INTERVAL 1e-7

#ifndef ALP_PERF
#define ALP_MEMSIZE(n) (n+1+ALP_CHK*n)
#else
#define ALP_MEMSIZE(n) (n+1+2*ALP_CHK*n)
#endif

extern double (*alp_integrator)(APPROX_FUNC1, void*, double, double);

typedef struct function_context {
    double a, b, l;
    e3d_t* f;
    solvecfg_t* cfg;
} alpfctx_t;

typedef struct arc_length_parameterization {
	double *u;

#ifdef ALP_PERF
	double* a_der;
	e3d_t* f;
#endif
	union {
		double* a;
		double  s;
	} c;
	int n;
} alp_t;

#ifdef ALP_PERF
double alp_eval_v(alp_t* alp, double l);
#endif

double alp_eval(alp_t* alp, double l);
alp_t alp_init(double* mem, int n_upper);
double alp_fit(e3d_t* f, alp_t* alp, int n_upper, double tolerance);

#endif
