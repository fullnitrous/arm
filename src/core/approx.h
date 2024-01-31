#ifndef _APPROX_H_
#define _APPROX_H_

#include "eval.h"
#include "util.h"
#include "matrix.h"

#define APPROX_CHEBYSHEV_MAX_K               30
#define APPROX_CHEBYSHEV_MEMSIZE(n)          (n)
#define APPROX_CHEBYSHEV_MEMSIZE_WITH_INT(n) (2*n)

#define APPROX_FUNC1 double (*f)(void*, double)
#define APPROX_FUNC2 vec2_t (*f)(void*, double)

typedef struct chebyshev {
	double a, b, *c, *cint;
	int n;
} chebyshev_t;

typedef struct approx_solve_config {
	double a, b, tol;
	int bi_guess_n, bi_nmax, only_bisect, nr_n;
} solvecfg_t;

chebyshev_t approx_chebyshev_init(double* m, double a, double b, int n);
void   approx_chebyshev(APPROX_FUNC1, void* f_, chebyshev_t* apx);
void   approx_chebyshev_int(chebyshev_t* apx);
double approx_eval_chebyshev(chebyshev_t* apx, int der, double x);

double approx_glq32(APPROX_FUNC1, void* f_, double a, double b);
double approx_solve(APPROX_FUNC2, void* f_, solvecfg_t* cfg, double y);


#endif
