#include "vector.h"
#include "matrix.h"
#include <stdint.h>
#include <assert.h>

#ifndef _EVAL_H_
#define _EVAL_H_

#define EVAL_BEZIER_MAX_DEG  10
#define EVAL_BSPLINE_MAX_DEG 10
#define EVAL_NURBS_MAX_DEG   10

#define EVAL_PPOLY1_MEMSIZE(n) (2*n)

#define EVAL_PPOLY3_MEMSIZE(n)       (4*n)
#define EVAL_PPOLY3_EQSYS_MEMSIZE(n) (4*(n-1)*4*(n-1))
#define EVAL_PPOLY3_PIVOT_MEMSIZE(n) (4*(n-1))
#define EVAL_PPOLY3_B_MEMSIZE(n)     (4*(n-1))

#define EVAL_PPOLY5_MEMSIZE(n)       (6*n)
#define EVAL_PPOLY5_EQSYS_MEMSIZE(n) (6*(n-1)*6*(n-1))
#define EVAL_PPOLY5_PIVOT_MEMSIZE(n) (6*(n-1))
#define EVAL_PPOLY5_B_MEMSIZE(n)     (6*(n-1))

#define EVAL_BEZIER_MEMSIZE(n) (n+1)

#ifdef __cplusplus
extern "C" {
#endif

extern int eval_bico_lut[EVAL_BEZIER_MAX_DEG+1][EVAL_BEZIER_MAX_DEG+1];

typedef enum e3d_types {
	E3D_PPOLY1,
	E3D_PPOLY3,
	E3D_PPOLY5,
	E3D_ARC,
	E3D_BEZIER,
	E3D_BSPLINE,
	E3D_NURBS
} e3dt_t;

typedef struct evaluation {
	vec3_t p, v, a;
	double m;
} eval_t;

typedef struct piecewise_polynomial {
	vec3_t*  c;
	double*  u;
	int      n;
} ppoly_t;

typedef struct arc {
	vec3_t n, c, r;
	double a;
} arc_t;

typedef struct bezier {
	vec3_t* p;
	int     n;
} bezier_t;

typedef struct bspline {
	vec3_t* p;
	double* u;
	int d, n;
} bspline_t;

typedef struct nurbs {
	vec4_t* p;
	double* u;
	int d, n;
} nurbs_t;

typedef struct evaluator_3d {
	eval_t (*evaluator)(void*, double);
	eval_t (*__evaluator)(void*, int, double);
	vec3_t (*__evaluator_fast)(void*, int, double);
	int k;
	void* function;
	e3dt_t type;
} e3d_t;

int eval_findinterval(double* u, double t, int n);
int eval_findspan(double* u, double t, int n, int p);

void eval_bezier_lut_init(void);
int __eval_fac(int n);
int __eval_bico(int n, int i);
int eval_bico(int n, int i);
double eval_bernstein_basis(int i, int n, double t);

ppoly_t   eval_ppoly_init(vec3_t* mem0, double* mem1, int n);
luctx_t   eval_ppoly_luctx_init(int* pivotmem, double* matrixmem, double* bmem);
bezier_t  eval_bezier_init(vec3_t* mem, int n);
bspline_t eval_bspline_init(vec3_t* mem0, double* mem1, int d, int n);
nurbs_t   eval_nurbs_init(vec4_t* mem0, double* mem1, int d, int n);

eval_t eval_ppoly1(void* f_, double t);
eval_t __eval_ppoly1(void* f_, int k, double t);
vec3_t __eval_ppoly1_fast(void* f_, int k, double t);

eval_t eval_ppoly3(void* f_, double t);
eval_t __eval_ppoly3(void* f_, int k, double t);
vec3_t __eval_ppoly3_fast(void* f_, int k, double t);

eval_t eval_ppoly5(void* f_, double t);
eval_t __eval_ppoly5(void* f_, int k, double t);
vec3_t __eval_ppoly5_fast(void* f_, int k, double t);

eval_t eval_arc(void* f_, double t);
vec3_t __eval_arc_fast(void* f_, int k, double t);

eval_t eval_bezier(void* f_, double t);
vec3_t __eval_bezier_fast(void* f_, int k, double t);

eval_t eval_bspline(void* f_, double t);
eval_t __eval_bspline(void* f_, int k, double t);
vec3_t __eval_bspline_fast(void* f_, int k, double t);

eval_t eval_nurbs(void* f_, double t);
eval_t __eval_nurbs(void* f_, int k, double t);
vec3_t __eval_nurbs_fast(void* f_, int k, double t);

int    eval_cmp(eval_t tst, eval_t ref, double tol);
eval_t eval_cfd(e3d_t* e3d, double h, double t);

double eval_maxv(eval_t* e, double maxat);

e3d_t  eval_e3d_make(void* function, e3dt_t type);

eval_t eval_e3d(e3d_t* e3d, double t);
eval_t __eval_e3d(e3d_t* e3d, double t);
vec3_t __eval_e3d_fast(e3d_t* e3d, double t);

int    eval_e3dinterval(e3d_t* e3d, double t);
double eval_e3dbound(e3d_t* e3d, int k);
int    eval_e3dnintervals(e3d_t* e3d);

#ifdef __cplusplus
}
#endif

#endif
