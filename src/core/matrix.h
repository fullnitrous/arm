#include <math.h>

#ifndef _MATRIX_H_
#define _MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

#define APPROACHES_ZERO(x) (fabs(x) < 1e-6)

typedef struct lu_context {
	int n, *p;
	double *m, *b;
} luctx_t;

int mat_lud(luctx_t* ctx);
void mat_lusolve(luctx_t* ctx);

#ifdef __cplusplus
}
#endif

#endif
