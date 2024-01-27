#include "eval.h"
#include "matrix.h"
#include <stddef.h>

#ifndef _INTRPL_H_
#define _INTRPL_H_

#ifdef __cplusplus
extern "C" {
#endif

void intrpl_parametrize(vec3_t* p, int n, double* u);
int  intrpl_ppoly1(vec3_t* p, int n, ppoly_t* ppoly);
int  intrpl_ppoly3(vec3_t* p, int n, vec3_t* v0, vec3_t* v1, luctx_t* ctx, ppoly_t* ppoly);
int  intrpl_ppoly5(vec3_t* p, int n, vec3_t* v0, vec3_t* v1, vec3_t* a0, vec3_t* a1, luctx_t* ctx, ppoly_t* ppoly);
int  intrpl_arc(vec3_t v0, vec3_t v1, vec3_t c, double ao, arc_t* arc);

#ifdef __cplusplus
}
#endif

#endif
