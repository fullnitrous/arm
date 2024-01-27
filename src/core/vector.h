#include <math.h>
#include <stddef.h>
#include <stdint.h>

#ifndef _VECTOR_H_
#define _VECTOR_H_

typedef struct vec2 { double x, y;       } vec2_t;
typedef struct vec3 { double x, y, z;    } vec3_t;
typedef struct vec4 { double x, y, z, w; } vec4_t;

extern const int __vec3_member_offsets[3];

#define VEC3GET(v, i) ((double*)((uint8_t*)(v) + __vec3_member_offsets[i]))

#ifdef __cplusplus
extern "C" {
#endif

vec3_t vecadd(vec3_t u, vec3_t v);
vec3_t vecsub(vec3_t u, vec3_t v);
vec3_t vecmul(vec3_t u, double s);
vec3_t vecdiv(vec3_t u, double s);
vec3_t vecabs(vec3_t u);
vec3_t vecrotx(vec3_t u, double a);
vec3_t vecroty(vec3_t u, double a);
vec3_t vecrotz(vec3_t u, double a);
vec3_t veccross(vec3_t u, vec3_t v);
double vecdot(vec3_t u, vec3_t v);
double vecmag(vec3_t u);
vec3_t vecproj(vec3_t u, vec3_t v);
vec3_t vecdir(vec3_t u, vec3_t v);
vec3_t vecnorm(vec3_t u);
double vecangle(vec3_t u, vec3_t v);
vec3_t veczeros(void);
vec3_t vecrotn(vec3_t v, vec3_t n, double a);
int    veccmp(vec3_t tst, vec3_t ref, double tol);

#ifdef __cplusplus
}
#endif

#endif
