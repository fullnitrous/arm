#include <math.h>

#ifndef _MISC_H_
#define _MISC_H_

#ifdef M_PI
#undef M_PI
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define M_PI 3.14159265358979323846264338327
#define IVAL(a, x, b) (a <= x && x <= b)
#define FEQ(a, b, tol) (-tol <= b-a && b-a <= tol)
#define FZ(a, tol) (-tol <= a && a <= tol)

double range(double start, double stop, int i, int res);
double min(double a, double b);
double max(double a, double b);

#ifdef __cplusplus
}
#endif

#endif
