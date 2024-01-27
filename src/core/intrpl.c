#include <stdio.h>
#include "intrpl.h"

void intrpl_parameterize(vec3_t* p, int n, double* u) {
	double sum;
	int i;

	sum = 0.0;

	u[0]   = 0;
	u[n-1] = 1;
	
	for(i = 1; i < n; i++) {
		sum += vecmag(vecsub(p[i], p[i-1]));
	}
	
	for(i = 1; i < n-1; i++) {
		u[i] = vecmag(vecsub(p[i], p[i-1])) / sum + u[i-1];
	}
	
	return;
}

int intrpl_ppoly1(vec3_t* p, int n, ppoly_t* ppoly) {
	vec3_t* c;
	double* u;
	int i;

	if(n < 2) { return 0; }
	
	c = ppoly->c;
	u = ppoly->u;

	intrpl_parameterize(p, n, u);

	for(i = 0; i < n-1; i++) {
		c[2*i] = vecdiv(vecsub(p[i+1], p[i]), u[i+1]-u[i]);
		c[2*i+1] = vecsub(p[i], vecmul(c[2*i], u[i]));
	}

	ppoly->n = n;
	
	return 1;
}

int intrpl_spline_eqsys(double* u, int d, int n, double* matrix) {
	int uk, sz, k, i, j, z, q, derivs;
	int deriv, coeff, lhs0, lhs1, lhs, rhs;
	double t, e;

	uk = d + 1;
	sz = uk*(n-1);
	k  = 0;
	e  = 1;
		
	for(i = 0; i < sz*sz; i++) {
		matrix[i] = 0.0;
	}

	for(i = 0; i < n-1; i++) {
		t = u[i];
		e = 1;
		
		lhs0 = 2*i*sz;
		lhs1 = (2*i+1)*sz;

		for(j = uk-1; j >= 0; j--) {
			matrix[lhs0 + k+j] = e;
			e *= t;
		}
	
		t = u[i+1];
		e = 1;

		for(j = uk-1; j >= 0; j--) {
			matrix[lhs1 + k+j] = e;
			e *= t;
		}

		k += uk;
	}

	derivs = d - 1;
	k = 0;

	for(i = 0; i < n-2; i++) {
		t = u[i+1];
		e = 1;

		for(deriv = 0; deriv < derivs; deriv++) {
			z = deriv+1;

			for(j = uk-2-deriv; j >= 0; j--) {
				coeff = z;

				for(q = 0; q < deriv; q++) {
					coeff *= z - q - 1;
				}

				lhs = (2*(n-1)+derivs*i+deriv)*sz;
				rhs = lhs + uk;
				
				matrix[lhs + k+j] =  coeff*e;
				matrix[rhs + k+j] = -coeff*e;

				e *= t;
				z++;

			}
		}

		k += uk;
	}

	return sz;
}

int intrpl_ppoly3(vec3_t* p, int n, vec3_t* v0, vec3_t* v1, luctx_t* ctx, ppoly_t* ppoly) {
	vec3_t* c;
	double* u;
	int i, component, sz;

	if(n < 2) { return 0; }

	c = ppoly->c;
	u = ppoly->u;

	intrpl_parameterize(p, n, u);

	sz = intrpl_spline_eqsys(u, 3, n, ctx->m);
	ctx->n = sz;
	
	if(v0 == NULL) {
		ctx->m[sz*(sz-2) + 0] = 6*u[0];
		ctx->m[sz*(sz-2) + 1] = 2;
	} else {
		ctx->m[sz*(sz-2) + 0] = 3*u[0]*u[0];
		ctx->m[sz*(sz-2) + 1] = 2*u[0];
		ctx->m[sz*(sz-2) + 2] = 1;
	}	
	if(v1 == NULL) {
		ctx->m[sz*(sz-1) + sz-4] = 6*u[n-1];
		ctx->m[sz*(sz-1) + sz-3] = 2;
	} else {
		ctx->m[sz*(sz-1) + sz-4] = 3*u[n-1]*u[n-1];
		ctx->m[sz*(sz-1) + sz-3] = 2*u[n-1];
		ctx->m[sz*(sz-1) + sz-2] = 1;
	}

	if(!mat_lud(ctx)) { return 0; }
	
	for(component = 0; component < 3; component++) {
		for(i = 0; i < sz; i++) { ctx->b[i] = 0.0; }
		for(i = 0; i < n-1; i++) {
			ctx->b[2*i + 0] = *VEC3GET(p+i+0, component);
			ctx->b[2*i + 1] = *VEC3GET(p+i+1, component);
		}
		
		if(v0 != NULL) { ctx->b[sz-2] = *VEC3GET(v0, component); }
		if(v1 != NULL) { ctx->b[sz-1] = *VEC3GET(v1, component); }

		mat_lusolve(ctx);

		for(i = 0; i < sz; i++) {
			*VEC3GET(c+i, component) = ctx->b[i];
		}
	}

	return 1;	
}

int intrpl_ppoly5(vec3_t* p, int n, vec3_t* v0, vec3_t* v1, vec3_t* a0, vec3_t* a1, luctx_t* ctx, ppoly_t* ppoly) {
	vec3_t* c; 
	double* u;
	int i, sz, component;
	
	if(n < 2) { return 0; }

	c = ppoly->c;
	u = ppoly->u;

	intrpl_parameterize(p, n, u);

	sz = intrpl_spline_eqsys(u, 5, n, ctx->m);
	ctx->n = sz;

	if(v0 == NULL) {
		ctx->m[sz*(sz-4) + 0] = 60*u[0]*u[0];
		ctx->m[sz*(sz-4) + 1] = 24*u[0];
		ctx->m[sz*(sz-4) + 2] = 6;
	} else {
		ctx->m[sz*(sz-4) + 0] = 5*u[0]*u[0]*u[0]*u[0];
		ctx->m[sz*(sz-4) + 1] = 4*u[0]*u[0]*u[0];
		ctx->m[sz*(sz-4) + 2] = 3*u[0]*u[0];
		ctx->m[sz*(sz-4) + 3] = 2*u[0];
		ctx->m[sz*(sz-4) + 4] = 1;
	}

	ctx->m[sz*(sz-3) + 0] = 20*u[0]*u[0]*u[0];
	ctx->m[sz*(sz-3) + 1] = 12*u[0]*u[0];
	ctx->m[sz*(sz-3) + 2] = 6*u[0];
	ctx->m[sz*(sz-3) + 3] = 2;

	if(v1 == NULL) {
		ctx->m[sz*(sz-2) + sz-6] = 60*u[n-1]*u[n-1];
		ctx->m[sz*(sz-2) + sz-5] = 24*u[n-1];
		ctx->m[sz*(sz-2) + sz-4] = 6;
	} else {
		ctx->m[sz*(sz-2) + sz-6] = 5*u[n-1]*u[n-1]*u[n-1]*u[n-1];
		ctx->m[sz*(sz-2) + sz-5] = 4*u[n-1]*u[n-1]*u[n-1];
		ctx->m[sz*(sz-2) + sz-4] = 3*u[n-1]*u[n-1];
		ctx->m[sz*(sz-2) + sz-3] = 2*u[n-1];
		ctx->m[sz*(sz-2) + sz-2] = 1;
	}

	ctx->m[sz*(sz-1) + sz-6] = 20*u[n-1]*u[n-1]*u[n-1];
	ctx->m[sz*(sz-1) + sz-5] = 12*u[n-1]*u[n-1];
	ctx->m[sz*(sz-1) + sz-4] = 6*u[n-1];
	ctx->m[sz*(sz-1) + sz-3] = 2;

	if(!mat_lud(ctx)) { return 0; }		

	for(component = 0; component < 3; component++) {
		for(i = 0; i < sz; i++) { ctx->b[i] = 0.0; }
		for(i = 0; i < n-1; i++) {
			ctx->b[2*i + 0] = *VEC3GET(p+i+0, component);
			ctx->b[2*i + 1] = *VEC3GET(p+i+1, component);
		}
		
		if(v0 != NULL) { ctx->b[sz-4] = *VEC3GET(v0, component); }
		if(a0 != NULL) { ctx->b[sz-3] = *VEC3GET(a0, component); }
		if(v1 != NULL) { ctx->b[sz-2] = *VEC3GET(v1, component); }
		if(a1 != NULL) { ctx->b[sz-1] = *VEC3GET(a1, component); }

		mat_lusolve(ctx);

		for(i = 0; i < sz; i++) {
			*VEC3GET(c+i, component) = ctx->b[i];
		}
	}

	return 1;
}

int intrpl_arc(vec3_t v0, vec3_t v1, vec3_t c, double ao, arc_t* arc) {
	vec3_t u, v;

	u = vecsub(v0, c);
	v = vecsub(v1, c);

	if(!APPROACHES_ZERO(vecmag(u) - vecmag(v))) { return 0; }
	if(APPROACHES_ZERO(vecmag(veccross(u, v)))) { return 0; }

	arc->n = vecnorm(veccross(u, v));
	arc->c = c;
	arc->r = u;
	arc->a = ao > 0.0 ? ao : vecangle(u, v);

	return 1;
}
