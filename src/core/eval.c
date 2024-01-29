#include "eval.h"

int eval_bico_lut[EVAL_BEZIER_MAX_DEG+1][EVAL_BEZIER_MAX_DEG+1];

int eval_findinterval(double* u, double t, int n) {
	int lo, hi, mi;

	if(t >= u[n-1]) { return n - 2; }
	if(t <= u[0])   { return 0;     }
	
	lo = 0;
	hi = n-1;
	mi = (lo + hi) >> 1;

	while(t < u[mi] || t >= u[mi+1]) {
		if(t < u[mi]) { hi = mi; }
		else          { lo = mi; }
		
		mi = (lo + hi) >> 1;
	}

	return mi;
}

int eval_findspan(double* u, double t, int n, int d) {
	int lo, hi, mi;

    if(t >= u[n]) { return n - 1; }
    if(t <= u[d]) { return d;     }

	lo = d;
	hi = n;
	mi = (lo + hi) >> 1;
    
	while(t < u[mi] || t >= u[mi+1]) {
        if(t < u[mi]) { hi = mi; }
        else          { lo = mi; }
        mi = (lo + hi) >> 1;
    }   
    return mi;
}

ppoly_t eval_ppoly_init(vec3_t* mem0, double* mem1, int n) {
	ppoly_t ppoly;
	ppoly.c = mem0;
	ppoly.u = mem1;
	ppoly.n = n;
	return ppoly;
}

luctx_t eval_ppoly_luctx_init(int* pivotmem, double* matrixmem, double* bmem) {
	luctx_t ctx;
	ctx.n = 0;
	ctx.p = pivotmem;
	ctx.m = matrixmem;
	ctx.b = bmem;
	return ctx;
}

bezier_t eval_bezier_init(vec3_t* mem, int n) {
	bezier_t bezier;
	bezier.p = mem;
	bezier.n = n;
	return bezier;
}

bspline_t eval_bspline_init(vec3_t* mem0, double* mem1, int d, int n) {
	bspline_t bspline;
	bspline.d = d;
	bspline.n = n;
	bspline.p = mem0;
	bspline.u = mem1;
	return bspline;
}

nurbs_t eval_nurbs_init(vec4_t* mem0, double* mem1, int d, int n) {
	nurbs_t nurbs;
	nurbs.d = d;
	nurbs.n = n;
	nurbs.p = mem0;
	nurbs.u = mem1;
	return nurbs;
}

eval_t eval_ppoly1(void* f_, double t) {
	ppoly_t* f = (ppoly_t*)f_;
	int k = eval_findinterval(f->u, t, f->n);
	return __eval_ppoly1(f, k, t);
}

eval_t __eval_ppoly1(void* f_, int k, double t) {
	vec3_t* c = ((ppoly_t*)f_)->c;
	eval_t e;
	
	int aidx = (k << 1);
	int bidx = aidx + 1;

	e.p.x = c[aidx].x*t + c[bidx].x;
	e.p.y = c[aidx].y*t + c[bidx].y;
	e.p.z = c[aidx].z*t + c[bidx].z; 
	
	e.v.x = c[aidx].x;
	e.v.y = c[aidx].y;
	e.v.z = c[aidx].z;

	e.a = veczeros();
	e.m = vecmag(e.v);

	return e;
}

vec3_t __eval_ppoly1_fast(void* f_, int k, double t) {
	vec3_t* c = ((ppoly_t*)f_)->c;
	vec3_t e;
	
	int aidx = (k << 1);
	int bidx = aidx + 1;

	e.x = c[aidx].x*t + c[bidx].x;
	e.y = c[aidx].y*t + c[bidx].y;
	e.z = c[aidx].z*t + c[bidx].z; 
	
	return e;
}

eval_t eval_ppoly3(void* f_, double t) {
	ppoly_t* f = (ppoly_t*)f_;
	int k = eval_findinterval(f->u, t, f->n);
	return __eval_ppoly3(f, k, t);
}

eval_t __eval_ppoly3(void* f_, int k, double t) {
	ppoly_t* f = (ppoly_t*)f_;

	vec3_t* c = f->c;
	
	eval_t e;
	
	int aidx = (k << 2);
	int bidx = aidx + 1;
	int cidx = aidx + 2;
	int didx = aidx + 3;
		
	double t2 = t*t;
	double t3 = t2*t;

	e.p.x = c[aidx].x*t3 + c[bidx].x*t2 + c[cidx].x*t + c[didx].x;
	e.p.y = c[aidx].y*t3 + c[bidx].y*t2 + c[cidx].y*t + c[didx].y;
	e.p.z = c[aidx].z*t3 + c[bidx].z*t2 + c[cidx].z*t + c[didx].z;

	e.v.x = 3.0*c[aidx].x*t2 + 2.0*c[bidx].x*t + c[cidx].x;
	e.v.y = 3.0*c[aidx].y*t2 + 2.0*c[bidx].y*t + c[cidx].y;
	e.v.z = 3.0*c[aidx].z*t2 + 2.0*c[bidx].z*t + c[cidx].z;

	e.a.x = 6.0*c[aidx].x*t + 2.0*c[bidx].x;
	e.a.y = 6.0*c[aidx].y*t + 2.0*c[bidx].y; 
	e.a.z = 6.0*c[aidx].z*t + 2.0*c[bidx].z; 

	e.m = vecmag(e.v);

	return e;
}

vec3_t __eval_ppoly3_fast(void* f_, int k, double t) {
	ppoly_t* f = (ppoly_t*)f_;

	vec3_t* c = f->c;

	vec3_t e;
	
	int aidx = (k << 2);
	int bidx = aidx + 1;
	int cidx = aidx + 2;
	int didx = aidx + 3;
		
	double t2 = t*t;
	double t3 = t2*t;

	e.x = c[aidx].x*t3 + c[bidx].x*t2 + c[cidx].x*t + c[didx].x;
	e.y = c[aidx].y*t3 + c[bidx].y*t2 + c[cidx].y*t + c[didx].y;
	e.z = c[aidx].z*t3 + c[bidx].z*t2 + c[cidx].z*t + c[didx].z;

	return e;
}

eval_t eval_ppoly5(void* f_, double t) {
	ppoly_t* f = (ppoly_t*)f_;
	int k = eval_findinterval(f->u, t, f->n);
	return __eval_ppoly5(f, k, t);
}

eval_t __eval_ppoly5(void* f_, int k, double t) {
	ppoly_t* f = (ppoly_t*)f_;

	vec3_t* c = f->c;

	eval_t e;
	
	int aidx = 6 * k;
	int bidx = aidx + 1;
	int cidx = aidx + 2;
	int didx = aidx + 3;
	int eidx = aidx + 4;
	int fidx = aidx + 5;
	
	double t2 = t*t;
	double t3 = t2*t;
	double t4 = t3*t;
	double t5 = t4*t;

	e.p.x = c[aidx].x*t5 + c[bidx].x*t4 + c[cidx].x*t3 + c[didx].x*t2 + c[eidx].x*t + c[fidx].x;
	e.p.y = c[aidx].y*t5 + c[bidx].y*t4 + c[cidx].y*t3 + c[didx].y*t2 + c[eidx].y*t + c[fidx].y;
	e.p.z = c[aidx].z*t5 + c[bidx].z*t4 + c[cidx].z*t3 + c[didx].z*t2 + c[eidx].z*t + c[fidx].z;

	e.v.x = 5.0*c[aidx].x*t4 + 4.0*c[bidx].x*t3 + 3.0*c[cidx].x*t2 + 2.0*c[didx].x*t + c[eidx].x;
	e.v.y = 5.0*c[aidx].y*t4 + 4.0*c[bidx].y*t3 + 3.0*c[cidx].y*t2 + 2.0*c[didx].y*t + c[eidx].y;
	e.v.z = 5.0*c[aidx].z*t4 + 4.0*c[bidx].z*t3 + 3.0*c[cidx].z*t2 + 2.0*c[didx].z*t + c[eidx].z;

	e.a.x = 20.0*c[aidx].x*t3 + 12.0*c[bidx].x*t2 + 6.0*c[cidx].x*t + 2.0*c[didx].x;
	e.a.y = 20.0*c[aidx].y*t3 + 12.0*c[bidx].y*t2 + 6.0*c[cidx].y*t + 2.0*c[didx].y;
	e.a.z = 20.0*c[aidx].z*t3 + 12.0*c[bidx].z*t2 + 6.0*c[cidx].z*t + 2.0*c[didx].z;

	e.m = vecmag(e.v);

	return e;
}

vec3_t __eval_ppoly5_fast(void* f_, int k, double t) {
	ppoly_t* f = (ppoly_t*)f_;

	vec3_t* c = f->c;

	vec3_t e;
	
	int aidx = 6 * k;
	int bidx = aidx + 1;
	int cidx = aidx + 2;
	int didx = aidx + 3;
	int eidx = aidx + 4;
	int fidx = aidx + 5;
	
	double t2 = t*t;
	double t3 = t2*t;
	double t4 = t3*t;
	double t5 = t4*t;

	e.x = c[aidx].x*t5 + c[bidx].x*t4 + c[cidx].x*t3 + c[didx].x*t2 + c[eidx].x*t + c[fidx].x;
	e.y = c[aidx].y*t5 + c[bidx].y*t4 + c[cidx].y*t3 + c[didx].y*t2 + c[eidx].y*t + c[fidx].y;
	e.z = c[aidx].z*t5 + c[bidx].z*t4 + c[cidx].z*t3 + c[didx].z*t2 + c[eidx].z*t + c[fidx].z;

	return e;
}


eval_t eval_arc(void* f_, double t) {
	arc_t* f = (arc_t*)f_;
	
	vec3_t r = f->r;
	vec3_t n = f->n;
	vec3_t c = f->c;
	double a = t * f->a;

	eval_t e;

	double cosa           = cos(a);
	double one_minus_cosa = 1.0 - cosa;
	double sina           = sin(a);
	double minus_sina     = -sina;
	double minus_cosa     = -cosa;
	double rdotn          = r.x*n.x + r.y*n.y + r.z*n.z;

	double ncrossr_x = n.y*r.z - n.z*r.y;
	double ncrossr_y = n.z*r.x - n.x*r.z;
	double ncrossr_z = n.x*r.y - n.y*r.x;
	
	e.p.x = r.x*cosa + one_minus_cosa*rdotn*n.x + sina*ncrossr_x + c.x;
    e.p.y = r.y*cosa + one_minus_cosa*rdotn*n.y + sina*ncrossr_y + c.y;
    e.p.z = r.z*cosa + one_minus_cosa*rdotn*n.z + sina*ncrossr_z + c.z;

    e.v.x = f->a*(r.x*minus_sina + sina*rdotn*n.x + cosa*ncrossr_x);
    e.v.y = f->a*(r.y*minus_sina + sina*rdotn*n.y + cosa*ncrossr_y);
    e.v.z = f->a*(r.z*minus_sina + sina*rdotn*n.z + cosa*ncrossr_z);
	
	e.a.x = f->a*f->a*(r.x*minus_cosa + cosa*rdotn*n.x + minus_sina*ncrossr_x);
	e.a.y = f->a*f->a*(r.y*minus_cosa + cosa*rdotn*n.y + minus_sina*ncrossr_y);
	e.a.z = f->a*f->a*(r.z*minus_cosa + cosa*rdotn*n.z + minus_sina*ncrossr_z);

	e.m = vecmag(e.v);

	return e;
}

vec3_t __eval_arc_fast(void* f_, int k, double t) {
	arc_t* f = (arc_t*)f_;
	
	vec3_t r = f->r;
	vec3_t n = f->n;
	vec3_t c = f->c;
	double a = t * f->a;

	vec3_t e;

	double cosa           = cos(a);
	double one_minus_cosa = 1.0 - cosa;
	double sina           = sin(a);
	double rdotn          = r.x*n.x + r.y*n.y + r.z*n.z;

	double ncrossr_x = n.y*r.z - n.z*r.y;
	double ncrossr_y = n.z*r.x - n.x*r.z;
	double ncrossr_z = n.x*r.y - n.y*r.x;
	
	e.x = r.x*cosa + one_minus_cosa*rdotn*n.x + sina*ncrossr_x + c.x;
    e.y = r.y*cosa + one_minus_cosa*rdotn*n.y + sina*ncrossr_y + c.y;
    e.z = r.z*cosa + one_minus_cosa*rdotn*n.z + sina*ncrossr_z + c.z;

	return e;
}

void eval_bezier_lut_init(void) {
	int n, i;

	for(n = 0; n <= EVAL_BEZIER_MAX_DEG; n++) {
		for(i = 0; i <= n; i++) {
			eval_bico_lut[n][i] = __eval_bico(n, i);
		}
	}

	return;
}

int __eval_fac(int n) {
	int y, i;
	for(i = y = 1; i <= n; y *= (i++));
	return y;
}

int __eval_bico(int n, int i) {
	return __eval_fac(n) / (__eval_fac(i)*__eval_fac(n-i));
}

int eval_bico(int n, int i) {
	if(n <= EVAL_BEZIER_MAX_DEG) { return eval_bico_lut[n][i]; }
	return __eval_bico(n, i);
}

double eval_bernstein_basis(int i, int n, double t) {
	return eval_bico(n, i) * pow(t, i) * pow(1.0 - t, n - i);
}

eval_t eval_bezier(void* f_, double t) {
	bezier_t* f = (bezier_t*)f_;
	
	vec3_t* p = f->p;
	int n = f->n;

	eval_t e;
	int i;
	double x, y, z, b;
	
	for(i = 0, x = y = z = 0.0; i <= n; i++) {
		b = eval_bernstein_basis(i, n, t);
		x += b*p[i].x;
		y += b*p[i].y;
		z += b*p[i].z;
	}
	
	e.p.x = x;
	e.p.y = y;
	e.p.z = z;

	for(i = 0, x = y = z = 0.0; i <= n - 1; i++) {
		b = eval_bernstein_basis(i, n - 1, t);
		x += b*(p[i+1].x - p[i].x);
		y += b*(p[i+1].y - p[i].y);
		z += b*(p[i+1].z - p[i].z);
	}

	e.v.x = n*x;
	e.v.y = n*y;
	e.v.z = n*z;

	e.m = vecmag(e.v);

	for(i = 0, x = y = z = 0.0; i <= n - 2; i++) {
		b = eval_bernstein_basis(i, n - 2, t);
		x += b*(p[i+2].x - 2*p[i+1].x + p[i].x);
		y += b*(p[i+2].y - 2*p[i+1].y + p[i].y);
		z += b*(p[i+2].z - 2*p[i+1].z + p[i].z);
	}

	b = (n-1)*n;
	
	e.a.x = b*x;
	e.a.y = b*y;
	e.a.z = b*z;

	return e;
}

vec3_t __eval_bezier_fast(void* f_, int k, double t) {
	bezier_t* f = (bezier_t*)f_;
	
	vec3_t* p = f->p;
	int n = f->n;

	vec3_t e;
	int i;
	double b;
	
	e.x = e.y = e.z = 0.0;

	for(i = 0; i <= n; i++) {
		b = eval_bernstein_basis(i, n, t);
		e.x += b*p[i].x;
		e.y += b*p[i].y;
		e.z += b*p[i].z;
	}
	
	return e;
}

eval_t eval_bspline(void* f_, double t) {
	bspline_t* f;
	double t_, *u;
	int d, n, k;
		
	f  = (bspline_t*)f_;
    u  = f->u;
    d  = f->d;
    n  = f->n;
	t_ = u[d] + t * (u[n] - u[d]);
    k  = eval_findspan(u, t_, n, d);

	return __eval_bspline(f_, k, t);
}

eval_t __eval_bspline(void* f_, int k, double t) {
	bspline_t* f;
	eval_t e;
	double alpha, tmp, *u;
	int i, j, r, d, n, left, right;

	vec3_t p[EVAL_BSPLINE_MAX_DEG+1];
	vec3_t v[EVAL_BSPLINE_MAX_DEG+0];
	vec3_t a[EVAL_BSPLINE_MAX_DEG-1];

	f = (bspline_t*)f_;
	u = f->u;
	d = f->d;
	n = f->n;
	t = u[d] + t * (u[n] - u[d]);

	
    for(i = 0; i <= d; i++) {
        p[i] = f->p[i + k - d];
    }
	
	for(i = 0; i <= d-1; i++) {
		tmp    = d / (u[i+k+1] - u[i+k-d+1]);
		v[i].x = tmp*(p[i+1].x - p[i].x);
		v[i].y = tmp*(p[i+1].y - p[i].y);
		v[i].z = tmp*(p[i+1].z - p[i].z);
	}
	
	for(i = 0; i <= d - 2; i++) {
		tmp    = (d-1) / (u[i+k+1] - u[i+k-d+2]);
		a[i].x = tmp*(v[i+1].x - v[i].x);
		a[i].y = tmp*(v[i+1].y - v[i].y);
		a[i].z = tmp*(v[i+1].z - v[i].z);
	}
		
	for(r = 1; r <= d; r++) {
		for(j = d; j >= r; j--) {
			right = j+1+k-r;
			left  = j+k-d;
			alpha = (t-u[left])/(u[right]-u[left]);
			left++;
			
			p[j].x = (1.0-alpha)*p[j-1].x + alpha*p[j].x;
			p[j].y = (1.0-alpha)*p[j-1].y + alpha*p[j].y;
			p[j].z = (1.0-alpha)*p[j-1].z + alpha*p[j].z;

			if(r <= d-1 && j <= d-1) {
				alpha  = (t-u[left])/(u[right]-u[left]);
				left++;
				v[j].x = (1.0-alpha)*v[j-1].x + alpha*v[j].x;
				v[j].y = (1.0-alpha)*v[j-1].y + alpha*v[j].y;
				v[j].z = (1.0-alpha)*v[j-1].z + alpha*v[j].z;
			}
			
			if(r <= d-2 && j <= d-2) {
				alpha  = (t-u[left])/(u[right]-u[left]);
				a[j].x = (1.0-alpha)*a[j-1].x + alpha*a[j].x;
				a[j].y = (1.0-alpha)*a[j-1].y + alpha*a[j].y;
				a[j].z = (1.0-alpha)*a[j-1].z + alpha*a[j].z;
			}
		}
    }

    e.p = p[d];
    
	if(d >= 1) { e.v = v[d-1];     }
	else       { e.v = veczeros(); }
	
	if(d >= 2) { e.a = a[d-2];     }
	else       { e.a = veczeros(); }

    e.m = vecmag(e.v);
    
	return e;
}

vec3_t __eval_bspline_fast(void* f_, int k, double t) {
	bspline_t* f;
	double alpha, *u;
	int i, j, r, d, n, left, right;

	vec3_t p[EVAL_BSPLINE_MAX_DEG+1];

	f = (bspline_t*)f_;
	u = f->u;
	d = f->d;
	n = f->n;
	t = u[d] + t * (u[n] - u[d]);
	
    for(i = 0; i <= d; i++) {
        p[i] = f->p[i + k - d];
    }
		
	for(r = 1; r <= d; r++) {
		for(j = d; j >= r; j--) {
			right = j+1+k-r;
			left  = j+k-d;
			alpha = (t-u[left])/(u[right]-u[left]);
			left++;
			
			p[j].x = (1.0-alpha)*p[j-1].x + alpha*p[j].x;
			p[j].y = (1.0-alpha)*p[j-1].y + alpha*p[j].y;
			p[j].z = (1.0-alpha)*p[j-1].z + alpha*p[j].z;
		}
    }

	return p[d];
}

eval_t eval_nurbs(void* f_, double t) {
	nurbs_t* f;
	double t_, *u;
	int d, n, k;
		
	f  = (nurbs_t*)f_;
    u  = f->u;
    d  = f->d;
    n  = f->n;
	t_ = u[d] + t * (u[n] - u[d]);
    k  = eval_findspan(u, t_, n, d);

	return __eval_nurbs(f_, k, t);
}

eval_t __eval_nurbs(void* f_, int k, double t) {
	nurbs_t* f;
	eval_t e;
	double alpha, tmp, *u;
	int i, j, r, d, n, left, right;

	vec4_t p[EVAL_NURBS_MAX_DEG+1];
    vec4_t v[EVAL_NURBS_MAX_DEG+0];
	vec4_t a[EVAL_NURBS_MAX_DEG-1];

	f = (nurbs_t*)f_;
	u = f->u;
	d = f->d;
	n = f->n;
	t = u[d] + t * (u[n] - u[d]);

    for(i = 0; i <= d; i++) {
		vec4_t v = f->p[i+k-d];
        p[i].x   = v.x * v.w;
		p[i].y   = v.y * v.w;
		p[i].z   = v.z * v.w;
		p[i].w   = v.w;
    }

	for(i = 0; i <= d-1; i++) {
		tmp    = d / (u[i+k+1] - u[i+k-d+1]);
		v[i].x = tmp*(p[i+1].x - p[i].x);
		v[i].y = tmp*(p[i+1].y - p[i].y);
		v[i].z = tmp*(p[i+1].z - p[i].z);
		v[i].w = tmp*(p[i+1].w - p[i].w);
	}
	
	for(i = 0; i <= d - 2; i++) {
		tmp    = (d-1) / (u[i+k+1] - u[i+k-d+2]);
		a[i].x = tmp*(v[i+1].x - v[i].x);
		a[i].y = tmp*(v[i+1].y - v[i].y);
		a[i].z = tmp*(v[i+1].z - v[i].z);
		a[i].w = tmp*(v[i+1].w - v[i].w);
	}
	
	for(r = 1; r <= d; r++) {
		for(j = d; j >= r; j--) {
			right = j+1+k-r;
			left  = j+k-d;
			alpha = (t-u[left])/(u[right]-u[left]);
			left++;
			
			p[j].x = (1.0-alpha)*p[j-1].x + alpha*p[j].x;
			p[j].y = (1.0-alpha)*p[j-1].y + alpha*p[j].y;
			p[j].z = (1.0-alpha)*p[j-1].z + alpha*p[j].z;
			p[j].w = (1.0-alpha)*p[j-1].w + alpha*p[j].w;

			if(r <= d-1 && j <= d-1) {
				alpha  = (t-u[left])/(u[right]-u[left]);
				left++;
				v[j].x = (1.0-alpha)*v[j-1].x + alpha*v[j].x;
				v[j].y = (1.0-alpha)*v[j-1].y + alpha*v[j].y;
				v[j].z = (1.0-alpha)*v[j-1].z + alpha*v[j].z;
				v[j].w = (1.0-alpha)*v[j-1].w + alpha*v[j].w;
			}
			
			if(r <= d-2 && j <= d-2) {
				alpha  = (t-u[left])/(u[right]-u[left]);
				a[j].x = (1.0-alpha)*a[j-1].x + alpha*a[j].x;
				a[j].y = (1.0-alpha)*a[j-1].y + alpha*a[j].y;
				a[j].z = (1.0-alpha)*a[j-1].z + alpha*a[j].z;
				a[j].w = (1.0-alpha)*a[j-1].w + alpha*a[j].w;
			}
		}
    }

	p[d].x /= p[d].w;
	p[d].y /= p[d].w;
	p[d].z /= p[d].w;

	e.p.x = p[d].x;
	e.p.y = p[d].y;
	e.p.z = p[d].z;		

	if(d >= 1) {
		v[d-1].x = (v[d-1].x - v[d-1].w*p[d].x) / p[d].w; 
		v[d-1].y = (v[d-1].y - v[d-1].w*p[d].y) / p[d].w; 
		v[d-1].z = (v[d-1].z - v[d-1].w*p[d].z) / p[d].w;
		e.v.x    = v[d-1].x; 
		e.v.y    = v[d-1].y;
		e.v.z    = v[d-1].z;
	} else {
		e.v = veczeros();
	}
	
	if(d >= 2) {
		e.a.x = (a[d-2].x - 2.0*v[d-1].w*v[d-1].x - a[d-2].w*p[d].x) / p[d].w; 
		e.a.y = (a[d-2].y - 2.0*v[d-1].w*v[d-1].y - a[d-2].w*p[d].y) / p[d].w; 
		e.a.z = (a[d-2].z - 2.0*v[d-1].w*v[d-1].z - a[d-2].w*p[d].z) / p[d].w; 
	} else {
		e.v = veczeros();
	}

	e.m = vecmag(e.v);

	return e;
}

vec3_t __eval_nurbs_fast(void* f_, int k, double t) {
	nurbs_t* f;
	vec3_t e;
	double alpha, *u;
	int i, j, r, d, n, left, right;

	vec4_t p[EVAL_NURBS_MAX_DEG+1];

	f = (nurbs_t*)f_;
	u = f->u;
	d = f->d;
	n = f->n;
	t = u[d] + t * (u[n] - u[d]);

    for(i = 0; i <= d; i++) {
		vec4_t v = f->p[i+k-d];
        p[i].x   = v.x * v.w;
		p[i].y   = v.y * v.w;
		p[i].z   = v.z * v.w;
		p[i].w   = v.w;
    }
	
	for(r = 1; r <= d; r++) {
		for(j = d; j >= r; j--) {
			right = j+1+k-r;
			left  = j+k-d;
			alpha = (t-u[left])/(u[right]-u[left]);
			left++;
			
			p[j].x = (1.0-alpha)*p[j-1].x + alpha*p[j].x;
			p[j].y = (1.0-alpha)*p[j-1].y + alpha*p[j].y;
			p[j].z = (1.0-alpha)*p[j-1].z + alpha*p[j].z;
			p[j].w = (1.0-alpha)*p[j-1].w + alpha*p[j].w;

		}
    }

	e.x = p[d].x / p[d].w;
	e.y = p[d].y / p[d].w;
	e.z = p[d].z / p[d].w;

	return e;
}

int eval_cmp(eval_t tst, eval_t ref, double tol) {
	if(veccmp(tst.p, ref.p, tol)) { return 1; }
	if(veccmp(tst.v, ref.v, tol)) { return 2; }
	if(veccmp(tst.a, ref.a, tol)) { return 3; }
	return 0;
}

eval_t eval_cfd(e3d_t* e3d, double h, double t) {
	eval_t e;
	eval_t tmp;
	double dt;
	
	const static double coeffs_1[9] = {
		-1.0/280.0,
		4.0/105.0,
		-1.0/5.0,
		4.0/5.0,
		0.0,
		-4.0/5.0,
		1.0/5.0,
		-4.0/105.0,
		1.0/280.0
	};
	
	const static double coeffs_2[9] = {
		-1.0/560.0,
		8.0/315.0,
		-1.0/5.0,
		8.0/5.0,
		-205.0/72.0,
		8.0/5.0,
		-1.0/5.0,
		8.0/315.0,
		-1.0/560.0
	};
	
	int i;

	double h_fac = 1.0;

	if(e3d->type == E3D_BSPLINE) {
		bspline_t* bspline = (bspline_t*)e3d->function;
		double* u = bspline->u;
		h_fac = u[bspline->n] - u[bspline->d];
	} else if(e3d->type == E3D_NURBS) {
		nurbs_t* nurbs = (nurbs_t*)e3d->function;
		double* u = nurbs->u;
		h_fac = u[nurbs->n] - u[nurbs->d];
	}
	
	dt  = h / 9.0;
	e.p = eval_e3d(e3d, t).p;
	e.v = veczeros();
	e.a = veczeros();
	
	for(i = 0; i < 9; i++) {
		tmp = eval_e3d(e3d, t + (4-i)*dt);
		e.v = vecadd(e.v, vecmul(tmp.p, coeffs_1[i]));
		e.a = vecadd(e.a, vecmul(tmp.p, coeffs_2[i]));
	}

	e.v = vecdiv(e.v, h_fac*dt);
	e.a = vecdiv(e.a, h_fac*h_fac*dt*dt);
	e.m = vecmag(e.v);

	return e;
}

e3d_t eval_e3d_make(void* function, e3dt_t type) {
	e3d_t e;
	e.type     = type;
	e.function = function;
	switch(type) {
		case E3D_PPOLY1:
			e.evaluator = &eval_ppoly1;
			e.__evaluator = &__eval_ppoly1;
			e.__evaluator_fast = &__eval_ppoly1_fast;
			break;
		case E3D_PPOLY3:
			e.evaluator = &eval_ppoly3;
			e.__evaluator = &__eval_ppoly3;
			e.__evaluator_fast = &__eval_ppoly3_fast;
			break;
		case E3D_PPOLY5:
			e.evaluator = &eval_ppoly5;
			e.__evaluator = &__eval_ppoly5;
			e.__evaluator_fast = &__eval_ppoly5_fast;
			break;
		case E3D_ARC:
			e.evaluator = &eval_arc;
			e.__evaluator = NULL;
			e.__evaluator_fast = &__eval_arc_fast;
			break;
		case E3D_BEZIER:
			e.evaluator = &eval_bezier;
			e.__evaluator = NULL;
			e.__evaluator_fast = &__eval_bezier_fast;
			break;
		case E3D_BSPLINE:
			e.evaluator = &eval_bspline;
			e.__evaluator = &__eval_bspline;
			e.__evaluator_fast = &__eval_bspline_fast;
			break;
		case E3D_NURBS:
			e.evaluator = &eval_nurbs;
			e.__evaluator = &__eval_nurbs;
			e.__evaluator_fast = &__eval_nurbs_fast;
			break;
		default:
			break;
	}
	return e;
}

eval_t eval_e3d(e3d_t* e3d, double t) {
	return e3d->evaluator(e3d->function, t);
}

eval_t __eval_e3d(e3d_t* e3d, double t) {
	return e3d->__evaluator(e3d->function, e3d->k, t);
}

vec3_t __eval_e3d_fast(e3d_t* e3d, double t) {
	return e3d->__evaluator_fast(e3d->function, e3d->k, t);
}

int eval_e3dnintervals(e3d_t* e3d) {
	switch(e3d->type) {
		case E3D_PPOLY1:
		case E3D_PPOLY3:
		case E3D_PPOLY5:
			return ((ppoly_t*)e3d->function)->n;
		case E3D_BSPLINE:
			return ((bspline_t*)e3d->function)->n - ((bspline_t*)e3d->function)->d;
		case E3D_NURBS:
			return ((nurbs_t*)e3d->function)->n - ((nurbs_t*)e3d->function)->d;
		default:
			return -1;
	}
}

int eval_e3dinterval(e3d_t* e3d, double t) {
	switch(e3d->type) {
		case E3D_PPOLY1:
		case E3D_PPOLY3:
		case E3D_PPOLY5:
			return eval_findinterval(
				((ppoly_t*)(e3d->function))->u,
				t,
				((ppoly_t*)(e3d->function))->n
			);
		case E3D_BSPLINE:
			return eval_findspan(
				((bspline_t*)(e3d->function))->u,
				t,
				((bspline_t*)(e3d->function))->n,
				((bspline_t*)(e3d->function))->d
			);
		case E3D_NURBS:
			return eval_findspan(
				((nurbs_t*)(e3d->function))->u,
				t,
				((nurbs_t*)(e3d->function))->n,
				((nurbs_t*)(e3d->function))->d
			);
		default:
			return -1;
	}
}

double eval_e3dbound(e3d_t* e3d, int k) {
	switch(e3d->type) {
		case E3D_PPOLY1:
		case E3D_PPOLY3:
		case E3D_PPOLY5:
			return ((ppoly_t*)e3d->function)->u[k];
		case E3D_BSPLINE:
			return ((bspline_t*)e3d->function)->u[((bspline_t*)e3d->function)->d + k];
		case E3D_NURBS:
			return ((nurbs_t*)e3d->function)->u[((nurbs_t*)e3d->function)->d + k];
		default:
			return -1;
	}
}

double eval_maxv(eval_t* e, double maxat) {
	double k;
	k = vecmag(veccross(e->v, e->a)) / pow(vecmag(e->v), 3.0);	
	return sqrt(maxat / k);
}
