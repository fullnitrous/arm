#include "approx.h"

chebyshev_t approx_chebyshev_init(double* mem, double a, double b, int n) {
	chebyshev_t ch;
	ch.a    = a;
	ch.b    = b;
	ch.n    = n;
	ch.c    = mem; mem += n;
	ch.cint = mem; mem += n;
	return ch;
}

void approx_chebyshev(APPROX_FUNC1, void* f_, chebyshev_t* apx) {
	
	double vals[APPROX_CHEBYSHEV_MAX_K];
	double bma, bpa, fac, sum, a, b;
	int i, k, n;
	
	n = apx->n;
	a = apx->a;
	b = apx->b;

	bma = 0.5*(b-a);
	bpa = 0.5*(b+a);
	fac = 2.0/n;
	
	for(k = 0; k < n; k++) {
		vals[k] = f(f_, cos(M_PI*(k+0.5)/n)*bma+bpa);
	}

	for(i = 0; i < n; i++) {
		sum = 0.0;
		for(k = 0; k < n; k++) {
			sum += vals[k]*cos(M_PI*i*(k+0.5)/n);
		}
		apx->c[i] = sum * fac;
	}
	return;
}

void approx_chebyshev_int(chebyshev_t* apx) {
	double a, b, con, fac, sum;
	int i, n;

	n = apx->n;
	a = apx->a;
	b = apx->b;

	con = 0.25*(b-a);
	fac = 1.0;
	sum = 0.0;
	for(i = 1; i <= n-2; i++) {
		apx->cint[i] = con * (apx->c[i-1] - apx->c[i+1]) / i;
		sum += fac * apx->cint[i];
		fac = - fac;
	}
	apx->cint[n-1] = con * apx->c[n-2] / (n-1);
	sum += fac * apx->cint[n-1];
	apx->cint[0] = 2.0 * sum;
	return;
}

double approx_eval_chebyshev(chebyshev_t* apx, int der, double x) {
	double y0, y1, d0, d1, tmp, a, b, *c;
	int j, n;

	assert(apx->a <= x && x <= apx->b);

	a = apx->a;
	b = apx->b;
	n = apx->n;

	switch(der) {
		case -1: c = apx->cint; break;
		default: c = apx->c; break;
	}

	y0 = (2.0*x-a-b)*(1.0/(b-a));
	y1 = 2.0*y0;
	d0 = c[n-1];
	d1 = 0.0;

	for(j = n-2; j > 0; j--) {
		tmp = d0;
		d0  = y1*d0 - d1 + c[j];
		d1  = tmp;
	}

	return y0*d0 - d1 + 0.5*c[0];
}

double approx_glq32(APPROX_FUNC1, void* f_, double a, double b) {
	double d, s, q;

	d = b - a;
	s = a + b;
	q = 0.0;

	q += 0.0965400885147278 * f(f_, 0.5*(d*-0.0483076656877383 + s));
	q += 0.0965400885147278 * f(f_, 0.5*(d* 0.0483076656877383 + s));
	q += 0.0956387200792749 * f(f_, 0.5*(d*-0.1444719615827965 + s));
	q += 0.0956387200792749 * f(f_, 0.5*(d* 0.1444719615827965 + s));
	q += 0.0938443990808046 * f(f_, 0.5*(d*-0.2392873622521371 + s));
	q += 0.0938443990808046 * f(f_, 0.5*(d* 0.2392873622521371 + s));
	q += 0.0911738786957639 * f(f_, 0.5*(d*-0.3318686022821277 + s));
	q += 0.0911738786957639 * f(f_, 0.5*(d* 0.3318686022821277 + s));
	q += 0.0876520930044038 * f(f_, 0.5*(d*-0.4213512761306353 + s));
	q += 0.0876520930044038 * f(f_, 0.5*(d* 0.4213512761306353 + s));
	q += 0.0833119242269467 * f(f_, 0.5*(d*-0.5068999089322294 + s));
	q += 0.0833119242269467 * f(f_, 0.5*(d* 0.5068999089322294 + s));
	q += 0.0781938957870703 * f(f_, 0.5*(d*-0.5877157572407623 + s));
	q += 0.0781938957870703 * f(f_, 0.5*(d* 0.5877157572407623 + s));
	q += 0.0723457941088485 * f(f_, 0.5*(d*-0.6630442669302152 + s));
	q += 0.0723457941088485 * f(f_, 0.5*(d* 0.6630442669302152 + s));
	q += 0.0658222227763618 * f(f_, 0.5*(d*-0.7321821187402897 + s));
	q += 0.0658222227763618 * f(f_, 0.5*(d* 0.7321821187402897 + s));
	q += 0.0586840934785355 * f(f_, 0.5*(d*-0.7944837959679424 + s));
	q += 0.0586840934785355 * f(f_, 0.5*(d* 0.7944837959679424 + s));
	q += 0.0509980592623762 * f(f_, 0.5*(d*-0.8493676137325700 + s));
	q += 0.0509980592623762 * f(f_, 0.5*(d* 0.8493676137325700 + s));
	q += 0.0428358980222267 * f(f_, 0.5*(d*-0.8963211557660521 + s));
	q += 0.0428358980222267 * f(f_, 0.5*(d* 0.8963211557660521 + s));
	q += 0.0342738629130214 * f(f_, 0.5*(d*-0.9349060759377397 + s));
	q += 0.0342738629130214 * f(f_, 0.5*(d* 0.9349060759377397 + s));
	q += 0.0253920653092621 * f(f_, 0.5*(d*-0.9647622555875064 + s));
	q += 0.0253920653092621 * f(f_, 0.5*(d* 0.9647622555875064 + s));
	q += 0.0162743947309057 * f(f_, 0.5*(d*-0.9856115115452684 + s));
	q += 0.0162743947309057 * f(f_, 0.5*(d* 0.9856115115452684 + s));
	q += 0.0070186100094701 * f(f_, 0.5*(d*-0.9972638618494816 + s));
	q += 0.0070186100094701 * f(f_, 0.5*(d* 0.9972638618494816 + s));

	return q * 0.5 * d;
}

double approx_solve(APPROX_FUNC2, void* f_, solvecfg_t* cfg, double y) {
	vec2_t o;
	double a, b, x, fx, fa, tol; 
	int i, n, exit;
	
	x   = 0.0;
	a   = cfg->a;
	b   = cfg->b;
	tol = cfg->tol;
	
	if(cfg->only_bisect) {
		n = cfg->bi_nmax;
		exit = 1;
	} else {
		n = cfg->bi_guess_n;
		exit = 0;
	}

BISECT:
	
	for(i = 0; i < n; i++) {
		x  = (a + b) / 2.0;
		fx = f(f_, x).x - y;
		fa = f(f_, a).x - y;

		if(fx == 0.0 || (b-a)/2.0 < tol) {
			return x;
		} else if(fx*fa < 0) {
			b = x;
		} else {
			a = x;
		}
	}
	
	if(exit) { return x; }
		
	n = cfg->nr_n;
	
	for(i = 0; i < n-1; i++) {
		o  = f(f_, x);
		x -= (o.x - y) / o.y;
	}

	/* reuse fx to store next interation */
	/* reuse fa to store difference */
	o  = f(f_, x);
	fx = x - (o.x - y) / o.y;
	fa = fabs(fx-x);
		
	if(fa == 0.0 || fa < tol) {
		return fx;
	} else {
		n = cfg->bi_nmax - cfg->bi_guess_n;
		exit = 1;
		goto BISECT;
	}
}
