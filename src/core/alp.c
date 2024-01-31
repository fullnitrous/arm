#include "alp.h"

double (*alp_integrator)(APPROX_FUNC1, void*, double, double) = &approx_glq32;

double alp_speed(void* f_, double t) {
	return eval_e3d((e3d_t*)f_, t).m;
}

vec2_t ttls(void* data, double t) {
	vec2_t o;
	alpfctx_t* ctx;
	ctx = (alpfctx_t*)data;
	o.x = alp_integrator(&alp_speed, ctx->f, ctx->a, t);
	o.y = alp_speed(ctx->f, t);
	return o;
}

double ltt(void* data, double l) {
	alpfctx_t* ctx;
	
	ctx = (alpfctx_t*)data;

	if(l == 0.0) { return ctx->a; }
	if(l == 1.0) { return ctx->b; }
	
	ctx->cfg->a = ctx->a;
	ctx->cfg->b = ctx->b;
		
	return approx_solve(&ttls, data, ctx->cfg, l);
}

double lttder(void* data, double l) {
	return 1.0/alp_speed(((alpfctx_t*)data)->f, ltt(data, l));
}

alp_t alp_init(double* mem, int n_upper) {
	alp_t alp;

	if(mem == NULL) {
		alp.u = alp.c.a = NULL;
		return alp;
	}
	alp.u   = mem; mem += n_upper + 1;
	alp.c.a = mem;

#ifdef ALP_PERF
	mem += ALP_CHK*n_upper;
	alp.a_der = mem;
#endif

	return alp;
}

#ifdef ALP_PERF
double alp_eval_v(alp_t* alp, double l) {
	double t, tder;
	int k;
	chebyshev_t ch;

	if(alp->u == NULL) {
		t = l / alp->c.s;
		tder = 1.0 / alp->c.s;
		return alp_speed(alp->f, t) * tder;
	}
	
	k = eval_findinterval(alp->u, l, alp->n);
	
	ch.n    = ALP_CHK;
	ch.a    = 0.0;
	ch.b    = alp->u[k+1] - alp->u[k];
	ch.c    = alp->a_der + ALP_CHK * k;
	ch.cint = alp->c.a + ALP_CHK * k;
	
	t    = approx_eval_chebyshev(&ch, -1, l - alp->u[k]);
	tder = approx_eval_chebyshev(&ch,  0, l - alp->u[k]);
	
	return alp_speed(alp->f, t) * tder;
}
#endif

double alp_eval(alp_t* alp, double l) {
	int k;
	chebyshev_t ch;

	if(alp->u == NULL) { return l / alp->c.s; }
	
	k = eval_findinterval(alp->u, l, alp->n);
	
	ch.n = ALP_CHK;
	ch.a = 0.0;
	ch.b = alp->u[k+1] - alp->u[k];
	ch.c = alp->c.a + ALP_CHK * k;

	return approx_eval_chebyshev(&ch, 0, l - alp->u[k]);
}

double alp_fit(e3d_t* f, alp_t* alp, int n_upper, double tolerance) {
	int i, j;
	double l, t, tder, vdev, arc_length, maxdev, guess;
	alpfctx_t fctx;
	solvecfg_t solvecfg;
	chebyshev_t ch;
	double chmem[APPROX_CHEBYSHEV_MEMSIZE_WITH_INT(ALP_CHK)];

#ifdef ALP_PERF
	alp->f = f;
#endif

	ch         = approx_chebyshev_init(chmem, 0.0, 0.0, ALP_CHK);
	arc_length = 0.0;

	if(f->type == E3D_PPOLY1) {
		j = eval_e3dnintervals(f);
		for(i = 0; i < j-1; i++) {
			f->k = i;
			arc_length += (eval_e3dbound(f, i+1)-eval_e3dbound(f, i))*__eval_e3d(f, 0.0).m;
		}
		alp->u = NULL;
		alp->c.s = arc_length;
		return arc_length;
	}

	if(f->type == E3D_ARC) {
		arc_length = eval_e3d(f, 0.5).m;
		alp->u = NULL;
		alp->c.s = arc_length;
		return arc_length;
	}

	if(f->type == E3D_BEZIER) {
		guess = 1.0;
	} else {
		guess = ALP_GUESS_FAC / eval_e3dnintervals(f);
		guess = guess < 1.0 ? guess : 1.0;
	}

	solvecfg.tol         = ALP_SOLVE_TOL;
	solvecfg.bi_guess_n  = ALP_SOLVE_GUESS_N;
	solvecfg.bi_nmax     = ALP_SOLVE_NMAX;
	solvecfg.only_bisect = 0;
	solvecfg.nr_n        = ALP_SOLVE_NR_N;

	fctx.cfg = &solvecfg;
	fctx.f   = f;
	
	*(alp->u) = 0.0;
	j = 0;
	fctx.a = 0.0;
	fctx.b = guess;
	
	while(fctx.a < 1.0) {
		maxdev = 0.0;
		fctx.l = alp_integrator(&alp_speed, f, fctx.a, fctx.b);
		ch.b   = fctx.l;
	
		approx_chebyshev(&lttder, (void*)&fctx, &ch);
		approx_chebyshev_int(&ch);

		for(i = 0; i < ALP_CHECKRES; i++) {
			l    = range(ch.a, ch.b, i, ALP_CHECKRES);
			t    = fctx.a + approx_eval_chebyshev(&ch, -1, l);
			tder = approx_eval_chebyshev(&ch, 0, l);
	
			vdev = fabs(alp_speed(f, t) * tder - 1.0);
			if(vdev > maxdev) { maxdev = vdev; }
		}

		if(fctx.b - fctx.a < ALP_MIN_INTERVAL) {
			return -2.0;
		}
			
		if(maxdev > tolerance) {
			fctx.b = fctx.a + (fctx.b - fctx.a) / 2.0;
		} else if(j < n_upper) {	
			*(alp->u+1+j) = arc_length + fctx.l;
			for(i = 0; i < ALP_CHK; i++) {
				*(alp->c.a+j*ALP_CHK + i) = *(ch.cint + i);
#ifdef ALP_PERF
				*(alp->a_der+j*ALP_CHK + i) = *(ch.c + i);
#endif
			}
			*(alp->c.a+j*ALP_CHK) += 2.0 * fctx.a;

			arc_length += fctx.l;

			fctx.a = fctx.b;
			fctx.b = fctx.b + guess < 1.0 ? fctx.b + guess : 1.0;

			j++;
		} else {
			return -1.0;
		}
	}
	
	alp->n = j + 1;

	return arc_length;
}
