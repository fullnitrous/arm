#include "matrix.h"

int mat_lud(luctx_t* ctx) {
	int i, j, k, n, n_minus_1, *pivot;
	double dtmp1, *dptr1, *dptr2, *matrix;

	pivot     = ctx->p;
	matrix    = ctx->m;
	n         = ctx->n;
	n_minus_1 = n-1;	

	for(i = 0; i < n_minus_1; i++) {
		*(pivot + i) = i;
	}
	*(pivot + n_minus_1) = 1;

	for(i = 0; i < n_minus_1; i++) {
		k = i;
		dptr1 = dptr2 = matrix + i*(n+1);
		for(j = i+1; j < n; j++) {
			dptr2 += n;
			if(fabs(*dptr2) > fabs(*dptr1)) {
				dptr1 = dptr2;
				k     = j;
			}
		}
		if(APPROACHES_ZERO(*dptr1)) { return 0; }
		if(k != i) {
			*(pivot + i) = k;
			*(pivot + n_minus_1) = -*(pivot + n_minus_1);
			dptr1 = matrix + i*n;
			dptr2 = matrix + k*n;
			for(j = 0; j < n; j++, dptr1++, dptr2++) {
				dtmp1  = *dptr1;
				*dptr1 = *dptr2;
				*dptr2 = dtmp1;
			}
		}
		for(j = i+1; j < n; j++) {
			dtmp1 = *(matrix + j*n + i) / *(matrix + i*n + i);
			dptr1 = matrix + j*n + i + 1;
			dptr2 = matrix + i*n + i + 1;
			for(k = i+1; k < n; k++, dptr1++, dptr2++) {
				*dptr1 -= dtmp1 * *dptr2;
			}
			*(matrix + j*n + i) = dtmp1;
		}
	}

	return 1;
}

void mat_lusolve(luctx_t* ctx) {
	int i, j, n, n_minus_1, *pivot;
	double dtmp1, *dptr1, *matrix, *b;

	pivot     = ctx->p;
	matrix    = ctx->m;
	b         = ctx->b;
	n         = ctx->n;
	n_minus_1 = n-1;
	
	for(i = 0; i < n_minus_1; i++) {
		j = *(pivot + i);
		if(j != i) {
			dtmp1  = *(b+i);
			*(b+i) = *(b+j);
			*(b+j) = dtmp1;
		}
	}

	for(i = 0; i < n; i++) {
		dtmp1 = *(b+i);
		dptr1 = matrix + i*n;
		for(j = 0; j < i; j++, dptr1++) {
			dtmp1 -= *dptr1 * *(b+j);
		}
		*(b+i) = dtmp1;
	}

	for(i = n_minus_1; i >= 0; i--) {
		dtmp1 = *(b+i);
		dptr1 = matrix + i*n + n_minus_1;
		for(j = n_minus_1; j > i; j--) {
			dtmp1 -= *dptr1 * *(b+j);
			dptr1--;
		}
		*(b+i) = dtmp1 / *dptr1;
	}
	
	return;
}
