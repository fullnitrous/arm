#include <eval.h>
#include <intrpl.h>
#include <vector.h>
#include <util.h>

#include <stdlib.h>
#include <stdio.h>

int main(void) {
	eval_bezier_lut_init();

	FILE* pointsf      = fopen("csv/points.csv", "w");
	FILE* arcpf        = fopen("csv/arcp.csv", "w");
	FILE* linearf      = fopen("csv/linear.csv", "w");
	FILE* cubic0f      = fopen("csv/cubic0.csv", "w");
	FILE* cubic1f      = fopen("csv/cubic1.csv", "w");
	FILE* quintic0f    = fopen("csv/quintic0.csv", "w");
	FILE* quintic1f    = fopen("csv/quintic1.csv", "w");
	FILE* arcf         = fopen("csv/arc.csv", "w");
	FILE* bezierf      = fopen("csv/bezier.csv", "w");
	FILE* bezierpolyf  = fopen("csv/bezierpoly.csv", "w");
	FILE* bsplinef     = fopen("csv/bspline.csv", "w");
	FILE* bsplinepolyf = fopen("csv/bsplinepoly.csv", "w");
	FILE* nurbsf       = fopen("csv/nurbs.csv", "w");
	FILE* nurbspolyf   = fopen("csv/nurbspoly.csv", "w");
	
	int n = 12;
	vec3_t p[n];

	int i;

	for(i = 0; i < n; i++) {
		double t = range(0, M_PI, i, n);
		p[i].x = t;
		p[i].y = sin(3*t);
		p[i].z = cos(3*t);

		fprintf(pointsf, "%f, %f, %f\n", p[i].x, p[i].y, p[i].z);
	}
	
	vec3_t linear_mem0[EVAL_PPOLY1_MEMSIZE(n)];
	double linear_mem1[n];
	ppoly_t linear_ = eval_ppoly_init(linear_mem0, linear_mem1, n);
	
	vec3_t cubic0_mem0[EVAL_PPOLY3_MEMSIZE(n)];
	double cubic0_mem1[n];
	ppoly_t cubic0_ = eval_ppoly_init(cubic0_mem0, cubic0_mem1, n);
	
	vec3_t cubic1_mem0[EVAL_PPOLY3_MEMSIZE(n)];
	double cubic1_mem1[n];
	ppoly_t cubic1_ = eval_ppoly_init(cubic1_mem0, cubic1_mem1, n);


	int pivot3_mem[EVAL_PPOLY3_PIVOT_MEMSIZE(n)];
	double matrix3_mem[EVAL_PPOLY3_EQSYS_MEMSIZE(n)];
	double b3_mem[EVAL_PPOLY3_B_MEMSIZE(n)];
	
	vec3_t quintic0_mem0[EVAL_PPOLY5_MEMSIZE(n)];
	double quintic0_mem1[n];
	ppoly_t quintic0_ = eval_ppoly_init(quintic0_mem0, quintic0_mem1, n);
	
	vec3_t quintic1_mem0[EVAL_PPOLY5_MEMSIZE(n)];
	double quintic1_mem1[n];
	ppoly_t quintic1_ = eval_ppoly_init(quintic1_mem0, quintic1_mem1, n);

	int pivot5_mem[EVAL_PPOLY5_PIVOT_MEMSIZE(n)];
	double matrix5_mem[EVAL_PPOLY5_EQSYS_MEMSIZE(n)];
	double b5_mem[EVAL_PPOLY5_B_MEMSIZE(n)];

	luctx_t ctx3 = eval_ppoly_luctx_init(pivot3_mem, matrix3_mem, b3_mem);
	luctx_t ctx5 = eval_ppoly_luctx_init(pivot5_mem, matrix5_mem, b5_mem);

	vec3_t v0 = {0, 0, -30};
	vec3_t v1 = {0, 0, -30};
	
	printf("=== INTERPOLATION ===\n");

	int status1 = intrpl_ppoly1(p, n, &linear_);
	int status3 = intrpl_ppoly3(p, n, NULL, NULL, &ctx3, &cubic0_);
	status3 *= intrpl_ppoly3(p, n, &v0, &v1, &ctx3, &cubic1_);
	
	int status5 = intrpl_ppoly5(p, n, NULL, NULL, NULL, NULL, &ctx5, &quintic0_);
	status5 *= intrpl_ppoly5(p, n, &v0, &v1, NULL, NULL, &ctx5, &quintic1_);


	arc_t arc_;

	v0 = (vec3_t){1, 0, 0};
	v1 = (vec3_t){0, 0, 1};
	vec3_t c  = {0, 0, 0};

	int statusa = intrpl_arc(v0, v1, c, (4.0/3.0)*M_PI, &arc_);

	fprintf(arcpf, "%f, %f, %f\n", v0.x, v0.y, v0.z);
	fprintf(arcpf, "%f, %f, %f\n", v1.x, v1.y, v1.z);
	fprintf(arcpf, "%f, %f, %f\n", c.x, c.y, c.z);

	printf("linear spline:  %s\n", status1 ? "\e[1;32msuccess\e[m" : "\e[0;31mfailure\e[m");
	printf("cubic spline:   %s\n", status3 ? "\e[1;32msuccess\e[m" : "\e[0;31mfailure\e[m");
	printf("quintic spline: %s\n", status5 ? "\e[1;32msuccess\e[m" : "\e[0;31mfailure\e[m");
	printf("arc:            %s\n", statusa ? "\e[1;32msuccess\e[m" : "\e[0;31mfailure\e[m");
	
	int bez_n = 3;
	vec3_t bezier_mem[EVAL_BEZIER_MEMSIZE(bez_n)];
	bezier_t bezier_ = eval_bezier_init(bezier_mem, bez_n);
	
	bezier_.p[0] = (vec3_t){0.0, 0.0, 0.0};
	bezier_.p[1] = (vec3_t){0.0, 0.5, 1.0};
	bezier_.p[2] = (vec3_t){1.0, 0.0, 1.0};
	bezier_.p[3] = (vec3_t){1.0, 0.5, 2.0};

	for(i = 0; i <= bez_n; i++) {
		vec3_t v = bezier_.p[i];
		fprintf(bezierpolyf, "%f, %f, %f\n", v.x, v.y, v.z); 
	}

	vec3_t bspline_mem0[6];
	double bspline_mem1[10];
	bspline_t bspline_ = eval_bspline_init(bspline_mem0, bspline_mem1, 3, 6);

	bspline_.p[0] = (vec3_t){-4.0, -4.0,  0.0};
	bspline_.p[1] = (vec3_t){ 1.0, -1.7,  0.1};
	bspline_.p[2] = (vec3_t){-4.2,  1.3, -2.0};
	bspline_.p[3] = (vec3_t){ 4.6, -1.0, -3.0};
	bspline_.p[4] = (vec3_t){-2.5,  3.5, -4.7};
	bspline_.p[5] = (vec3_t){ 5.4,  5.5, -3.2};

	bspline_.u[0] = 0.0;
	bspline_.u[1] = 0.0;
	bspline_.u[2] = 0.0;
	bspline_.u[3] = 0.0;
	bspline_.u[4] = 0.444;
	bspline_.u[5] = 0.556;
	bspline_.u[6] = 1.0;
	bspline_.u[7] = 1.0;
	bspline_.u[8] = 1.0;
	bspline_.u[9] = 1.0;
	
	for(i = 0; i < bspline_.n; i++) {
		vec3_t v = bspline_.p[i];
		fprintf(bsplinepolyf, "%f, %f, %f\n", v.x, v.y, v.z); 
	}

	vec4_t nurbs_mem0[51];
	double nurbs_mem1[55];
	nurbs_t nurbs_ = eval_nurbs_init(nurbs_mem0, nurbs_mem1, 3, 51);

	nurbs_.p[0]  = (vec4_t){54.4930, 52.1390, 0.0000, 1.0000};
	nurbs_.p[1]  = (vec4_t){55.5070, 52.1390, 0.0000, 1.0000};
	nurbs_.p[2]  = (vec4_t){56.0820, 49.6150, 0.0000, 1.0000};
	nurbs_.p[3]  = (vec4_t){56.7800, 44.9710, 0.0000, 1.2000};
	nurbs_.p[4]  = (vec4_t){69.5750, 51.3580, 0.0000, 1.0000};
	nurbs_.p[5]  = (vec4_t){77.7860, 58.5730, 0.0000, 1.0000};
	nurbs_.p[6]  = (vec4_t){90.5260, 67.0810, 0.0000, 1.0000};
	nurbs_.p[7]  = (vec4_t){105.9730, 63.8010, 0.0000, 1.0000};
	nurbs_.p[8]  = (vec4_t){100.4000, 47.3260, 0.0000, 1.0000};
	nurbs_.p[9]  = (vec4_t){94.5670, 39.9130, 0.0000, 1.0000};
	nurbs_.p[10] = (vec4_t){92.3690, 30.4850, 0.0000, 1.0000};
	nurbs_.p[11] = (vec4_t){83.4400, 33.7570, 0.0000, 2.0000};
	nurbs_.p[12] = (vec4_t){91.8920, 28.5090, 0.0000, 1.0000};
	nurbs_.p[13] = (vec4_t){89.4440, 20.3930, 0.0000, 1.0000};
	nurbs_.p[14] = (vec4_t){83.2180, 15.4460, 0.0000, 5.0000};
	nurbs_.p[15] = (vec4_t){87.6210, 4.8300, 0.0000, 3.0000};
	nurbs_.p[16] = (vec4_t){80.9450, 9.2670, 0.0000, 1.0000};
	nurbs_.p[17] = (vec4_t){79.8340, 14.5350, 0.0000, 1.1000};
	nurbs_.p[18] = (vec4_t){76.0740, 8.5220, 0.0000, 1.0000};
	nurbs_.p[19] = (vec4_t){70.1830, 12.5500, 0.0000, 1.0000};
	nurbs_.p[20] = (vec4_t){64.1710, 16.8650, 0.0000, 1.0000};
	nurbs_.p[21] = (vec4_t){59.9930, 22.1220, 0.0000, 1.0000};
	nurbs_.p[22] = (vec4_t){55.6800, 36.3590, 0.0000, 1.0000};
	nurbs_.p[23] = (vec4_t){56.9250, 24.9950, 0.0000, 1.0000};
	nurbs_.p[24] = (vec4_t){59.7650, 19.8280, 0.0000, 1.0000};
	nurbs_.p[25] = (vec4_t){54.4930, 14.9400, 0.0000, 1.0000};
	nurbs_.p[26] = (vec4_t){49.2200, 19.8280, 0.0000, 1.0000};
	nurbs_.p[27] = (vec4_t){52.0600, 24.9940, 0.0000, 1.0000};
	nurbs_.p[28] = (vec4_t){53.3050, 36.3590, 0.0000, 1.0000};
	nurbs_.p[29] = (vec4_t){48.9920, 22.1220, 0.0000, 1.0000};
	nurbs_.p[30] = (vec4_t){44.8140, 16.8650, 0.0000, 1.0000};
	nurbs_.p[31] = (vec4_t){38.8020, 12.5510, 0.0000, 1.0000};
	nurbs_.p[32] = (vec4_t){32.9110, 8.5210, 0.0000, 1.0000};
	nurbs_.p[33] = (vec4_t){29.1520, 14.5350, 0.0000, 1.1000};
	nurbs_.p[34] = (vec4_t){28.0400, 9.2670, 0.0000, 1.0000};
	nurbs_.p[35] = (vec4_t){21.3640, 4.8300, 0.0000, 3.0000};
	nurbs_.p[36] = (vec4_t){25.7680, 15.4470, 0.0000, 5.0000};
	nurbs_.p[37] = (vec4_t){19.5390, 20.3910, 0.0000, 1.0000};
	nurbs_.p[38] = (vec4_t){17.0970, 28.5120, 0.0000, 1.0000};
	nurbs_.p[39] = (vec4_t){25.5370, 33.7500, 0.0000, 2.0000};
	nurbs_.p[40] = (vec4_t){16.6020, 30.4960, 0.0000, 1.0000};
	nurbs_.p[41] = (vec4_t){14.1990, 39.8030, 0.0000, 1.0000};
	nurbs_.p[42] = (vec4_t){8.6680, 47.4080, 0.0000, 1.0000};
	nurbs_.p[43] = (vec4_t){3.0000, 63.7940, 0.0000, 1.0000};
	nurbs_.p[44] = (vec4_t){18.4650, 67.0840, 0.0000, 1.0000};
	nurbs_.p[45] = (vec4_t){31.1970, 58.5720, 0.0000, 1.0000};
	nurbs_.p[46] = (vec4_t){39.4110, 51.3580, 0.0000, 1.0000};
	nurbs_.p[47] = (vec4_t){52.2040, 44.9710, 0.0000, 1.2000};
	nurbs_.p[48] = (vec4_t){52.9040, 49.6140, 0.0000, 1.0000};
	nurbs_.p[49] = (vec4_t){53.4780, 52.1390, 0.0000, 1.0000};
	nurbs_.p[50] = (vec4_t){54.4920, 52.1390, 0.0000, 1.0000};
	
	nurbs_.u[0]  = 0.0000;
	nurbs_.u[1]  = 0.0000;
	nurbs_.u[2]  = 0.0000;
	nurbs_.u[3]  = 0.0000;
	nurbs_.u[4]  = 0.0083;
	nurbs_.u[5]  = 0.0150;
	nurbs_.u[6]  = 0.0361;
	nurbs_.u[7]  = 0.0855;
	nurbs_.u[8]  = 0.1293;
	nurbs_.u[9]  = 0.1509;
	nurbs_.u[10] = 0.1931;
	nurbs_.u[11] = 0.2273;
	nurbs_.u[12] = 0.2435;
	nurbs_.u[13] = 0.2561;
	nurbs_.u[14] = 0.2692;
	nurbs_.u[15] = 0.2889;
	nurbs_.u[16] = 0.3170;
	nurbs_.u[17] = 0.3316;
	nurbs_.u[18] = 0.3482;
	nurbs_.u[19] = 0.3553;
	nurbs_.u[20] = 0.3649;
	nurbs_.u[21] = 0.3837;
	nurbs_.u[22] = 0.4005;
	nurbs_.u[23] = 0.4269;
	nurbs_.u[24] = 0.4510;
	nurbs_.u[25] = 0.4660;
	nurbs_.u[26] = 0.4891;
	nurbs_.u[27] = 0.5000;
	nurbs_.u[28] = 0.5109;
	nurbs_.u[29] = 0.5340;
	nurbs_.u[30] = 0.5489;
	nurbs_.u[31] = 0.5731;
	nurbs_.u[32] = 0.5994;
	nurbs_.u[33] = 0.6163;
	nurbs_.u[34] = 0.6351;
	nurbs_.u[35] = 0.6447;
	nurbs_.u[36] = 0.6518;
	nurbs_.u[37] = 0.6683;
	nurbs_.u[38] = 0.6830;
	nurbs_.u[39] = 0.7111;
	nurbs_.u[40] = 0.7307;
	nurbs_.u[41] = 0.7439;
	nurbs_.u[42] = 0.7565;
	nurbs_.u[43] = 0.7729;
	nurbs_.u[44] = 0.8069;
	nurbs_.u[45] = 0.8491;
	nurbs_.u[46] = 0.8707;
	nurbs_.u[47] = 0.9145;
	nurbs_.u[48] = 0.9639;
	nurbs_.u[49] = 0.9850;
	nurbs_.u[50] = 0.9917;
	nurbs_.u[51] = 1.0000;
	nurbs_.u[52] = 1.0000;
	nurbs_.u[53] = 1.0000;
	nurbs_.u[54] = 1.0000;

	for(i = 0; i < nurbs_.n; i++) {
		vec4_t v = nurbs_.p[i];
		fprintf(nurbspolyf, "%f, %f, %f\n", v.x, v.y, v.z); 
	}

	e3d_t linear   = eval_e3d_make(&linear_,   E3D_PPOLY1);
	e3d_t cubic0   = eval_e3d_make(&cubic0_,   E3D_PPOLY3);
	e3d_t cubic1   = eval_e3d_make(&cubic1_,   E3D_PPOLY3);
	e3d_t quintic0 = eval_e3d_make(&quintic0_, E3D_PPOLY5);
	e3d_t quintic1 = eval_e3d_make(&quintic1_, E3D_PPOLY5);
	e3d_t arc      = eval_e3d_make(&arc_,      E3D_ARC);
	e3d_t bezier   = eval_e3d_make(&bezier_,   E3D_BEZIER);
	e3d_t bspline  = eval_e3d_make(&bspline_,  E3D_BSPLINE);
	e3d_t nurbs    = eval_e3d_make(&nurbs_,    E3D_NURBS);

	for(i = 0; i < 1000; i++) {
		double t = range(0.0, 1.0, i, 1000);
		
		eval_t e = eval_e3d(&linear, t);
		fprintf(linearf, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);

		e = eval_e3d(&cubic0, t);
		fprintf(cubic0f, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);
		
		e = eval_e3d(&cubic1, t);
		fprintf(cubic1f, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);
		
		e = eval_e3d(&quintic0, t);
		fprintf(quintic0f, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);
		
		e = eval_e3d(&quintic1, t);
		fprintf(quintic1f, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);
		
		e = eval_e3d(&arc, t);
		fprintf(arcf, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);
		
		e = eval_e3d(&bezier, t);
		fprintf(bezierf, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);
		
		e = eval_e3d(&bspline, t);
		fprintf(bsplinef, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);

		e = eval_e3d(&nurbs, t);
		fprintf(nurbsf, "%f, %f, %f\n", e.p.x, e.p.y, e.p.z);

	}

	double tol = 0.5;
	double h   = 1e-3;
	
	int cubic_passed   = 0;
	int quintic_passed = 0;
	int arc_passed     = 0;
	int bezier_passed  = 0;
	int bspline_passed = 0;
	int nurbs_passed   = 0;

	for(i = 0; i < 100; i++) {
		double t = range(h, 1.0-h, i, 100);
		
		eval_t e = eval_e3d(&cubic0, t);
		eval_t c = eval_cfd(&cubic0, h, t);
		switch(eval_cmp(e, c, tol)) {
			case 1: cubic_passed = 1; break;
			case 2: cubic_passed = 2; break;
			case 3: cubic_passed = 3; break;
			default: break;
		}
		e = eval_e3d(&quintic0, t);
		c = eval_cfd(&quintic0, h, t);
		switch(eval_cmp(e, c, tol)) {
			case 1: quintic_passed = 1; break;
			case 2: quintic_passed = 2; break;
			case 3: quintic_passed = 3; break;
			default: break;
		}
		e = eval_e3d(&arc, t);
		c = eval_cfd(&arc, h, t);
		switch(eval_cmp(e, c, tol)) {
			case 1: arc_passed = 1; break;
			case 2: arc_passed = 2; break;
			case 3: arc_passed = 3; break;
			default: break;
		}
		e = eval_e3d(&bezier, t);
		c = eval_cfd(&bezier, h, t);
		switch(eval_cmp(e, c, tol)) {
			case 1: bezier_passed = 1; break;
			case 2: bezier_passed = 2; break;
			case 3: bezier_passed = 3; break;
			default: break;
		}
		e = eval_e3d(&bspline, t);
		c = eval_cfd(&bspline, h, t);
		switch(eval_cmp(e, c, tol)) {
			case 1: bspline_passed = 1; break;
			case 2: bspline_passed = 2; break;
			case 3: bspline_passed = 3; break;
			default: break;
		}
		e = eval_e3d(&nurbs, t);
		c = eval_cfd(&nurbs, h, t);
		switch(eval_cmp(e, c, tol)) {
			case 1: nurbs_passed = 1; break;
			case 2: nurbs_passed = 2; break;
			case 3: nurbs_passed = 3; break;
			default: break;
		}
	}
	
	const char* errors[] = {
		"\e[1;32msuccess\e[m", 
		"\e[0;31mposition failure\e[m",
		"\e[0;31mvelocity failure\e[m",
		"\e[0;31macceleration failure\e[m"
	};

	printf("\n=== EVALUATION === \n");

	printf("3d evaluator unit test (differentiation inequalities):\n");
	printf("linear:  \e[0;33mn/a\e[m\n");
	printf("cubic:   %s\n", errors[cubic_passed]);
	printf("quintic: %s\n", errors[quintic_passed]);
	printf("arc:     %s\n", errors[arc_passed]);
	printf("bezier:  %s\n", errors[bezier_passed]);
	printf("bspline: %s\n", errors[bspline_passed]);
	printf("nurbs:   %s\n", errors[nurbs_passed]);

	fclose(pointsf);
	fclose(linearf);
	fclose(cubic0f);
	fclose(cubic1f);
	fclose(quintic0f);
	fclose(quintic1f);
	fclose(arcf);
	fclose(arcpf);
	fclose(bezierf);
	fclose(bsplinef);
	fclose(bsplinepolyf);
	fclose(bezierpolyf);
	fclose(nurbsf);
	fclose(nurbspolyf);

	/* cba, uncomment this later when i copy alp from old repo and post about it

	int alp_n_upper = 50;
	double alpmem_ppoly3[ALP_MEMSIZE(alp_n_upper)];
	double alpmem_ppoly5[ALP_MEMSIZE(alp_n_upper)];
	double alpmem_bezier[ALP_MEMSIZE(alp_n_upper)];
	double alpmem_bspline[ALP_MEMSIZE(alp_n_upper)];
	double alpmem_nurbs[ALP_MEMSIZE(alp_n_upper)];
	
	alp_t alp_ppoly1  = alp_init(NULL, alp_n_upper);
	alp_t alp_ppoly3  = alp_init(alpmem_ppoly3, alp_n_upper);
	alp_t alp_ppoly5  = alp_init(alpmem_ppoly5, alp_n_upper);
	alp_t alp_arc     = alp_init(NULL, alp_n_upper);
	alp_t alp_bezier  = alp_init(alpmem_bezier, alp_n_upper);
	alp_t alp_bspline = alp_init(alpmem_bspline, alp_n_upper);
	alp_t alp_nurbs   = alp_init(alpmem_nurbs, alp_n_upper);

	double length_ppoly1, length_ppoly3, length_ppoly5;
	double length_arc, length_bezier, length_bspline, length_nurbs;
	
	length_ppoly1  = alp_fit(&linear,  &alp_ppoly1, alp_n_upper, 0.001);
	length_ppoly3  = alp_fit(&cubic,   &alp_ppoly3, alp_n_upper, 0.001);
	length_ppoly5  = alp_fit(&quintic, &alp_ppoly5, alp_n_upper, 0.001);
	length_arc     = alp_fit(&arc,     &alp_arc, alp_n_upper, 0.001);
	length_bezier  = alp_fit(&bezier,  &alp_bezier, alp_n_upper, 0.001);
	length_bspline = alp_fit(&bspline,   &alp_bspline, alp_n_upper, 0.001);
	length_nurbs   = alp_fit(&nurbs, &alp_nurbs, alp_n_upper, 0.001);

	printf("\n=== ARC LENGTH PARAMETERIZATION ===\n");

	printf("arc length ppoly1:  %f\n", length_ppoly1);
	printf("arc length ppoly3:  %f\n", length_ppoly3);
	printf("arc length ppoly5:  %f\n", length_ppoly5);
	printf("arc length arc:     %f\n", length_arc);
	printf("arc length bezier:  %f\n", length_bezier);
	printf("arc length bspline: %f\n", length_bspline);
	printf("arc length nurbs:   %f\n", length_nurbs);
	
	FILE* alp_ppoly1f  = fopen("csv/alp_ppoly1.csv", "w");
	FILE* alp_ppoly3f  = fopen("csv/alp_ppoly3.csv", "w");
	FILE* alp_ppoly5f  = fopen("csv/alp_ppoly5.csv", "w");
	FILE* alp_arcf     = fopen("csv/alp_arc.csv", "w");
	FILE* alp_bezierf  = fopen("csv/alp_bezier.csv", "w");
	FILE* alp_bsplinef = fopen("csv/alp_bspline.csv", "w");
	FILE* alp_nurbsf   = fopen("csv/alp_nurbs.csv", "w");

	FILE* alp_ppoly1tf  = fopen("csv/alp_ppoly1t.csv", "w");
	FILE* alp_ppoly3tf  = fopen("csv/alp_ppoly3t.csv", "w");
	FILE* alp_ppoly5tf  = fopen("csv/alp_ppoly5t.csv", "w");
	FILE* alp_arctf     = fopen("csv/alp_arct.csv", "w");
	FILE* alp_beziertf  = fopen("csv/alp_beziert.csv", "w");
	FILE* alp_bsplinetf = fopen("csv/alp_bsplinet.csv", "w");
	FILE* alp_nurbstf   = fopen("csv/alp_nurbst.csv", "w");
	
	FILE* alp_ppoly3tpf  = fopen("csv/alp_ppoly3tp.csv", "w");
	FILE* alp_ppoly5tpf  = fopen("csv/alp_ppoly5tp.csv", "w");
	FILE* alp_beziertpf  = fopen("csv/alp_beziertp.csv", "w");
	FILE* alp_bsplinetpf = fopen("csv/alp_bsplinetp.csv", "w");
	FILE* alp_nurbstpf   = fopen("csv/alp_nurbstp.csv", "w");

	for(i = 0; i < alp_ppoly3.n; i++) {
		double l = alp_ppoly3.u[i];
		double t = alp_eval(&alp_ppoly3, l);
		fprintf(alp_ppoly3tpf, "%0.30f, %0.30f\n", l, t);
	}

	for(i = 0; i < alp_ppoly5.n; i++) {
		double l = alp_ppoly5.u[i];
		double t = alp_eval(&alp_ppoly5, l);
		fprintf(alp_ppoly5tpf, "%0.30f, %0.30f\n", l, t);
	}

	for(i = 0; i < alp_bezier.n; i++) {
		double l = alp_bezier.u[i];
		double t = alp_eval(&alp_bezier, l);
		fprintf(alp_beziertpf, "%0.30f, %0.30f\n", l, t);
	}


	for(i = 0; i < alp_bspline.n; i++) {
		double l = alp_bspline.u[i];
		double t = alp_eval(&alp_bspline, l);
		fprintf(alp_bsplinetpf, "%0.30f, %0.30f\n", l, t);
	}

	for(i = 0; i < alp_nurbs.n; i++) {
		double l = alp_nurbs.u[i];
		double t = alp_eval(&alp_nurbs, l);
		fprintf(alp_nurbstpf, "%0.30f, %0.30f\n", l, t);
	}
	
	for(i = 0; i < 100000; i++) {
		double l = range(0, 1.0, i, 100000);
		
		double t_ppoly1  = alp_eval(&alp_ppoly1,  l*length_ppoly1);
		double t_ppoly3  = alp_eval(&alp_ppoly3,  l*length_ppoly3);
		double t_ppoly5  = alp_eval(&alp_ppoly5,  l*length_ppoly5);
		double t_arc     = alp_eval(&alp_arc,     l*length_arc);
		double t_bezier  = alp_eval(&alp_bezier,  l*length_bezier);
		double t_bspline = alp_eval(&alp_bspline, l*length_bspline);
		double t_nurbs   = alp_eval(&alp_nurbs,   l*length_nurbs);
		
		double v_ppoly1  = alp_eval_v(&alp_ppoly1, l*length_ppoly1);
		double v_ppoly3  = alp_eval_v(&alp_ppoly3, l*length_ppoly3);
		double v_ppoly5  = alp_eval_v(&alp_ppoly5, l*length_ppoly5);
		double v_arc     = alp_eval_v(&alp_arc, l*length_arc);
		double v_bezier  = alp_eval_v(&alp_bezier, l*length_bezier);
		double v_bspline = alp_eval_v(&alp_bspline, l*length_bspline);
		double v_nurbs   = alp_eval_v(&alp_nurbs, l*length_nurbs);

		fprintf(alp_ppoly1f,  "%0.30f, %0.30f\n", l*length_ppoly1, v_ppoly1);
		fprintf(alp_ppoly3f,  "%0.30f, %0.30f\n", l*length_ppoly3, v_ppoly3);
		fprintf(alp_ppoly5f,  "%0.30f, %0.30f\n", l*length_ppoly5, v_ppoly5);
		fprintf(alp_arcf,     "%0.30f, %0.30f\n", l*length_arc, v_arc);
		fprintf(alp_bezierf,  "%0.30f, %0.30f\n", l*length_bezier, v_bezier);
		fprintf(alp_bsplinef, "%0.30f, %0.30f\n", l*length_bspline, v_bspline);
		fprintf(alp_nurbsf,   "%0.30f, %0.30f\n", l*length_nurbs, v_nurbs);

		fprintf(alp_ppoly1tf,  "%0.30f, %0.30f\n", l*length_ppoly1, t_ppoly1);
		fprintf(alp_ppoly3tf,  "%0.30f, %0.30f\n", l*length_ppoly3, t_ppoly3);
		fprintf(alp_ppoly5tf,  "%0.30f, %0.30f\n", l*length_ppoly5, t_ppoly5);
		fprintf(alp_arctf,     "%0.30f, %0.30f\n", l*length_arc, t_arc);
		fprintf(alp_beziertf,  "%0.30f, %0.30f\n", l*length_bezier, t_bezier);
		fprintf(alp_bsplinetf, "%0.30f, %0.30f\n", l*length_bspline, t_bspline);
		fprintf(alp_nurbstf,   "%0.30f, %0.30f\n", l*length_nurbs, t_nurbs);
	}

	fclose(alp_ppoly1f);
	fclose(alp_ppoly3f);
	fclose(alp_ppoly5f);
	fclose(alp_arcf);
	fclose(alp_bsplinef);
	fclose(alp_nurbsf);
	fclose(alp_bezierf);

	fclose(alp_ppoly1tf);
	fclose(alp_ppoly3tf);
	fclose(alp_ppoly5tf);
	fclose(alp_arctf);
	fclose(alp_bsplinetf);
	fclose(alp_nurbstf);
	fclose(alp_beziertf);
	
	fclose(alp_ppoly3tpf);
	fclose(alp_ppoly5tpf);
	fclose(alp_bsplinetpf);
	fclose(alp_nurbstpf);
	
	FILE* ppoly3_curvaturef = fopen("csv/ppoly3_curvature.csv", "w");
	
	h = 1e-3;

	for(i = 0; i < 1000; i++) {
		double l = range(0.0+h, 1.0-h, i, 1000);
		
		double l0 = l - h;
		double l1 = l;
		double l2 = l + h;
		
		double t0 = alp_eval(&alp_ppoly3, l0*length_ppoly3);
		double t1 = alp_eval(&alp_ppoly3, l1*length_ppoly3);
		double t2 = alp_eval(&alp_ppoly3, l2*length_ppoly3);

		eval_t e0 = eval_e3d(&cubic, t0);
		eval_t e1 = eval_e3d(&cubic, t1);
		eval_t e2 = eval_e3d(&cubic, t2);
		
		double maxv0 = eval_maxv(&e0, 5.0);
		double maxv1 = eval_maxv(&e1, 5.0);
		double maxv2 = eval_maxv(&e2, 5.0);
		
		double d2 = (maxv2 - 2.0 * maxv1 + maxv0) / (h*h);

		fprintf(ppoly3_curvaturef, "%0.30f, %0.30f\n", l*length_ppoly3, maxv0);
	}

	fclose(ppoly3_curvaturef);
	*/

	return 0;
}
