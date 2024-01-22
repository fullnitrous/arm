#include "opw.h"

void __matmul4x4fast(double a[16], double b[16], double out[16]) {
	int r, c, i;
	
	for(r = 0; r < 4; r++) {
		for(c = 0; c < 4; c++) {
			out[4*r+c] = 0;
			for(i = 0; i < 4; i++) {
				out[4*r+c] += a[4*r+i]*b[4*i+c];
			}
		}
	}
	
	return;
}

rotm_t opw_euler(double roll, double pitch, double yaw) {
	rotm_t mat;
	double cos_pitch, sin_pitch, cos_roll, sin_roll, sin_yaw, cos_yaw;

	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
	cos_roll = cos(roll);
	sin_roll = sin(roll);
	sin_yaw = sin(yaw);
	cos_yaw = cos(yaw);

	mat.r11 = cos_pitch*cos_roll;
	mat.r12 = -sin_roll*cos_pitch;
	mat.r13 = sin_pitch;
	mat.r21 = sin_pitch*sin_yaw*cos_roll+sin_roll*cos_yaw;
	mat.r22 = -sin_pitch*sin_roll*sin_yaw+cos_roll*cos_yaw;
	mat.r23 = -sin_yaw*cos_pitch;
	mat.r31 = -sin_pitch*cos_roll*cos_yaw+sin_roll*sin_yaw;
	mat.r32 = sin_pitch*sin_roll*cos_yaw+sin_yaw*cos_roll;
	mat.r33 = cos_pitch*cos_yaw;

	return mat;
}


void opw_dh_transform(dh_t dh, double* m) {
	double sin_theta, cos_theta;
	double sin_alpha, cos_alpha;
	
	sin_theta = sin(dh.theta);
	cos_theta = cos(dh.theta);
	sin_alpha = sin(dh.alpha);
	cos_alpha = cos(dh.alpha); 

    m[0] = cos_theta;
	m[1] = -sin_theta*cos_alpha;
	m[2] = sin_theta*sin_alpha;
	m[3] = dh.a*cos_theta;
    
	m[4] = sin_theta;
	m[5] = cos_theta*cos_alpha;
	m[6] = -cos_theta*sin_alpha;
	m[7] = dh.a*sin_theta;

    m[8]  = 0;
	m[9]  = sin_alpha;
	m[10] = cos_alpha;
	m[11] = dh.d;
    
	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;

    return;
}

kinstate_t opw_kinstate(opw_t params) {
	kinstate_t state;

	state.linkage[0].theta = 0;
	state.linkage[0].d     = params.c1;
	state.linkage[0].a     = params.a1;
	state.linkage[0].alpha = -M_PI/2;

    state.linkage[1].theta = 0;
	state.linkage[1].d     = params.b;
	state.linkage[1].a     = params.c2;
	state.linkage[1].alpha = 0;

    state.linkage[2].theta = 0;
	state.linkage[2].d     = 0;
	state.linkage[2].a     = -params.a2;
	state.linkage[2].alpha = -M_PI/2;

    state.linkage[3].theta = 0;
	state.linkage[3].d     = params.c3;
	state.linkage[3].a     = 0;
	state.linkage[3].alpha = M_PI/2;

    state.linkage[4].theta = 0;
	state.linkage[4].d     = 0;
	state.linkage[4].a     = 0;
	state.linkage[4].alpha = -M_PI/2;
	
	state.linkage[5].theta = 0;
	state.linkage[5].d     = params.c4;
	state.linkage[5].a     = 0;
	state.linkage[5].alpha = 0;

	return state;
}

void opw_update_kinstate(kinstate_t* kinstate, double* theta) {
	kinstate->linkage[0].theta = theta[0];
	kinstate->linkage[1].theta = theta[1] - M_PI / 2;
	kinstate->linkage[2].theta = theta[2] - M_PI / 2;
	kinstate->linkage[3].theta = theta[3];
	kinstate->linkage[4].theta = theta[4];
	kinstate->linkage[5].theta = theta[5] + M_PI;

	return;
}

void opw_forward(kinstate_t state, double frames[6*4*4]) {
	int i;
	double transform_mat[16];

	opw_dh_transform(state.linkage[0], frames);
	
	for(i = 1; i < 6; i++) {
		opw_dh_transform(state.linkage[i], transform_mat);
		__matmul4x4fast(frames+16*(i-1), transform_mat, frames+16*i);
	}

    return;
}

#define c1 (param.c1)
#define c2 (param.c2)
#define c3 (param.c3)
#define c4 (param.c4)
#define a1 (param.a1)
#define a2 (param.a2)
#define b  (param.b)
#define ux0 (pos.x)
#define uy0 (pos.y)
#define uz0 (pos.z)
#define e11 (rot.r11)
#define e12 (rot.r12)
#define e13 (rot.r13)
#define e21 (rot.r21)
#define e22 (rot.r22)
#define e23 (rot.r23)
#define e31 (rot.r31)
#define e32 (rot.r32)
#define e33 (rot.r33)

double opw_theta6_correction(double theta1, rotm_t rot) {
	double x, y, ux, uy;

	ux = -sin(theta1);
	uy = cos(theta1);
	x  = uy*e33*e11 + ux*e21 + e13*e31;
	y  = -ux*e33*e11 + uy*e21 + e23*e31;

	return atan2(y, x);
}

void opw_inverse(opw_t param, vec3_t pos, rotm_t rot, double sols[6*8]) {
	double cx0, cy0, cz0;
	double nx1, s1_sqr, s2_sqr, k_sqr, s1, s2, k;
	double tmp1, tmp2, tmp3, tmp4;
	double m_i, m_ii, m_iii, m_iv;
	double theta1_i, theta1_ii;
	double theta2_i, theta2_ii, theta2_iii, theta2_iv;
	double theta3_i, theta3_ii, theta3_iii, theta3_iv;
	double theta4_i, theta4_ii, theta4_iii, theta4_iv;
	double theta4_v, theta4_vi, theta4_vii, theta4_viii;
	double theta5_i, theta5_ii, theta5_iii, theta5_iv;
	double theta5_v, theta5_vi, theta5_vii, theta5_viii;
	double theta6_i, theta6_ii, theta6_iii, theta6_iv;
	double theta6_v, theta6_vi, theta6_vii, theta6_viii;
	double _s1[4], _c1[4], _s23[4], _c23[4];
	double zero_threshold;

	zero_threshold = 1e-6;	

	cx0 = ux0 - c4*e13;
	cy0 = uy0 - c4*e23;
	cz0 = uz0 - c4*e33;

	nx1    = sqrt(cx0*cx0 + cy0*cy0 - b*b) - a1;
	s1_sqr = nx1*nx1 + (cz0-c1)*(cz0-c1);
	s2_sqr = (nx1+2*a1)*(nx1+2*a1) + (cz0-c1)*(cz0-c1);
	k_sqr  = a2*a2 + c3*c3;
	s1     = sqrt(s1_sqr);
	s2     = sqrt(s2_sqr);
	k      = sqrt(k_sqr);
	
	tmp1 = atan2(cy0, cx0);
	tmp2 = atan2(b, nx1 + a1);
	
	theta1_i  = tmp1 - tmp2;
	theta1_ii = tmp1 + tmp2 - M_PI;

	tmp1 = acos((s1_sqr+c2*c2-k_sqr) / (2*s1*c2));
	tmp2 = atan2(nx1, cz0-c1);
	tmp3 = acos((s2_sqr+c2*c2-k_sqr) / (2*s2*c2));
	tmp4 = atan2(nx1+2*a1, cz0-c1);

	theta2_i   = -tmp1 + tmp2;
	theta2_ii  =  tmp1 + tmp2;
	theta2_iii = -tmp3 - tmp4;
	theta2_iv  =  tmp3 - tmp4;

	tmp1       =  acos((s1_sqr-c2*c2-k_sqr) / (2*c2*k));
	tmp2       =  acos((s2_sqr-c2*c2-k_sqr) / (2*c2*k));
	tmp3       =  atan2(a2, c3);
	theta3_i   =  tmp1 - tmp3;
	theta3_ii  = -tmp1 - tmp3;
	theta3_iii =  tmp2 - tmp3;
	theta3_iv  = -tmp2 - tmp3;

	_s1[0] = _s1[1] = sin(theta1_i);
	_s1[2] = _s1[3] = sin(theta1_ii);

	_c1[0] = _c1[1] = cos(theta1_i);
	_c1[2] = _c1[3] = cos(theta1_ii);

	_s23[0] = sin(theta2_i   + theta3_i);
	_s23[1] = sin(theta2_ii  + theta3_ii);
	_s23[2] = sin(theta2_iii + theta3_iii);
	_s23[3] = sin(theta2_iv  + theta3_iv);

	_c23[0] = cos(theta2_i   + theta3_i);
	_c23[1] = cos(theta2_ii  + theta3_ii);
	_c23[2] = cos(theta2_iii + theta3_iii);
	_c23[3] = cos(theta2_iv  + theta3_iv);

	m_i   = e13*_s23[0]*_c1[0] + e23*_s23[0]*_s1[0] + e33*_c23[0];
	m_ii  = e13*_s23[1]*_c1[1] + e23*_s23[1]*_s1[1] + e33*_c23[1];
	m_iii = e13*_s23[2]*_c1[2] + e23*_s23[2]*_s1[2] + e33*_c23[2];
	m_iv  = e13*_s23[3]*_c1[3] + e23*_s23[3]*_s1[3] + e33*_c23[3];

	theta5_i   = atan2(sqrt(1 - m_i*m_i),     m_i);
	theta5_ii  = atan2(sqrt(1 - m_ii*m_ii),   m_ii);
	theta5_iii = atan2(sqrt(1 - m_iii*m_iii), m_iii);
	theta5_iv  = atan2(sqrt(1 - m_iv*m_iv),   m_iv);

	theta5_v    = -theta5_i;
	theta5_vi   = -theta5_ii;
	theta5_vii  = -theta5_iii;
	theta5_viii = -theta5_iv;

	if(fabs(theta5_i) < zero_threshold) {
		theta4_i = 0;
		theta6_i = opw_theta6_correction(theta1_i, rot);
	} else {
		tmp1 =  e23*_c1[0] - e13*_s1[0];
		tmp2 =  e13*_c23[0]*_c1[0] + e23*_c23[0]*_s1[0] - e33*_s23[0];
		tmp3 =  e12*_s23[0]*_c1[0] + e22*_s23[0]*_s1[0] + e32*_c23[0];
		tmp4 = -e11*_s23[0]*_c1[0] - e21*_s23[0]*_s1[0] - e31*_c23[0];
		theta4_i = atan2(tmp1, tmp2);
		theta6_i = atan2(tmp3, tmp4);
	}

	if(fabs(theta5_ii) < zero_threshold) {
		theta4_ii = 0;
		theta6_ii = opw_theta6_correction(theta1_i, rot);
	} else {
		tmp1 =  e23*_c1[1] - e13*_s1[1];
		tmp2 =  e13*_c23[1]*_c1[1] + e23*_c23[1]*_s1[1] - e33*_s23[1];
		tmp3 =  e12*_s23[1]*_c1[1] + e22*_s23[1]*_s1[1] + e32*_c23[1];
		tmp4 = -e11*_s23[1]*_c1[1] - e21*_s23[1]*_s1[1] - e31*_c23[1];
		theta4_ii = atan2(tmp1, tmp2);
		theta6_ii = atan2(tmp3, tmp4);
	}

	if(fabs(theta5_iii) < zero_threshold) {
		theta4_iii = 0;
		theta6_iii = opw_theta6_correction(theta1_ii, rot);
	} else {
		tmp1 =  e23*_c1[2] - e13*_s1[2];
		tmp2 =  e13*_c23[2]*_c1[2] + e23*_c23[2]*_s1[2] - e33*_s23[2];
		tmp3 =  e12*_s23[2]*_c1[2] + e22*_s23[2]*_s1[2] + e32*_c23[2];
		tmp4 = -e11*_s23[2]*_c1[2] - e21*_s23[2]*_s1[2] - e31*_c23[2];
		theta4_iii = atan2(tmp1, tmp2);
		theta6_iii = atan2(tmp3, tmp4);
	}

	if(fabs(theta5_iv) < zero_threshold) {
		theta4_iv = 0;
		theta6_iv = opw_theta6_correction(theta1_ii, rot);
	} else {
		tmp1 =  e23*_c1[3] - e13*_s1[3];
		tmp2 =  e13*_c23[3]*_c1[3] + e23*_c23[3]*_s1[3] - e33*_s23[3];
		tmp3 =  e12*_s23[3]*_c1[3] + e22*_s23[3]*_s1[3] + e32*_c23[3];
		tmp4 = -e11*_s23[3]*_c1[3] - e21*_s23[3]*_s1[3] - e31*_c23[3];
		theta4_iv = atan2(tmp1, tmp2);
		theta6_iv = atan2(tmp3, tmp4);
	}

	theta4_v    = theta4_i   + M_PI;
	theta4_vi   = theta4_ii  + M_PI;
	theta4_vii  = theta4_iii + M_PI;
	theta4_viii = theta4_iv  + M_PI;

	theta6_v    = theta6_i   - M_PI;
	theta6_vi   = theta6_ii  - M_PI;
	theta6_vii  = theta6_iii - M_PI;
	theta6_viii = theta6_iv  - M_PI;

	sols[6*0+0] = sols[6*1+0] = sols[6*4+0] = sols[6*5+0] = theta1_i;
	sols[6*0+1] = sols[6*4+1] = theta2_i;
	sols[6*0+2] = sols[6*4+2] = theta3_i;
	sols[6*0+3] = theta4_i;
	sols[6*0+4] = theta5_i;
	sols[6*0+5] = theta6_i;
	sols[6*1+1] = sols[6*5+1] = theta2_ii;
	sols[6*1+2] = sols[6*5+2] = theta3_ii;
	sols[6*1+3] = theta4_ii;
	sols[6*1+4] = theta5_ii;
	sols[6*1+5] = theta6_ii;
	sols[6*2+0] = sols[6*3+0] = theta1_ii;
	sols[6*2+1] = sols[6*6+1] = theta2_iii;
	sols[6*2+2] = sols[6*6+2] = theta3_iii;
	sols[6*2+3] = theta4_iii;
	sols[6*2+4] = theta5_iii;
	sols[6*2+5] = theta6_iii;
	sols[6*3+1] = sols[6*7+1] = theta2_iv;
	sols[6*3+2] = sols[6*7+2] = theta3_iv;
	sols[6*3+3] = theta4_iv;
	sols[6*3+4] = theta5_iv;
	sols[6*3+5] = theta6_iv;
	sols[6*4+3] = theta4_v;
	sols[6*4+4] = theta5_v;
	sols[6*4+5] = theta6_v;
	sols[6*5+3] = theta4_vi;
	sols[6*5+4] = theta5_vi;
	sols[6*5+5] = theta6_vi;
	sols[6*6+0] = theta1_ii;
	sols[6*7+0] = theta1_ii;
	sols[6*6+3] = theta4_vii;
	sols[6*6+4] = theta5_vii;
	sols[6*6+5] = theta6_vii;
	sols[6*7+3] = theta4_viii;
	sols[6*7+4] = theta5_viii;
	sols[6*7+5] = theta6_viii;

	return;	
}

#undef c1
#undef c2
#undef c3
#undef c4
#undef a1
#undef a2
#undef b
#undef ux0
#undef uy0
#undef uz0
#undef e11
#undef e12
#undef e13
#undef e21
#undef e22
#undef e23
#undef e31
#undef e32
#undef e33
