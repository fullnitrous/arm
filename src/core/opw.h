#include <stdio.h>
#include "vector.h"
#include "misc.h"

#ifndef _OPW_KINEMATICS_H_
#define _OPW_KINEMATICS_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct opw_arm_parameters {
	double c1, c2, c3, c4; 
	double a1, a2; 
	double b;
} opw_t;

typedef struct rot_matrix {
	double r11, r12, r13;
	double r21, r22, r23;
	double r31, r32, r33;
} rotm_t;

typedef struct denavit_hartenberg {
	double theta, d, a, alpha;
} dh_t;

typedef struct kinematic_state {
	dh_t linkage[6];
} kinstate_t;

rotm_t opw_euler(double roll, double pitch, double yaw);

kinstate_t opw_kinstate(opw_t params);
void opw_update_kinstate(kinstate_t* kinstate, double* theta);

void opw_forward(kinstate_t state, double frames[6*4*4]);
void opw_inverse(opw_t params, vec3_t pos, rotm_t rot, double sols[6*8]);

#ifdef __cplusplus
}
#endif

#endif
