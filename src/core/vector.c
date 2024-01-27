#include "vector.h"

const int __vec3_member_offsets[] = {offsetof(struct vec3, x), offsetof(struct vec3, y), offsetof(struct vec3, z)};

vec3_t vecadd(vec3_t u, vec3_t v) {
	
	u.x += v.x;
	u.y += v.y;
	u.z += v.z;
	
	return u;
}

vec3_t vecsub(vec3_t u, vec3_t v) {
	
	u.x -= v.x;
	u.y -= v.y;
	u.z -= v.z;
	
	return u;
}

vec3_t vecmul(vec3_t u, double s) {
	
	u.x *= s; 
	u.y *= s; 
	u.z *= s; 
	
	return u;
}

vec3_t vecdiv(vec3_t u, double s) {
	
	u.x /= s; 
	u.y /= s; 
	u.z /= s; 
	
	return u;
}

vec3_t vecabs(vec3_t u) {
	
	u.x = fabs(u.x);
	u.y = fabs(u.y);
	u.z = fabs(u.z);
	
	return u;
}

vec3_t vecrotx(vec3_t u, double a) {
	vec3_t o;
	
	o.x = u.x;
	o.y = cos(a)*u.y - sin(a)*u.z;
	o.z = sin(a)*u.y + cos(a)*u.z;
	
	return o;
}

vec3_t vecroty(vec3_t u, double a) {
	vec3_t o;
	
	o.x = cos(a)*u.x + sin(a)*u.z;
	o.y = u.y;
	o.z = cos(a)*u.z - sin(a)*u.x;
	
	return o;
}

vec3_t vecrotz(vec3_t u, double a) {
	vec3_t o;
	
	o.x = cos(a)*u.x - sin(a)*u.y;
	o.y = sin(a)*u.x + cos(a)*u.y;
	o.z = u.z;
	
	return o;
}

vec3_t veccross(vec3_t u, vec3_t v) {
	vec3_t o;
	
	o.x = u.y*v.z - u.z*v.y;
	o.y = u.z*v.x - u.x*v.z;
	o.z = u.x*v.y - u.y*v.x; 
	
	return o;
}

double vecdot(vec3_t u, vec3_t v) {
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

double vecmag(vec3_t u) {
	return sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
}

vec3_t vecproj(vec3_t u, vec3_t v) {
	double tmp;
	
	tmp  = u.x*v.x + u.y*v.y + u.z*v.z;
	tmp /= v.x*v.x + v.y*v.y + v.z*v.z;
	v.x *= tmp;
	v.y *= tmp;
	v.z *= tmp;

	return v;
}

vec3_t vecdir(vec3_t u, vec3_t v) {
	vec3_t tmp0;
	double tmp1;

	tmp0.x  = u.x - v.x;
	tmp0.y  = u.y - v.y;
	tmp0.z  = u.z - v.z;
	tmp1    = sqrt(tmp0.x*tmp0.x + tmp0.y*tmp0.y + tmp0.z*tmp0.z);
	tmp0.x /= tmp1;
	tmp0.y /= tmp1;
	tmp0.z /= tmp1;

	return tmp0;
}

vec3_t vecnorm(vec3_t u) {
	double tmp;

	tmp = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
	
	u.x /= tmp;
	u.y /= tmp;
	u.z /= tmp;
	
	return u;
}

double vecangle(vec3_t u, vec3_t v) {
	double tmp;

	tmp  = u.x*v.x + u.y*v.y + u.z*v.z;
	tmp /= sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
	tmp /= sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

	return acos(tmp);
}

vec3_t veczeros(void) {
	vec3_t o;
	
	o.x = 0.0;
	o.y = 0.0;
	o.z = 0.0;
	
	return o;
}

vec3_t vecrotn(vec3_t v, vec3_t n, double a) {
	vec3_t o;
	double tmp;

	tmp  = cos(a);
	o.x  = v.x*tmp;
	o.y  = v.y*tmp;
	o.z  = v.z*tmp;
	tmp  = (1.0 - tmp)*(v.x*n.x + v.y*n.y + v.z*n.z);
	o.x += tmp*n.x;
	o.y += tmp*n.y;
	o.z += tmp*n.z;
	tmp  = sin(a);
	o.x += tmp*(n.y*v.z - n.z*v.y);
	o.y += tmp*(n.z*v.x - n.x*v.z);
	o.z += tmp*(n.x*v.y - n.y*v.x);

	return o;
}

int veccmp(vec3_t tst, vec3_t ref, double tol) {
	int i;
	double tv;
	double rv;

	for(i = 0; i < 3; i++) {
		tv = *VEC3GET(&tst, i);
		rv = *VEC3GET(&ref, i);		
		if((rv+tol) < tv || tv < (rv-tol)) {
			return 1;
		}
	}

	return 0;
}
