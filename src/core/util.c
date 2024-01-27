#include "util.h"

double range(double start, double stop, int i, int res) {
    double t = i/(double)(res-1);
    double v = start + t*(stop-start);
    return v;
}

double min(double a, double b) {
	return a < b ? a : b;
}

double max(double a, double b) {
	return a > b ? a : b;
}
