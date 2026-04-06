#pragma once
#define _GNU_SOURCE
#include <complex.h>
#include <stdatomic.h>
#include <time.h>

#define atomic(T) _Atomic(T)
#define swap(a,b) do { auto _tmp = a; a = b; b = _tmp; } while(0)
#define for_each(array)	for (auto it = array; it < array + (sizeof(array)/sizeof(array)[0]); ++it)


static inline int min(int a, int b) { return (a < b) ? a : b; }
static inline int max(int a, int b) { return (a > b) ? a : b; }

static inline float cnormf(complex float c) {
	float r = crealf(c);
	float i = cimagf(c);
	return r*r + i*i;
}

static inline double cnorm(complex double c) {
	double r = creal(c);
	double i = cimag(c);
	return r*r + i*i;
}

double /*seconds*/ get_time_from_os() {
	struct timespec timestamp;
	clock_gettime(CLOCK_REALTIME, &timestamp);

	return timestamp.tv_sec + 1.0e-9 * timestamp.tv_nsec;
}
