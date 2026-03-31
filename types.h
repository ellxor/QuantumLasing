#pragma once
#define _GNU_SOURCE
#include <assert.h>
#include <complex.h>
#include <stdint.h>
#include <stddef.h>
#include "parameters.h"

typedef float float32;
typedef double float64;

typedef uint32_t index_t;
typedef uint32_t sector_t;

typedef enum EnergyLevel EnergyLevel;

constexpr index_t InvalidIndex = 0;

enum EnergyLevel {
	GROUND1 = 0,
	GROUND2 = 1,
	EXCITED = 2,
	LEVELS,
};

struct GT {
	int8_t M[LEVELS*(LEVELS+1)/2]; // triangle number of elements
};

constexpr int Dimension = 2000 > (N*N*N/8) ? 2000 : (N*N*N/8);
typedef complex float WaveVector[Dimension];

static inline index_t indexof3(int n1, int n2, int n3, sector_t sector) {
	int nu1 =  sector        & 0xff;
	int nu2 = (sector >>  8) & 0xff;
	int nu3 = (sector >> 16) & 0xff;

	index_t index = (n1 - nu2)
	              + (n2 - nu3) * (nu1 - nu2 + 1)
	              + (n3 - nu3) * (nu1 - nu2 + 1)*(nu2 - nu3 + 1) + 1;

	return index;
}

static inline index_t indexof(struct GT W_nu, sector_t sector) {
	return indexof3(W_nu.M[1], W_nu.M[2], W_nu.M[0], sector);
}

static inline sector_t sectorof3(int nu1, int nu2, int nu3) {
	return nu1 | nu2 << 8 | nu3 << 16;
}

static inline sector_t sectorof(struct GT W_nu) {
	return sectorof3(W_nu.M[3], W_nu.M[4], W_nu.M[5]);
}
