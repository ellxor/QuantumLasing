#define _GNU_SOURCE
#include <assert.h>
#include <iso646.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "parameters.h"
#include "utility.h"


constexpr size_t Dimension = 2048 > (N*N*N/8) ? 2048 : (N*N*N/8);
constexpr size_t Levels = 3; // three-level-system

typedef enum EnergyLevel energy_level_t;
typedef enum BoxChange box_change_t;
typedef const int tau_t[Levels];

typedef uint32_t index_t;
typedef uint32_t sector_t;

enum EnergyLevel {
	Ground1,
	Ground2,
	Excited,
};

enum BoxChange {
	RemoveBox = -1,
	InsertBox = +1,
};

struct GT {
	int M[Levels*(Levels+1)/2];
};

struct HamiltonianTableEntry {
	index_t target;
	float   factor;
};

struct LindbladTableEntry {
	sector_t sector;
	index_t  target[2];
	float    factor[2];
};

struct Combination {
	uint8_t insert_index;
	uint8_t remove_index;
	uint8_t destination;
};



static inline sector_t sectorof3(int nu1, int nu2, int nu3) {
	static_assert(N < 256, "number of atoms does not fit into a single byte");
	return nu1 | nu2 << 8 | nu3 << 16;
}

static inline sector_t sectorof(struct GT W_nu) {
	return sectorof3(W_nu.M[3], W_nu.M[4], W_nu.M[5]);
}

static inline void read_sector(sector_t sector, int *nu1, int *nu2, int *nu3) {
	*nu1 =  sector        & 0xff;
	*nu2 = (sector >>  8) & 0xff;
	*nu3 = (sector >> 16) & 0xff;
}


static inline index_t indexof3(int n1, int n2, int n3, sector_t sector) {
	int nu1, nu2, nu3;
	read_sector(sector, &nu1, &nu2, &nu3);

	index_t index = (n1 - nu2)
	              + (n2 - nu3) * (nu1 - nu2 + 1)
	              + (n3 - n2)  * (nu1 - nu2 + 1)*(nu2 - nu3 + 1) + 1;

	assert(index < Dimension);
	return index;
}

static inline index_t indexof(struct GT W_nu, sector_t sector) {
	return indexof3(W_nu.M[1], W_nu.M[2], W_nu.M[0], sector);
}


thread_local struct HamiltonianTableEntry Ground1ToExcited[Dimension][2];
thread_local struct HamiltonianTableEntry Ground2ToExcited[Dimension][2];
thread_local struct HamiltonianTableEntry ExcitedToGround1[Dimension][2];
thread_local struct HamiltonianTableEntry ExcitedToGround2[Dimension][2];
thread_local struct HamiltonianTableEntry Ground1ToGround2[Dimension][1];
thread_local struct HamiltonianTableEntry Ground2ToGround1[Dimension][1];

thread_local struct LindbladTableEntry Ground1ToExcitedJumps[Dimension][8];
thread_local struct LindbladTableEntry Ground2ToExcitedJumps[Dimension][8];
thread_local struct LindbladTableEntry ExcitedToGround1Jumps[Dimension][8];
thread_local struct LindbladTableEntry ExcitedToGround2Jumps[Dimension][8];

constexpr size_t TotalTableBytes
	= sizeof Ground1ToExcited + sizeof Ground2ToExcited
        + sizeof ExcitedToGround1 + sizeof ExcitedToGround2
        + sizeof Ground1ToGround2 + sizeof Ground2ToGround1
        + sizeof Ground1ToExcitedJumps + sizeof Ground2ToExcitedJumps
        + sizeof ExcitedToGround1Jumps + sizeof ExcitedToGround2Jumps;


constexpr tau_t SymmetryChanges[Levels][6] = {
	[Ground1] = {{1,1,1}, {1,2,1}, {1,1,2}, {1,2,2}, {1,1,3}, {1,2,3}},
	[Ground2] = {{0,1,1}, {0,2,1}, {0,1,2}, {0,2,2}, {0,1,3}, {0,2,3}},
	[Excited] = {{0,0,1}, {0,0,2}, {0,0,3}},
};


static inline int p(struct GT W_mu, int j, int k) {
	auto index = k*(k-1)/2 + j - 1; // (k-1) triangle number + (j-1)
	return W_mu.M[index] + k - j;
}

static inline float A(struct GT W_mu, int l, tau_t tau) {
	int sign = (tau[l-2] >= tau[l-1]) ? 1 : -1;
	float numerator, denominator, factor = 1;

	for (int k = 1; k <= l; ++k) {
		if (k != tau[l-1]) {
			numerator = p(W_mu, tau[l-2], l-1) - p(W_mu, k, l) + 1;
			denominator = p(W_mu, tau[l-1], l) - p(W_mu, k, l);
			factor *= numerator / denominator;
		}
	}

	for (int k = 1; k <= l - 1; ++k) {
		if (k != tau[l-2]) {
			numerator = p(W_mu, tau[l-1], l) - p(W_mu, k, l-1);
			denominator = p(W_mu, tau[l-2],l-1) - p(W_mu, k, l-1) + 1;
			factor *= numerator / denominator;
		}
	}

	return sign * sqrt(factor);
}

static inline float zeta(struct GT W_mu, energy_level_t i, tau_t tau) {
	int numerator = 1;
	int denominator = 1;

	for (int k = 1; k <= i; ++k) {
		numerator *= p(W_mu, tau[i], i+1) - p(W_mu, k, i);
	}

	for (int k = 1; k <= i + 1; ++k) {
		if (k != tau[i]) {
			denominator *= p(W_mu, tau[i], i+1) - p(W_mu, k, i+1);
		}
	}

	float factor = sqrtf((float)numerator / (float)denominator);

	for (int l = i+2; l <= Levels; ++l) {
		factor *= A(W_mu, l, tau);
	}

	return factor;
}


static inline double log_f_factor(struct GT W_nu) {
	int nu1 = W_nu.M[3];
	int nu2 = W_nu.M[4];
	int nu3 = W_nu.M[5];

	double log_numerator = log((nu1 - nu2 + 1)*(nu2 - nu3 + 1)*(nu1 - nu3 + 2));
	double log_denominator = lgamma(nu1 + 3) + lgamma(nu2 + 2) + lgamma(nu3 + 1);

	return log_numerator - log_denominator;
}

static inline float ratio(struct GT W_nu, struct GT W_mu) {
	return exp(log_f_factor(W_mu) - log_f_factor(W_nu));
}

static inline struct GT compute_delta(struct GT W_nu, tau_t tau, box_change_t change) {
	struct GT W_mu = W_nu;

	for (int k = 0; k < Levels; ++k) {
		auto index = k*(k+1)/2 + tau[k] - 1;
		if (index >= 0) W_mu.M[index] += change;
	}

	return W_mu;
}

static inline bool check_valid_swt(struct GT W_nu) {
	int n1 = W_nu.M[1];
	int n2 = W_nu.M[2];
	int n3 = W_nu.M[0];
	int nu1 = W_nu.M[3];
	int nu2 = W_nu.M[4];
	int nu3 = W_nu.M[5];

	return (nu1 >= n1) and (n1 >= nu2) and (nu2 >= n2)
	   and (n1 >= n3) and (n3 >= n2) and (n2 >= nu3) and (nu3 >= 0);
}


void hamiltonian_clebsh_gordan(struct GT W_nu, energy_level_t a, energy_level_t b, double R[Levels],
                               const struct Combination combs[], struct HamiltonianTableEntry table[])
{
	constexpr int CombinationCount = 6;

	for (int k = 0; k < CombinationCount; ++k) {
		auto tau_a = SymmetryChanges[a][combs[k].insert_index];
		auto tau_b = SymmetryChanges[b][combs[k].remove_index];

		struct GT W_mu = compute_delta(W_nu, tau_b, RemoveBox);
		struct GT W_la = compute_delta(W_mu, tau_a, InsertBox);

		if (!check_valid_swt(W_mu)) continue;
		if (!check_valid_swt(W_la)) continue;

		table[combs[k].destination].target = indexof(W_la, sectorof(W_nu));
		table[combs[k].destination].factor += R[tau_b[Levels - 1] - 1] * zeta(W_mu, a, tau_a) * zeta(W_mu, b, tau_b);
	}
}


void lindblad_clebsh_gordan(struct GT W_nu, energy_level_t a, energy_level_t b, double R[Levels],
	                    const struct Combination combs[], struct LindbladTableEntry table[])
{
	constexpr int CombinationCount = 18;
	constexpr int UniqueSectors = 7;

	float ratios[CombinationCount] = {};
	float factors[CombinationCount] = {};

	for (int k = 0; k < CombinationCount; ++k) {
		auto tau_a = SymmetryChanges[a][combs[k].insert_index];
		auto tau_b = SymmetryChanges[b][combs[k].remove_index];

		struct GT W_mu = compute_delta(W_nu, tau_b, RemoveBox);
		struct GT W_la = compute_delta(W_mu, tau_a, InsertBox);

		if (!check_valid_swt(W_mu)) continue;
		if (!check_valid_swt(W_la)) continue;

		factors[k] = zeta(W_mu, a, tau_a) * zeta(W_mu, b, tau_b);
		ratios[k] = R[tau_b[Levels - 1] - 1];

		auto dest = combs[k].destination;
		size_t i = dest / 2, j = dest % 2;

		table[i].sector = sectorof(W_la);
		table[i].target[j] = indexof(W_la, table[i].sector);
		table[i].factor[j] += ratios[k] * pow(factors[k], 2);
	}

	// Solve for same symmetry sector case
	float first  = table[0].factor[0];
	float second = table[0].factor[1];
	float cross_term = 0;

	for (int k = 0; k < Levels; ++k) {
		cross_term += ratios[k] * factors[k] * factors[k+Levels];
	}

	if (cross_term == 0) {
		table[0].factor[0] = sqrt(first);
		table[0].factor[1] = sqrt(second);
	}

	else {
		// This condition MUST hold for the system to be solved with 2 jumps
		assert(pow(cross_term,2) <= table[0].factor[0]*table[0].factor[1]);

		table[0].factor[0] = sqrt(first);
		table[0].factor[1] = 0;
		table[UniqueSectors].factor[0] = cross_term / sqrt(first);
		table[UniqueSectors].factor[1] = sqrt(second - pow(cross_term,2)/first);

		table[UniqueSectors].target[0] = table[0].target[0];
		table[UniqueSectors].target[1] = table[0].target[1];
		table[UniqueSectors].sector = table[0].sector;
	}

	// Cross term signs for rest
	uint8_t sign_flips = 0;

	for (int k = 2*Levels; k < CombinationCount; k += 2) {
		auto cross_term = factors[k] * factors[k + 1];
		if (cross_term < 0) sign_flips |= 1 << (combs[k].destination / 2);
	}

	for (int i = 1; i < UniqueSectors; ++i) {
		table[i].factor[0] = sqrt(table[i].factor[0]);
		table[i].factor[1] = sqrt(table[i].factor[1]);
		if (sign_flips & (1 << i)) table[i].factor[1] *= -1.0f;
	}
}


void update_tables(sector_t sector, int excitations) {
	constexpr struct Combination GroundToGroundCombinations[] = {
		{0,0,0}, {1,1,0}, {2,2,0}, {3,3,0}, {4,4,0}, {5,5,0}
	};

	constexpr struct Combination GroundToExcitedCombinations[] = {
		{0,0,0}, {1,2,0}, {2,4,0}, {0,1,1}, {1,3,1}, {2,5,1},
		{0,2, 2}, {0,3, 3}, {0,4, 4}, {0,5, 5},
		{1,0, 6}, {1,1, 7}, {1,4, 8}, {1,5, 9},
		{2,0,10}, {2,1,11}, {2,2,12}, {2,3,13},
	};

	constexpr struct Combination ExcitedToGroundCombinations[18] = {
		{0,0,0}, {2,1,0}, {4,2,0}, {1,0,1}, {3,1,1}, {5,2,1},
		{2,0, 2}, {3,0, 3}, {4,0, 4}, {5,0, 5},
		{0,1, 6}, {1,1, 7}, {4,1, 8}, {5,1, 9},
		{0,2,10}, {1,2,11}, {2,2,12}, {3,2,13},
	};

	memset(Ground1ToExcited, 0, sizeof Ground1ToExcited);
	memset(Ground2ToExcited, 0, sizeof Ground2ToExcited);
	memset(ExcitedToGround1, 0, sizeof ExcitedToGround1);
	memset(ExcitedToGround2, 0, sizeof ExcitedToGround2);
	memset(Ground1ToGround2, 0, sizeof Ground1ToGround2);
	memset(Ground2ToGround1, 0, sizeof Ground2ToGround1);

	memset(Ground1ToExcitedJumps, 0, sizeof Ground1ToExcitedJumps);
	memset(Ground2ToExcitedJumps, 0, sizeof Ground2ToExcitedJumps);
	memset(ExcitedToGround1Jumps, 0, sizeof ExcitedToGround1Jumps);
	memset(ExcitedToGround2Jumps, 0, sizeof ExcitedToGround2Jumps);

	int nu1, nu2, nu3;
	read_sector(sector, &nu1, &nu2, &nu3);
	struct GT W_nu = {{ [3] = nu1, [4] = nu2, [5] = nu3 }};

	double ratios[Levels] = {};

	for (int k = 0; k < Levels; ++k) {
		constexpr size_t offset = Levels*(Levels-1)/2;
		struct GT W_mu = W_nu;

		W_mu.M[offset + k] -= 1;
		ratios[k] = ratio(W_nu, W_mu);
	}

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			int excited_atoms = N - n1 - n2;
			int photon_count = excitations - excited_atoms;
			if (photon_count < 0) continue; // OPTIMISATION: this state cannot not occupied: src[index] = 0

			for (int n3 = n2; n3 <= n1; ++n3) {
				W_nu.M[0] = n3;
				W_nu.M[1] = n1;
				W_nu.M[2] = n2;

				index_t index = indexof3(n1, n2, n3, sector);

				hamiltonian_clebsh_gordan(W_nu, Excited, Ground1, ratios, GroundToExcitedCombinations, Ground1ToExcited[index]);
				hamiltonian_clebsh_gordan(W_nu, Excited, Ground2, ratios, GroundToExcitedCombinations, Ground2ToExcited[index]);
				hamiltonian_clebsh_gordan(W_nu, Ground1, Excited, ratios, ExcitedToGroundCombinations, ExcitedToGround1[index]);
				hamiltonian_clebsh_gordan(W_nu, Ground2, Excited, ratios, ExcitedToGroundCombinations, ExcitedToGround1[index]);
				hamiltonian_clebsh_gordan(W_nu, Ground1, Ground1, ratios, GroundToGroundCombinations,  Ground1ToGround2[index]);
				hamiltonian_clebsh_gordan(W_nu, Ground1, Ground2, ratios, GroundToGroundCombinations,  Ground2ToGround1[index]);

				lindblad_clebsh_gordan(W_nu, Excited, Ground1, ratios, GroundToExcitedCombinations, Ground1ToExcitedJumps[index]);
				lindblad_clebsh_gordan(W_nu, Excited, Ground2, ratios, GroundToExcitedCombinations, Ground2ToExcitedJumps[index]);
				lindblad_clebsh_gordan(W_nu, Ground1, Excited, ratios, ExcitedToGroundCombinations, ExcitedToGround1Jumps[index]);
				lindblad_clebsh_gordan(W_nu, Ground2, Excited, ratios, ExcitedToGroundCombinations, ExcitedToGround2Jumps[index]);
			}
		}
	}
}
