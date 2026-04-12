#define _GNU_SOURCE
#include <assert.h>
#include <iso646.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "parameters.c"
#include "utility.c"

#define triangle(n) ((n)*(n+1)/2)


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
	int M[triangle(Levels)];
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


static inline sector_t sectorof(struct GT W_nu) {
	int nu1 = W_nu.M[3];
	int nu2 = W_nu.M[4];
	int nu3 = W_nu.M[5];

	static_assert(N < 256, "number of atoms does not fit into a single byte");
	return nu1 | nu2 << 8 | nu3 << 16;
}

static inline void read_sector(sector_t sector, int *nu1, int *nu2, int *nu3) {
	*nu1 =  sector        & 0xff;
	*nu2 = (sector >>  8) & 0xff;
	*nu3 = (sector >> 16) & 0xff;
}

static inline index_t indexof6(int n1, int n2, int n3, int nu1, int nu2, int nu3) {
	index_t index = (n1 - nu2)
	              + (n2 - nu3) * (nu1 - nu2 + 1)
	              + (n3 - n2)  * (nu1 - nu2 + 1)*(nu2 - nu3 + 1)
		      + 1;

	assert(index < Dimension);
	return index;
}

static inline index_t indexof(struct GT W_nu) {
	return indexof6(W_nu.M[1], W_nu.M[2], W_nu.M[0], W_nu.M[3], W_nu.M[4], W_nu.M[5]);
}


thread_local struct HamiltonianTableEntry Ground1ToExcited[Dimension][2];
thread_local struct HamiltonianTableEntry Ground2ToExcited[Dimension][2];
thread_local struct HamiltonianTableEntry ExcitedToGround1[Dimension][2];
thread_local struct HamiltonianTableEntry ExcitedToGround2[Dimension][2];
thread_local struct HamiltonianTableEntry Ground1ToGround2[Dimension][1];
thread_local struct HamiltonianTableEntry Ground2ToGround1[Dimension][1];

thread_local struct LindbladTableEntry Ground1ToExcitedJumps[Dimension][9];
thread_local struct LindbladTableEntry Ground2ToExcitedJumps[Dimension][9];
thread_local struct LindbladTableEntry ExcitedToGround1Jumps[Dimension][9];
thread_local struct LindbladTableEntry ExcitedToGround2Jumps[Dimension][9];

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
	auto index = triangle(k - 1) + j - 1;
	return W_mu.M[index] + k - j;
}

static inline double A(struct GT W_mu, int l, tau_t tau) {
	int sign = (tau[l-2] >= tau[l-1]) ? 1 : -1;
	double numerator, denominator, factor = 1;

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

static inline double zeta(struct GT W_mu, energy_level_t i, tau_t tau) {
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

	auto factor = sqrt((double)numerator / (double)denominator);

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

static inline double ratio(struct GT W_nu, struct GT W_mu) {
	return exp(log_f_factor(W_mu) - log_f_factor(W_nu));
}

static inline struct GT compute_delta(struct GT W_nu, tau_t tau, box_change_t change) {
	struct GT W_mu = W_nu;

	for (int k = 0; k < Levels; ++k) {
		if (!tau[k]) continue;

		auto index = triangle(k) + tau[k] - 1;
		W_mu.M[index] += change;
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

		table[combs[k].destination % 2].target = indexof(W_la);
		table[combs[k].destination % 2].factor += R[tau_b[Levels - 1] - 1] * zeta(W_mu, a, tau_a) * zeta(W_mu, b, tau_b);
	}
}


void lindblad_clebsh_gordan(struct GT W_nu, energy_level_t a, energy_level_t b, double R[Levels],
	                    const struct Combination combs[], struct LindbladTableEntry table[])
{
	constexpr int CombinationCount = 18;

	for (int k = 0; k < CombinationCount; ++k) {
		auto tau_a = SymmetryChanges[a][combs[k].insert_index];
		auto tau_b = SymmetryChanges[b][combs[k].remove_index];

		struct GT W_mu = compute_delta(W_nu, tau_b, RemoveBox);
		struct GT W_la = compute_delta(W_mu, tau_a, InsertBox);

		if (!check_valid_swt(W_mu)) continue;
		if (!check_valid_swt(W_la)) continue;

		auto dest = combs[k].destination;
		size_t i = dest / 2, j = dest % 2;

		table[i].sector = sectorof(W_la);
		table[i].target[j] = indexof(W_la);
		table[i].factor[j] = sqrt(R[tau_b[Levels - 1] - 1]) * zeta(W_mu, a, tau_a) * zeta(W_mu, b, tau_b);
	}
}

atomic(size_t) table_millis;

void update_tables(int nu1, int nu2, int nu3, int excitations) {
	auto start = get_time_from_os();

	constexpr struct Combination GroundToGroundCombinations[] = {
		{0,0,0}, {1,1,0}, {2,2,0}, {3,3,0}, {4,4,0}, {5,5,0}
	};

	constexpr struct Combination GroundToExcitedCombinations[] = {
		{0,0, 0}, {0,1, 1}, {1,2, 2}, {1,3, 3}, {2,4, 4}, {2,5, 5},
		{0,2, 6}, {0,3, 7}, {0,4, 8}, {0,5, 9}, {1,0,10}, {1,1,11},
		{1,4,12}, {1,5,13}, {2,0,14}, {2,1,15}, {2,2,16}, {2,3,17},
	};

	constexpr struct Combination ExcitedToGroundCombinations[] = {
		{0,0, 0}, {1,0, 1}, {2,1, 2}, {3,1, 3}, {4,2, 4}, {5,2, 5},
		{2,0, 6}, {3,0, 7}, {4,0, 8}, {5,0, 9}, {0,1,10}, {1,1,11},
		{4,1,12}, {5,1,13}, {0,2,14}, {1,2,15}, {2,2,16}, {3,2,17},
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

	struct GT W_nu = {{ [3] = nu1, [4] = nu2, [5] = nu3 }};
	double ratios[Levels] = {};

	for (int k = 0; k < Levels; ++k) {
		struct GT W_mu = W_nu;
		W_mu.M[triangle(Levels - 1) + k] -= 1;
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

				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);

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

	auto end = get_time_from_os();
	table_millis += 1000 * (end - start);
}
