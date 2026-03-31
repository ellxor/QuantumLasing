#include <stdatomic.h>
#define _GNU_SOURCE
// #include <fenv.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <threads.h>

#include "clebsh_gordan.c"
#include "parameters.h"
#include "random.h"
#include "types.h"
#include "utility.h"


void linear_hamiltonian_step(WaveVector dst, WaveVector src, sector_t sector, int excitations)
{
	memset(dst, 0, sizeof(WaveVector));

	int nu1 =  sector        & 0xff;
	int nu2 = (sector >>  8) & 0xff;
	int nu3 = (sector >> 16) & 0xff;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			int excited_atoms = N - n1 - n2;
			int photon_count = excitations - excited_atoms;
			if (photon_count < 0) continue; // OPTIMISATION: this state cannot not occupied: src[index] = 0

			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3);
				dst[index] -= src[index] * TimeStep/2 * (Kappa*photon_count + GammaDown*(2 * excited_atoms) + GammaUp*(n1 + n2));

				auto coeff = src[index] * I * TimeStep;
				dst[index] -= coeff * (OmegaC * photon_count + OmegaE * excited_atoms);

				for_each(Ground1ToExcited[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count);
				for_each(Ground2ToExcited[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count);
				for_each(ExcitedToGround1[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count + 1);
				for_each(ExcitedToGround2[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count + 1);

				static_assert(Phi == 0, "Microwave-drive phase-shift is not yet implemented.");
				for_each(Ground1ToGround2[index]) if (it->target) dst[it->target] -= it->factor * coeff * Omega;
				for_each(Ground2ToGround1[index]) if (it->target) dst[it->target] -= it->factor * coeff * Omega;
			}
		}
	}
}


void exponential_hamiltonian_step(WaveVector wave, sector_t sector, int excitations) {
	constexpr int RungeKuttaPoly = 4; // order of integration step

	// In this case of an exponential and linear Hamiltonian, the Runge-Kutta method
	// is identical to a Taylor series expansion, so this is performed for efficiency.
	static thread_local WaveVector _a, _b; // Create two temporary wave vectors as a double-buffering technique.
	auto a = wave; // Controlled by pointers which are cheap to swap.
	auto b = _b;

	int nu1 =  sector        & 0xff;
	int nu2 = (sector >>  8) & 0xff;
	int nu3 = (sector >> 16) & 0xff;

	int factorial = 1;

	for (int i = 1; i <= RungeKuttaPoly; ++i) {
		linear_hamiltonian_step(b, a, sector, excitations); // b now contains -i Heff dt a

		factorial *= i;
		float factor = 1.0f / factorial;

		for (int n1 = nu2; n1 <= nu1; ++n1) {
			for (int n2 = nu3; n2 <= nu2; ++n2) {
				for (int n3 = n2; n3 <= n1; ++n3) {
					index_t index = indexof3(n1, n2, n3);
					wave[index] += factor * b[index];
				}
			}
		}

		if (i == 1) a = _a; // a is temporarily set to wave for first iteration to avoid a copy
		swap(a, b); // perform double-buffering
	}
}

constexpr size_t TotalJumps = 33;
typedef typeof(Ground1ToExcitedJumps) LindBladTable;

size_t choose_next_jump(float jump_table[TotalJumps])
{
	size_t index = 0;

	float accumulator = 0;
	float r = random_uniform();

	while (index < TotalJumps and (accumulator += jump_table[index]) < r) ++index;
	return index;
}


void lindblad_jump(WaveVector wave, LindBladTable table, size_t choice, sector_t *sector)
{
	static thread_local WaveVector copy;
	memcpy(copy, wave, sizeof(WaveVector));
	memset(wave, 0, sizeof(WaveVector));

	int nu1 =  *sector        & 0xff;
	int nu2 = (*sector >>  8) & 0xff;
	int nu3 = (*sector >> 16) & 0xff;

	sector_t next_sector = 0;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3);

				if (table[index][choice].target[0]) wave[table[index][choice].target[0]] += copy[index] * table[index][choice].factor[0];
				if (table[index][choice].target[1]) wave[table[index][choice].target[1]] += copy[index] * table[index][choice].factor[1];
				next_sector |= table[index][choice].sector;
			}
		}
	}

	if (next_sector != *sector) {
		int nu1 = next_sector & 0xff;
		int nu2 = (next_sector >> 8) & 0xff;
		int nu3 = (next_sector >> 16) & 0xff;

		update_tables(nu1, nu2, nu3);
		*sector = next_sector;
	}
}


void lindblad_photon_annihilation(WaveVector wave, sector_t sector, int *excitations)
{
	int nu1 =  sector        & 0xff;
	int nu2 = (sector >>  8) & 0xff;
	int nu3 = (sector >> 16) & 0xff;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			int excited_atoms = N - n1 - n2;
			int photon_count = *excitations - excited_atoms;
			if (photon_count < 0) continue;

			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3);
				wave[index] = sqrtf(photon_count);
			}
		}
	}

	*excitations -= 1;
}


float compute_norm(WaveVector wave, sector_t sector) {
	int nu1 = sector & 0xff;
	int nu2 = (sector >> 8) & 0xff;
	int nu3 = (sector >> 16) & 0xff;

	float norm = 0;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3);
				norm += cnormf(wave[index]);
			}
		}
	}

	return norm;
}


atomic(float) average_photon_count[IntegrationSteps];
atomic(float) average_excited_atoms[IntegrationSteps];

void run_simulation(int thread_index)
{
	set_random_seed(thread_index);
	struct GT initial = {{ N/2, N/2,0, N,0,0 }}; // all atoms in first ground state

	sector_t sector = sectorof(initial);
	int excitations = N - initial.M[1] - initial.M[2];

	int nu1 = initial.M[3];
	int nu2 = initial.M[4];
	int nu3 = initial.M[5];

	static thread_local WaveVector wave = {};
	wave[indexof(initial)] = 1;

	float jump_table[TotalJumps] = {0};
	update_tables(nu1, nu2, nu3);

	for (int step = 0; step < IntegrationSteps; ++step) {
		size_t choice = choose_next_jump(jump_table);
		size_t table = choice / 8, index = choice % 8;

		switch (table) {
			case 0: lindblad_jump(wave, Ground1ToExcitedJumps, index, &sector); excitations += 1; break;
			case 1: lindblad_jump(wave, Ground2ToExcitedJumps, index, &sector); excitations += 1; break;
			case 2: lindblad_jump(wave, ExcitedToGround1Jumps, index, &sector); excitations = max(0, excitations - 1); break;
			case 3: lindblad_jump(wave, ExcitedToGround2Jumps, index, &sector); excitations = max(0, excitations - 1); break;
			case 4:
				if (index == 0) lindblad_photon_annihilation(wave, sector, &excitations);
				else exponential_hamiltonian_step(wave, sector, excitations);
		}

		float inner_product = compute_norm(wave, sector);
		float scale = 1.0f / sqrtf(inner_product);

		float expected_photon_count = 0;
		float expected_excited_atoms = 0;

		int nu1 = sector & 0xff;
		int nu2 = (sector >> 8) & 0xff;
		int nu3 = (sector >> 16) & 0xff;

		memset(jump_table, 0, sizeof jump_table);

		for (int n1 = nu2; n1 <= nu1; ++n1) {
			for (int n2 = nu3; n2 <= nu2; ++n2) {
				for (int n3 = n2; n3 <= n1; ++n3) {
					index_t index = indexof3(n1, n2, n3);
					wave[index] *= scale;

					auto norm = cnormf(wave[index]);

					int excited_atoms = N - n1 - n2;
					int photon_count = max(0, excitations - excited_atoms);
					expected_photon_count += photon_count * norm;
					expected_excited_atoms += excited_atoms * norm;

					size_t jump_index = 0;
					for_each(Ground1ToExcitedJumps[index]) jump_table[jump_index++] += TimeStep * GammaUp   * norm * (powf(it->factor[0], 2) + powf(it->factor[1], 2));
					for_each(Ground2ToExcitedJumps[index]) jump_table[jump_index++] += TimeStep * GammaUp   * norm * (powf(it->factor[0], 2) + powf(it->factor[1], 2));
					for_each(ExcitedToGround1Jumps[index]) jump_table[jump_index++] += TimeStep * GammaDown * norm * (powf(it->factor[0], 2) + powf(it->factor[1], 2));
					for_each(ExcitedToGround2Jumps[index]) jump_table[jump_index++] += TimeStep * GammaDown * norm * (powf(it->factor[0], 2) + powf(it->factor[1], 2));
				}
			}
		}

		average_photon_count[step] += expected_photon_count;
		average_excited_atoms[step] += expected_excited_atoms;
		jump_table[TotalJumps - 1] = TimeStep * Kappa * expected_photon_count;
	}

	// printf("# Tables: %zu KiB\n", TotalTableBytes >> 10);
	// printf("# Table Generation: %f ms\n", 1000.0 * (t2 - t1) / CLOCKS_PER_SEC);
	// printf("# Simulation: %f ms\n", 1000.0 * (t3 - t2) / CLOCKS_PER_SEC);
}

atomic(int) thread_pool;
atomic(int) threads_done;

int run_simulation_thread_wrapper(void *) {
	int next;

	while ((next = atomic_fetch_add(&thread_pool, -1)) > 0) {
		auto begin = get_time_from_os();
		run_simulation(next);
		auto end = get_time_from_os();

		auto complete = atomic_fetch_add(&threads_done, 1) + 1;
		fprintf(stderr, "Thread [%d/%d] completed in %.3f seconds.\n", complete, TrajectoryCount, end - begin);
	}

	return 0;
}


int main() {
	thread_pool = TrajectoryCount;
	threads_done = 0;

	thrd_t threads[ThreadCount];
	size_t index = 0;

	for_each(threads) thrd_create(it, run_simulation_thread_wrapper, (void*)index++);
	for_each(threads) thrd_join(*it, nullptr);

	for (int i = 0; i < IntegrationSteps; ++i) {
		printf("%g\t%g\n",
			average_photon_count[i] / TrajectoryCount,
			average_excited_atoms[i] / TrajectoryCount);
	}
}
