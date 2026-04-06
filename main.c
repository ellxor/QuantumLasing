#define _GNU_SOURCE
#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include "clebsh_gordan.c"
#include "parameters.h"
#include "random.h"
#include "utility.h"

typedef complex float wave_vector_t[Dimension];


void linear_hamiltonian_step(wave_vector_t dst, wave_vector_t src, sector_t sector, int excitations)
{
	memset(dst, 0, sizeof(wave_vector_t));

	int nu1, nu2, nu3;
	read_sector(sector, &nu1, &nu2, &nu3);

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			int excited_atoms = N - n1 - n2;
			int photon_count = excitations - excited_atoms;
			if (photon_count < 0) continue;

			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3, sector);
				dst[index] -= src[index] * TimeStep/2 * (Kappa*photon_count + GammaDown*(2 * excited_atoms) + GammaUp*(n1 + n2));

				auto coeff = src[index] * I * TimeStep;
				dst[index] -= coeff * (OmegaC * photon_count + OmegaE * excited_atoms);

				for_each(Ground1ToExcited[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count);
				for_each(Ground2ToExcited[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count);
				for_each(ExcitedToGround1[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count + 1);
				for_each(ExcitedToGround2[index]) if (it->target) dst[it->target] -= it->factor * coeff * GCoupling * sqrtf(photon_count + 1);

				constexpr auto c1 = cosf(Phi) + I*sinf(Phi);
				constexpr auto c2 = cosf(Phi) - I*sinf(Phi);

				for_each(Ground1ToGround2[index]) if (it->target) dst[it->target] -= it->factor * coeff * Omega * c1;
				for_each(Ground2ToGround1[index]) if (it->target) dst[it->target] -= it->factor * coeff * Omega * c2;
			}
		}
	}
}


void exponential_hamiltonian_step(wave_vector_t wave, sector_t sector, int excitations) {
	constexpr int RungeKuttaPoly = 4; // order of integration step

	// In this case of an exponential and linear Hamiltonian, the Runge-Kutta method
	// is identical to a Taylor series expansion, so this is performed for efficiency.
	static thread_local wave_vector_t _a, _b; // Create two temporary wave vectors as a double-buffering technique.
	auto a = wave; // Controlled by pointers which are cheap to swap.
	auto b = _b;

	int nu1, nu2, nu3;
	read_sector(sector, &nu1, &nu2, &nu3);
	int factorial = 1;

	for (int i = 1; i <= RungeKuttaPoly; ++i) {
		linear_hamiltonian_step(b, a, sector, excitations); // b now contains -i Heff dt a

		factorial *= i;
		float factor = 1.0f / factorial;

		for (int n1 = nu2; n1 <= nu1; ++n1) {
			for (int n2 = nu3; n2 <= nu2; ++n2) {
				for (int n3 = n2; n3 <= n1; ++n3) {
					index_t index = indexof3(n1, n2, n3, sector);
					wave[index] += factor * b[index];
				}
			}
		}

		if (i == 1) a = _a; // a is temporarily set to wave for first iteration to avoid a copy
		swap(a, b); // perform double-buffering
	}
}

constexpr size_t TotalJumps = 37;
typedef typeof(Ground1ToExcitedJumps) LindBladTable;

size_t choose_next_jump(float jump_table[TotalJumps])
{
	size_t index = 0;

	float accumulator = 0;
	float r = random_uniform();

	while (index < TotalJumps and (accumulator += jump_table[index]) < r) ++index;
	return index;
}


void lindblad_jump(wave_vector_t wave, LindBladTable table, size_t choice, sector_t *sector, int excitations)
{
	static thread_local wave_vector_t copy;
	memcpy(copy, wave, sizeof(wave_vector_t));
	memset(wave, 0, sizeof(wave_vector_t));

	int nu1, nu2, nu3;
	read_sector(*sector, &nu1, &nu2, &nu3);
	sector_t next_sector = 0;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3, *sector);

				wave[table[index][choice].target[0]] += copy[index] * table[index][choice].factor[0];
				wave[table[index][choice].target[1]] += copy[index] * table[index][choice].factor[1];
				next_sector |= table[index][choice].sector;
			}
		}
	}

	update_tables(next_sector, excitations);
	*sector = next_sector;
}


void lindblad_photon_annihilation(wave_vector_t wave, sector_t sector, int *excitations)
{
	int nu1, nu2, nu3;
	read_sector(sector, &nu1, &nu2, &nu3);

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			int excited_atoms = N - n1 - n2;
			int photon_count = *excitations - excited_atoms;
			photon_count = max(0, photon_count);

			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3, sector);
				wave[index] *= sqrtf(photon_count);
			}
		}
	}

	*excitations -= 1;
}


float compute_norm(wave_vector_t wave, sector_t sector) {
	int nu1, nu2, nu3;
	read_sector(sector, &nu1, &nu2, &nu3);
	float norm = 0;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof3(n1, n2, n3, sector);
				norm += cnormf(wave[index]);
			}
		}
	}

	return norm;
}

constexpr float RecordScaling = 1e6;
atomic(size_t) average_photon_count[IntegrationSteps];
atomic(size_t) average_excited_atoms[IntegrationSteps];

void run_simulation(int thread_index)
{
	set_random_seed(thread_index);
	constexpr struct GT initial = {{ (N+1)/2, N,0, N,0,0 }}; // half ground1, half ground2

	sector_t sector = sectorof(initial);
	int excitations = N - initial.M[1] - initial.M[2];

	static thread_local wave_vector_t wave;
	memset(wave, 0, sizeof(wave_vector_t));
	wave[indexof(initial, sector)] = 1;

	update_tables(sector, excitations);
	float jump_table[TotalJumps] = {0};

	for (int step = 0; step < IntegrationSteps; ++step) {
		size_t choice = choose_next_jump(jump_table);
		size_t table = choice / 9, index = choice % 9;

		switch (table) {
			case 0: lindblad_jump(wave, Ground1ToExcitedJumps, index, &sector, excitations += 1); break;
			case 1: lindblad_jump(wave, Ground2ToExcitedJumps, index, &sector, excitations += 1); break;
			case 2: lindblad_jump(wave, ExcitedToGround1Jumps, index, &sector, excitations -= 1); break;
			case 3: lindblad_jump(wave, ExcitedToGround2Jumps, index, &sector, excitations -= 1); break;
			case 4:
				if (index == 0) lindblad_photon_annihilation(wave, sector, &excitations);
				else exponential_hamiltonian_step(wave, sector, excitations);
		}

		assert(excitations >= 0);

		float inner_product = compute_norm(wave, sector);
		float scale = 1.0f / sqrtf(inner_product);

		float expected_photon_count = 0;
		float expected_excited_atoms = 0;

		int nu1, nu2, nu3;
		read_sector(sector, &nu1, &nu2, &nu3);
		memset(jump_table, 0, sizeof jump_table);

		for (int n1 = nu2; n1 <= nu1; ++n1) {
			for (int n2 = nu3; n2 <= nu2; ++n2) {
				int excited_atoms = N - n1 - n2;
				int photon_count = excitations - excited_atoms;
				if (photon_count < 0) continue;

				for (int n3 = n2; n3 <= n1; ++n3) {
					index_t index = indexof3(n1, n2, n3, sector);
					wave[index] *= scale;

					auto norm = cnormf(wave[index]);
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

		jump_table[TotalJumps - 1] = TimeStep * Kappa * expected_photon_count;

		average_photon_count[step] += (size_t)(expected_photon_count * RecordScaling);
		average_excited_atoms[step] += (size_t)(expected_excited_atoms * RecordScaling);
	}
}

atomic(int) thread_pool;
atomic(int) threads_done = 0;
atomic(size_t) millis = 0;

void *run_simulation_thread_wrapper(void *) {
	int next;

	while ((next = atomic_fetch_add(&thread_pool, -1)) > 0) {
		auto begin = get_time_from_os();
		run_simulation(next);
		auto end = get_time_from_os();

		millis += (size_t)((end - begin) * 1000);

		auto complete = atomic_fetch_add(&threads_done, 1) + 1;
		fprintf(stderr, "Trajectory [%d/%d] completed in %.3f seconds.\n", complete, TrajectoryCount, end - begin);
	}

	return nullptr;
}


int main() {
	fprintf(stderr, "Parameters: N = %d, GammaUp = %g\n", N, GammaUp);
	fprintf(stderr, "Allocating tables: %zu KB.\n", (TotalTableBytes * ThreadCount) >> 10);

	thread_pool = TrajectoryCount;
	threads_done = 0;

	pthread_t threads[ThreadCount];

	for_each(threads) pthread_create(it, nullptr, run_simulation_thread_wrapper, nullptr);
	for_each(threads) pthread_join(*it, nullptr);

	for (int i = 0; i < IntegrationSteps; ++i) {
		printf("%g\t%g\n",
			(double)average_photon_count[i]  / RecordScaling / TrajectoryCount,
			(double)average_excited_atoms[i] / RecordScaling / TrajectoryCount);
	}

	printf("# Time per trajectory: %zu ms\n", millis / TrajectoryCount);
}
