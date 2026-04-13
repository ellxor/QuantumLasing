#define _GNU_SOURCE
#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include "clebsh_gordan.c"
#include "parameters.c"
#include "random.c"
#include "utility.c"

typedef complex float wave_vector_t[Dimension];


void effective_hamiltonian_step(wave_vector_t dst, int nu1, int nu2, int nu3, int excitations)
{
	alignas(64) static thread_local wave_vector_t src;
	memcpy(src, dst, sizeof src);

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			int excited_atoms = N - n1 - n2;
			int photon_count = excitations - excited_atoms;
			if (photon_count < 0) continue;

			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
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


atomic(size_t) jump_count = 0;

void lindblad_jump(wave_vector_t wave, LindBladTable table, size_t choice, int *nu1, int *nu2, int *nu3, int excitations)
{
	++jump_count;

	alignas(64) static thread_local wave_vector_t copy;
	memcpy(copy, wave, sizeof(wave_vector_t));
	memset(wave, 0, sizeof(wave_vector_t));

	sector_t next_sector = 0;

	for (int n1 = *nu2; n1 <= *nu1; ++n1) {
		for (int n2 = *nu3; n2 <= *nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, *nu1, *nu2, *nu3);

				wave[table[index][choice].target[0]] += copy[index] * table[index][choice].factor[0];
				wave[table[index][choice].target[1]] += copy[index] * table[index][choice].factor[1];
				next_sector |= table[index][choice].sector;
			}
		}
	}

	read_sector(next_sector, nu1, nu2, nu3);
	update_tables(*nu1, *nu2, *nu3, excitations, 0, 0);
}


void lindblad_photon_annihilation(wave_vector_t wave, int nu1, int nu2, int nu3, int excitations)
{
	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			int excited_atoms = N - n1 - n2;
			int photon_count = max(0, excitations - excited_atoms + 1);

			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
				wave[index] *= sqrtf(photon_count);
			}
		}
	}
}


float compute_norm(wave_vector_t wave, int nu1, int nu2, int nu3) {
	float norm = 0;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
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

	int nu1 = initial.M[3];
	int nu2 = initial.M[4];
	int nu3 = initial.M[5];
	int excitations = N - initial.M[1] - initial.M[2];

	alignas(64) static thread_local wave_vector_t wave;
	memset(wave, 0, sizeof(wave_vector_t));
	wave[indexof(initial)] = 1;

	update_tables(nu1, nu2, nu3, excitations, 0, 0);
	float jump_table[5] = {0};

	for (int step = 0; step < IntegrationSteps; ++step) {
		float r = random_uniform();
		float t = 0; for_each(jump_table) t += *it;

		if (likely(r >= t)) {
			effective_hamiltonian_step(wave, nu1, nu2, nu3, excitations);
		}

		else if (r >= t - jump_table[4]) {
			lindblad_photon_annihilation(wave, nu1, nu2, nu3, excitations -= 1);
		}

		else {
			size_t choice = 0;
			float accumulator = 0;

			for_each(jump_table) {
				if ((accumulator + *it) >= r) break;
				accumulator += *it;
				++choice;
			}

			energy_level_t as[] = { Excited, Excited, Ground1, Ground2 };
			energy_level_t bs[] = { Ground1, Ground2, Excited, Excited };

			float factors[] = { GammaUp, GammaUp, GammaDown, GammaDown };
			float factor = TimeStep * factors[choice];

			int excitation_diffs[] = { +1, +1, -1, -1 };
			int excitation_diff = excitation_diffs[choice];

			auto a = as[choice];
			auto b = bs[choice];

			update_tables(nu1, nu2, nu3, excitations, a, b);
			float jump_sub_table[9] = {0};

			for (int n1 = nu2; n1 <= nu1; ++n1) {
				for (int n2 = nu3; n2 <= nu2; ++n2) {
					int excited_atoms = N - n1 - n2;
					int photon_count = excitations - excited_atoms;
					if (photon_count < 0) continue;

					for (int n3 = n2; n3 <= n1; ++n3) {
						index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
						auto norm = cnormf(wave[index]);

						size_t jump_index = 0;
						for_each(LindbladJumps[index]) jump_sub_table[jump_index++] += factor * norm * (powf(it->factor[0], 2) + powf(it->factor[1], 2));
					}
				}
			}

			r -= accumulator;
			choice = 0, accumulator = 0;

			for_each(jump_sub_table) {
				if ((accumulator + *it) >= r) break;
				accumulator += *it;
				++choice;
			}

			assert(choice < 9);
			lindblad_jump(wave, LindbladJumps, choice, &nu1, &nu2, &nu3, excitations += excitation_diff);
		}

		assert(excitations >= 0);

		float inner_product = compute_norm(wave, nu1, nu2, nu3);
		float scale = 1.0f / sqrtf(inner_product);

		float expected_photon_count = 0;
		float expected_excited_atoms = 0;
		float expected_ground1_atoms = 0;
		float expected_ground2_atoms = 0;

		for (int n1 = nu2; n1 <= nu1; ++n1) {
			for (int n2 = nu3; n2 <= nu2; ++n2) {
				int excited_atoms = N - n1 - n2;
				int photon_count = excitations - excited_atoms;
				if (photon_count < 0) continue;

				for (int n3 = n2; n3 <= n1; ++n3) {
					index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
					wave[index] *= scale;

					auto norm = cnormf(wave[index]);
					expected_photon_count += photon_count * norm;
					expected_excited_atoms += excited_atoms * norm;
					expected_ground1_atoms += n3 * norm;
				}
			}
		}

		expected_ground2_atoms = N - expected_ground1_atoms - expected_excited_atoms;

		jump_table[0] = TimeStep * GammaUp   * expected_ground1_atoms;
		jump_table[1] = TimeStep * GammaUp   * expected_ground2_atoms;
		jump_table[2] = TimeStep * GammaDown * expected_excited_atoms;
		jump_table[3] = TimeStep * GammaDown * expected_excited_atoms;
		jump_table[4] = TimeStep * Kappa     * expected_photon_count;

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
	fprintf(stderr, "Allocating tables: %zu KB (%zu : %zu).\n", (TotalTableBytes * ThreadCount) >> 10, (sizeof LindbladJumps) >> 10, (9 * sizeof(wave_vector_t)) >> 10);

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

	fprintf(stderr, "Time per trajectory: %zu ms (%zu+%zu)\n", millis / TrajectoryCount, table_millis[0] / TrajectoryCount, table_millis[1] / TrajectoryCount);
	fprintf(stderr, "Jumps: %zu/%d (%.2f%%)\n", jump_count, TrajectoryCount * IntegrationSteps, 100.0f * jump_count / (TrajectoryCount * IntegrationSteps));
}
