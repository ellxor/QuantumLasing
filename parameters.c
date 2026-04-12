#pragma once
#define _GNU_SOURCE
#include <math.h>

#ifndef OPTION_N
#  define OPTION_N 32
#endif

#ifndef OPTION_O
#  define OPTION_O 0.225f
#endif

constexpr int N = OPTION_N;

// Hamiltonian coefficients
constexpr float Omega = 1;				// energy of microwave driving frequency
constexpr float OmegaC = 0;				// cavity energy (of single photon)
constexpr float OmegaE = Omega;				// excited state energy
constexpr float GCoupling = 0.9f * Omega / sqrtf(N);	// cavity-atom coupling strength
constexpr float Phi = 0;				// phase shift of microwave drive

// Lindblad decay rates
constexpr float Kappa     = 0.8f * Omega;
constexpr float GammaUp   = OPTION_O * Omega;
constexpr float GammaDown = Omega/2 - GammaUp;

// Integration parameters
constexpr float TimeStep = 1e-4f;
constexpr int IntegrationSteps = 100'000;

// Misc.
constexpr int TrajectoryCount = 1;
constexpr int ThreadCount = 1;
