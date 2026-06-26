/*
 * signal_random.h - Random and Stochastic Signals (L2/L8)
 *
 * White noise generation, AWGN, colored noise
 * Uniform and Gaussian random number generation
 * Monte Carlo signal analysis
 * Random process basics: ensemble, stationarity, ergodicity checks
 *
 * Course: MIT 6.003 Ch.10 / Stanford EE102A Ch.12 / Berkeley EE123 Ch.10
 * Textbook: Proakis & Manolakis (2007) Ch.12; Papoulis & Pillai (2002)
 */
#ifndef SIGNAL_RANDOM_H
#define SIGNAL_RANDOM_H

#include "signal_basis.h"

/* L2: PRNG state (Linear Congruential Generator) */
typedef struct {
    uint64_t state;
} signal_rng_t;

/* L2: RNG initialization and basic generators */
signal_rng_t *signal_rng_create(uint64_t seed);
void signal_rng_free(signal_rng_t *rng);
double signal_rng_uniform(signal_rng_t *rng);
double signal_rng_gaussian(signal_rng_t *rng, double mean, double std);

/* L2: Noise signal generation */
int signal_fill_white_noise(signal_ct_t *s, signal_rng_t *rng, double std);
int signal_fill_awgn(signal_ct_t *s, signal_rng_t *rng, double noise_std);
int signal_fill_colored_noise(signal_ct_t *s, signal_rng_t *rng, double std, double alpha);

/* L2: Add noise to existing signal (AWGN channel model) */
int signal_add_awgn(signal_ct_t *s, signal_rng_t *rng, double snr_db);

/* L2: Signal-to-Noise Ratio computation
 * SNR_dB = 10 * log10(P_signal / P_noise) */
double signal_compute_snr(const signal_ct_t *signal, const signal_ct_t *noisy);

/* L2: Random binary sequence generation (for digital comm simulation) */
int signal_fill_random_bits(signal_ct_t *s, signal_rng_t *rng, double bit_duration, double amplitude);

/* L8: Monte Carlo estimation of signal detection probability
 * For a matched filter detector in AWGN:
 *   P_d = Q(Q^{-1}(P_fa) - sqrt(2*SNR))
 * This uses Monte Carlo trials to estimate P_d for given P_fa and SNR. */
double signal_monte_carlo_pd(const signal_ct_t *template_signal, double snr_db,
                              double threshold, int num_trials, signal_rng_t *rng);

/* L8: Ergodicity check - compare time average vs ensemble average
 * A WSS process is mean-ergodic if time average converges to ensemble mean. */
int signal_check_ergodicity_mean(const signal_ct_t **ensemble, int num_realizations,
                                   double *time_avg, double *ensemble_avg, double *error);

#endif /* SIGNAL_RANDOM_H */

