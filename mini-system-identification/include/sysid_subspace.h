/**
 * @file sysid_subspace.h
 * @brief Subspace state-space system identification (L5 Algorithms, L8 Advanced)
 *
 * Knowledge coverage:
 *   L5 - N4SID algorithm (Van Overschee & De Moor, 1994)
 *   L5 - MOESP algorithm (Verhaegen & Dewilde, 1992)
 *   L5 - CVA (Canonical Variate Analysis) method (Larimore, 1990)
 *   L5 - Block Hankel matrix construction from I/O data
 *   L3 - Oblique projection in subspace methods
 *   L3 - SVD for order estimation and state sequence recovery
 *   L8 - Closed-loop subspace identification (SSARX, PBSID)
 *   L8 - Recursive subspace identification via propagator method
 *
 * Subspace methods solve the state-space identification problem in one step,
 * without iterative optimization:
 *   1. Construct block Hankel matrices from I/O data
 *   2. Perform oblique projection to eliminate noise
 *   3. SVD to estimate state sequence and system order
 *   4. Least-squares to extract (A,B,C,D) matrices
 *
 * The fundamental equation (Van Overschee & De Moor 1996):
 *   Y_f = Γ_i X_f + H_i U_f + N_f
 * where Γ_i is the extended observability matrix and H_i is the Toeplitz matrix.
 *
 * References:
 *   Van Overschee, P. & De Moor, B. (1996) "Subspace Identification for Linear Systems"
 *   Verhaegen, M. & Dewilde, P. (1992) "Subspace model identification"
 *   Larimore, W.E. (1990) "Canonical variate analysis in identification"
 */

#ifndef SYSID_SUBSPACE_H
#define SYSID_SUBSPACE_H

#include "sysid_types.h"
#include "sysid_models.h"

/*---------------------------------------------------------------------------
 * Subspace algorithm options
 *---------------------------------------------------------------------------*/

typedef enum {
    SYSID_SUB_N4SID  = 0,  /**< Numerical Algorithms for Subspace State Space ID */
    SYSID_SUB_MOESP  = 1,  /**< Multivariable Output-Error State sPace */
    SYSID_SUB_CVA    = 2   /**< Canonical Variate Analysis */
} sysid_subspace_method_t;

/** Subspace identification options */
typedef struct {
    sysid_subspace_method_t method;  /**< Algorithm choice */
    int  i;              /**< Block rows (past/future horizon), typically 10-20 */
    int  nx_max;         /**< Maximum state dimension to consider */
    double svd_tol;      /**< SVD singular value truncation tolerance */
    int  use_weights;    /**< Use CVA or MOESP weighting matrices */
} sysid_subspace_opts_t;

/*---------------------------------------------------------------------------
 * L3: Block Hankel matrix construction
 *---------------------------------------------------------------------------*/

/**
 * Construct block Hankel matrix from a scalar signal sequence
 *
 *         [ s(0)   s(1)   ...  s(j-1) ]
 *         [ s(1)   s(2)   ...  s(j)   ]
 *   H =   [  ...    ...   ...   ...    ]
 *         [ s(i-1) s(i)   ...  s(i+j-2)]
 *
 * @param s       Signal vector [N]
 * @param N       Signal length
 * @param i       Number of block rows
 * @param j       Number of columns (typically j = N - 2*i + 1)
 * @param H       Output Hankel matrix [i × j], row-major
 */
int  sysid_block_hankel(const double *s, int N, int i, int j, double *H);

/** Construct block Hankel from multi-channel signal (nu outputs, row-major input) */
int  sysid_block_hankel_multi(const double *s, int N, int nu, int i, int j, double *H);

/** Build data equation matrices for subspace ID:
 *   U_p, U_f: past/future input Hankel (i×nu rows each)
 *   Y_p, Y_f: past/future output Hankel (i×ny rows each)
 */
int  sysid_subspace_data_matrices(const double *u, const double *y, int N,
                                   int nu, int ny, int i,
                                   double *Up, double *Uf,
                                   double *Yp, double *Yf, int *j);

/*---------------------------------------------------------------------------
 * L5: N4SID algorithm
 *
 * Algorithm (Van Overschee & De Moor 1994):
 *   1. Compute oblique projection O_i = Y_f /_{U_f} W_p (with W_p = [U_p; Y_p])
 *   2. SVD: O_i = U Σ Vᵀ using CVA weighting: W1 O_i W2
 *   3. Estimate order n from singular value gap
 *   4. Extended observability Γ_i = U Σ^{1/2}
 *   5. State sequence X_i = Γ_i^† O_i
 *   6. Least squares: [X_{i+1}; Y_i] = [A B; C D] [X_i; U_i]
 *---------------------------------------------------------------------------*/

/**
 * N4SID algorithm for state-space identification
 *
 * @param u        Input data [N × nu], row-major
 * @param y        Output data [N × ny], row-major
 * @param N        Number of samples
 * @param nu       Number of inputs
 * @param ny       Number of outputs
 * @param opts     Algorithm options
 * @param model    Output: identified state-space model
 * @param nx_out   Output: estimated state dimension (from SVD gap)
 * @param sv       Output: singular values [opts.i*ny], plot to see gap
 * @return         0 on success, -1 on error
 */
int  sysid_n4sid(const double *u, const double *y, int N, int nu, int ny,
                 const sysid_subspace_opts_t *opts,
                 sysid_ss_t *model, int *nx_out, double *sv);

/**
 * MOESP algorithm (Verhaegen & Dewilde 1992)
 * Uses RQ factorization for numerical efficiency.
 */
int  sysid_moesp(const double *u, const double *y, int N, int nu, int ny,
                 const sysid_subspace_opts_t *opts,
                 sysid_ss_t *model, int *nx_out, double *sv);

/**
 * CVA algorithm (Larimore 1990)
 * Uses canonical correlation weighting for optimal accuracy.
 */
int  sysid_cva(const double *u, const double *y, int N, int nu, int ny,
               const sysid_subspace_opts_t *opts,
               sysid_ss_t *model, int *nx_out, double *sv);

/** Generic subspace call with method selection */
int  sysid_subspace_id(const double *u, const double *y, int N, int nu, int ny,
                       const sysid_subspace_opts_t *opts,
                       sysid_ss_t *model, int *nx_out, double *sv);

/*---------------------------------------------------------------------------
 * L8: Advanced — Closed-loop subspace identification (SSARX)
 *
 * Uses high-order ARX to pre-estimate the plant, then reduces to state-space.
 * Consistent under feedback (Jansson, 2005).
 *---------------------------------------------------------------------------*/

/** SSARX: Subspace ID via ARX for closed-loop data */
int  sysid_ssarx(const double *u, const double *y, int N, int nu, int ny,
                 const sysid_subspace_opts_t *opts,
                 sysid_ss_t *model, int *nx_out);

/*---------------------------------------------------------------------------
 * Utility: Estimate system order from singular values (gap detection)
 *---------------------------------------------------------------------------*/

/** Detect order as index of maximum singular value gap ratio */
int  sysid_detect_order(const double *sv, int n_sv, double *gap_ratio);

/** Initialize subspace options with sensible defaults */
void sysid_subspace_opts_init(sysid_subspace_opts_t *opts);

#endif /* SYSID_SUBSPACE_H */
