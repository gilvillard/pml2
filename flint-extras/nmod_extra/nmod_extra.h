#ifndef __NMOD_EXTRA__H
#define __NMOD_EXTRA__H

/*------------------------------------------------------------*/
/*        TODO: safe to assume this exists?                   */
/*------------------------------------------------------------*/

#include <flint/mpn_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/fft_small.h>

#if (defined __SIZEOF_INT128__ && GMP_LIMB_BITS == 64)
#define HAS_INT128
#define mp_dlimb_t unsigned __int128
#endif

/*------------------------------------------------------------*/
/* TODO: find the proper way to test                          */
/*------------------------------------------------------------*/
#ifdef __AVX2__
#define HAS_AVX2
#endif

/*------------------------------------------------------------*/
/* TODO: find which flags to test for                         */
/*------------------------------------------------------------*/

#ifdef HAS_AVX2 // TODO: use flint flags for avx512/avx2?
#include <immintrin.h>
#endif

#if (GMP_LIMB_BITS == 64)
#define mp_hlimb_t uint32_t     // half-limb
#define mp_qlimb_t uint16_t     // quarter-limb
#define mp_hlimb_signed_t int32_t      // signed half-limb
#define mp_qlimb_signed_t int16_t       // signed quarter-limb
#endif


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*        A few extra functionalities for Fp                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* returns the smallest i such that 2^i >= x                  */
/*------------------------------------------------------------*/
int next_power_of_two(ulong x);

/*------------------------------------------------------------*/
/* returns 1/p mod 2^k, assuming p is odd                     */
/* ill-defined when p is even                                 */
/*------------------------------------------------------------*/
ulong inverse_mod_power_of_two(ulong p, int k);

/*------------------------------------------------------------*/
/* finds an element of order at least n                       */
/* returns 0 if not found                                     */
/*------------------------------------------------------------*/
ulong nmod_find_root(long n, nmod_t mod);


#ifdef HAS_INT128



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* 64 bit modular multiplication using a preconditionner      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* high word of a*b, assuming words are 64 bits               */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong mul_hi(ulong a, ulong b)
{
    mp_dlimb_t prod;
    prod = a * (mp_dlimb_t)b;
    return prod >> 64;
}

/*------------------------------------------------------------*/
/* returns floor (2^64*b)/p                                   */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong prep_mul_mod_precon(ulong b, ulong p)
{
    return ((mp_dlimb_t) b << 64) / p;
}

/*------------------------------------------------------------*/
/* preconditioned product                                     */
/* returns ab mod p                                           */
/* a is in [0..2^64), b is in [0..p)                          */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong mul_mod_precon(ulong a, ulong b, ulong p, ulong i_b)
{
    ulong q;
    long t;

    q = mul_hi(i_b, a);
    t = a*b - q*p - p;
    if (t < 0)
        return t + p;
    else
        return t;
}

/*------------------------------------------------------------*/
/* returns ab mod p in [0..2p)                                */
/* a is in [0..2^64), b is in [0..p)                          */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong mul_mod_precon_unreduced(ulong a, ulong b, ulong p, ulong i_b)
{
    ulong q;

    q = mul_hi(i_b, a);
    return a*b - q*p;
}

#endif



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* A few extra operations for the fft_small types             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* returns a + b mod n, assuming a,b reduced mod n            */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE vec1n vec1n_addmod(vec1n a, vec1n b, vec1n n)
{
    return n - b > a ? a + b : a + b - n;
}

/*------------------------------------------------------------*/
/* returns a + b mod n, assuming a,b reduced mod n            */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE vec1d vec1d_addmod(vec1d a, vec1d b, vec1d n)
{
    return a + b - n >= 0 ? a + b - n : a + b;
}

/*------------------------------------------------------------*/
/* returns a + b mod n, assuming a,b reduced mod n            */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE vec4d vec4d_addmod(vec4d a, vec4d b, vec4d n)
{
    return vec4d_reduce_2n_to_n(vec4d_add(a, b), n);
}

/*------------------------------------------------------------*/
/* TODO: if AVX512 supported, use cvtepi64_pd instead         */
/* loads a vec4n from a and converts it to double             */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE vec4d vec4d_load_unaligned_nn_ptr(nn_ptr a)
{
#ifdef FLINT_HAVE_AVX512
    return  _mm256_setr_m128d( _mm_cvtepi64_pd(_mm_loadu_si128((vec2n *) a)),
                               _mm_cvtepi64_pd(_mm_loadu_si128((vec2n *) (a + 2))) );
#else // TODO need to test if has avx2?
// when AVX2 is available, this is a bit slower than the solution above
    return vec4n_convert_limited_vec4d(vec4n_load_unaligned(a));
#endif
}

/*------------------------------------------------------------*/
/* TODO: if AVX512 supported, use cvtpd_epi64 instead         */
/* converts a vec4d to vec4n and stores it                    */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE void vec4d_store_unaligned_nn_ptr(nn_ptr dest, vec4d a)
{
    vec4n_store_unaligned(dest, vec4d_convert_limited_vec4n(a));
}


FLINT_FORCE_INLINE void vec4n_store_aligned(ulong* z, vec4n a)
{
    _mm256_store_si256((__m256i*) z, a);
}

/* reduce_to_pm1n(a, n, ninv): return a mod n in (-n,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_to_pm1no(vec2d a, vec2d n, vec2d ninv)
{                                                               
    return _mm_fnmadd_pd(_mm_round_pd(_mm_mul_pd(a, ninv), 4), n, a); 
}

/* reduce_pm1no_to_0n(a, n): return a mod n in [0,n) assuming a in (-n,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_pm1no_to_0n(vec2d a, vec2d n)
{
    return _mm_blendv_pd(a, _mm_add_pd(a, n), a); 
}

/* reduce_to_0n(a, n, ninv): return a mod n in [0,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_to_0n(vec2d a, vec2d n, vec2d ninv)
{ 
    return vec2d_reduce_pm1no_to_0n(vec2d_reduce_to_pm1no(a, n, ninv), n);
}

FLINT_FORCE_INLINE vec2d vec2d_set_d2(double a1, double a0)
{
    return _mm_set_pd(a0, a1);
}

#define vec4n_bit_shift_right_45(a) vec4n_bit_shift_right((a), 45)

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                    CRT and multimod                        */
/* -we only handle 1 to 4 hard-coded primes p0..pk            */
/* -we take as extra parameter another modulus p              */           
/* -multimod reduces vectors of ulong mod all pi's        */
/* -CRT does Chinese Remainders, and reduces the result mod p */
/* -use AVX2 for p < 2^50, long arithmetic otherwise          */
/* TODO: set a macro for testing large/small modulus          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// 4 primes should be enough for all practical purposes involving nmods
// TODO: make sure we test if so, just in case?
//
// PRIME0 < PRIME1 < PRIME2 < PRIME3 needed for the double implementation
#define PRIME0 659706976665601
#define PRIME1 910395627798529
#define PRIME2 1086317488242689
#define PRIME3 1108307720798209

typedef struct {
    ulong num_primes;
    ulong p;
    nmod_t mod;

    nn_ptr data;
    double pinv;
    nmod_t mod_primes[4];
    ulong primes[4];
    double primes_inv[4];
    ulong inverse_cofactors[4];
    double p0_red, p1_red, p0p1_red, p0p1p2_red, p0_red2, p0p1red_2, invp0_p1, invp0p1_p2, p0p1_red3, invp0p1p2_p3;
} nmod_multimod_CRT_struct;

typedef nmod_multimod_CRT_struct nmod_multimod_CRT_t[1];


/*------------------------------------------------------------*/
/* initializes all data in C                                  */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_init(nmod_multimod_CRT_t C, ulong p, ulong num_primes);

/*------------------------------------------------------------*/
/* clears all data in C                                       */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_clear(nmod_multimod_CRT_t C);

/*------------------------------------------------------------*/
/* out[i] = CRT(residues[j][i], j < num_primes) mod p, i < nb */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_CRT(nn_ptr out, nn_ptr *residues, ulong nb, nmod_multimod_CRT_t C);

/*------------------------------------------------------------*/
/* residues[j][i] = input[i] mod prime[j]                     */
/* for i < nb, j < num_primes                                 */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_reduce(nn_ptr *residues, nn_ptr input, ulong nb, nmod_multimod_CRT_t C);


#endif
