#include <flint/nmod_mat.h>
#include <flint/machine_vectors.h>
#include <flint/nmod_vec.h>

#include "nmod_extra.h"
#include "nmod_vec_extra.h"

/** matrix multiplication using AVX2 instructions for moduli less than 2^30 */
void nmod_mat_mul_small_modulus(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    // C aliases A or B: use a temporary
    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, A->r, B->c, A->mod.n);
        nmod_mat_mul_small_modulus(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return;
    }

    // precomputed values for computations modulo mod
    const double p = A->mod.n;
    const double pinv = 1 / (double) p;

    const vec2d p2 = vec2d_set_d2(p, p);
    const vec2d pinv2 = vec2d_set_d2(pinv, pinv);

    ulong power_two;
    NMOD_RED(power_two, UWORD(1) << 45, A->mod);
    ulong power_two_bis;
    NMOD_RED(power_two_bis, 1L<<DOT_SPLIT_BITS, A->mod);

    // transpose of B
    nmod_mat_t BT;
    nmod_mat_init(BT, B->c, B->r, B->mod.n);
    nmod_mat_transpose(BT, B);

    // let's go
    ulong res[2];
    slong i = 0;
    for (; i+1 < A->r; i += 2)
        for (slong j = 0; j < BT->r; j++)
        {
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
            _nmod_vec_dot2_small_modulus(res, A->rows[i], A->rows[i+1], BT->rows[j], A->c, power_two, p2, pinv2);
            C->rows[i][j] = res[0];
            C->rows[i+1][j] = res[1];
#else
            _nmod_vec_dot2_small_modulus(res,
                                         nmod_mat_entry_ptr(A, i, 0),
                                         nmod_mat_entry_ptr(A, i+1, 0),
                                         nmod_mat_entry_ptr(BT, j, 0),
                                         A->c, power_two, p2, pinv2);
            nmod_mat_entry(C, i, j) = res[0];
            nmod_mat_entry(C, i+1, j) = res[1];
#endif
        }

    for (; i < A->r; i++)
        for (slong j = 0; j < BT->r; j++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
            C->rows[i][j] = _nmod_vec_dot2_split(A->rows[i], BT->rows[j], A->c, A->mod, power_two_bis);
#else
            nmod_mat_entry(C, i, j) = _nmod_vec_dot2_split(nmod_mat_entry_ptr(A, i, 0),
                                                           nmod_mat_entry_ptr(BT, j, 0),
                                                           A->c, A->mod, power_two_bis);
#endif

    nmod_mat_clear(BT);
}

void nmod_mat_mul_nmod_vec_small_modulus(nn_ptr v, const nmod_mat_t A, nn_srcptr u, ulong len)
{
    // precomputed values for computations modulo mod
    const double p = A->mod.n;
    const double pinv = 1 / (double) p;

    const vec2d p2 = vec2d_set_d2(p, p);
    const vec2d pinv2 = vec2d_set_d2(pinv, pinv);

    ulong power_two;
    NMOD_RED(power_two, UWORD(1) << 45, A->mod);
    ulong power_two_bis;
    NMOD_RED(power_two_bis, 1L<<DOT_SPLIT_BITS, A->mod);

    // let's go
    ulong res[2];
    slong i = 0;
    for (; i+1 < A->r; i += 2)
    {
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        _nmod_vec_dot2_small_modulus(res, A->rows[i], A->rows[i+1], (nn_ptr)u, len, power_two, p2, pinv2);
#else
        _nmod_vec_dot2_small_modulus(res,
                                     nmod_mat_entry_ptr(A, i, 0),
                                     nmod_mat_entry_ptr(A, i+1, 0),
                                     (nn_ptr)u,
                                     len, power_two, p2, pinv2);
#endif
        v[i] = res[0];
        v[i+1] = res[1];
    }

    for (; i < A->r; i++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        v[i] = _nmod_vec_dot2_split(A->rows[i], u, len, A->mod, power_two_bis);
#else
        v[i] = _nmod_vec_dot2_split(nmod_mat_entry_ptr(A, i, 0), u,
                                    len, A->mod, power_two_bis);
#endif
}
