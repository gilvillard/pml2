#include <flint/nmod_mat.h>
#include "nmod_mat_poly.h"

void nmod_mat_poly_mul_coeff(nmod_mat_t coeff,
                             const nmod_mat_poly_t mat1,
                             const nmod_mat_poly_t mat2,
                             slong k)
{
    // only consider indices i such that:
    //     0 <= i <= k
    //     i < mat1->length
    //     k-i < mat2->length
    // so  i < min(k+1,mat1->length)
    // and i >= max(0, k + 1 - mat2->length)
    const slong ubound = FLINT_MIN(k+1, mat1->length);
    const slong lbound = FLINT_MAX(0, k+1 - mat2->length);

    // lbound >= ubound ==> coeff is zero, no term to consider
    if (lbound >= ubound)
    {
        nmod_mat_zero(coeff);
        return;
    }

    // now lbound < ubound
    // first handle i == lbound separately, to avoid wasting time zero-ing `coeff`
    nmod_mat_mul(coeff, mat1->coeffs + lbound, mat2->coeffs + (k - lbound));

    // `if` just here to avoid initializing temp for nothing
    if (lbound + 1 < ubound)
    {
        nmod_mat_t temp;
        nmod_mat_init(temp, mat1->r, mat2->c, mat1->mod.n);
        for (slong i = lbound+1; i < ubound; i++)
        {
            nmod_mat_mul(temp, mat1->coeffs + i, mat2->coeffs + (k - i));
            nmod_mat_add(coeff, coeff, temp);
        }
        nmod_mat_clear(temp);
    }
}

void nmod_mat_poly_evaluate_nmod(nmod_mat_t eval,
                                 const nmod_mat_poly_t matp,
                                 mp_limb_t pt)
{
    slong k = matp->length;

    if (k == 0)
    {
        nmod_mat_zero(eval);
        return;
    }

    if (k == 1 || pt == 0)
    {
        nmod_mat_set(eval, matp->coeffs + 0);
        return;
    }

    k--; // k == degree
    nmod_mat_set(eval, matp->coeffs + k);
    k--; // k == degree-1

    // Horner: eval = matp[k] + eval*pt, k = degree-1 ... 0
    for ( ; k >= 0; k--)
        nmod_mat_scalar_addmul_ui(eval, nmod_mat_poly_coeff(matp, k), eval, pt);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
