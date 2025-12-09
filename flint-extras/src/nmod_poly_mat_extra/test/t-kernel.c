/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_extra.h"  

// test one given input
int core_test_kernel(const nmod_poly_mat_t mat)
{

    slong m=mat->r;
    slong n=mat->c;
    
    // Flint rank 
    // ----------
    
    slong rkflint;

    rkflint=nmod_poly_mat_rank(mat);
    

    // PML nullspace
    // -------------

    nmod_poly_mat_t N; // Not initialized, to modify 
    slong nz,nbnull;

    int i;
    slong iz[m];

    for (i = 0; i < m; i++) 
    {
        iz[i]=0; 
    }

    slong tshift[n];

    nz=nmod_poly_mat_zls(N, tshift, mat, iz, 2);

    if (nz==0) 
        nbnull=nz;
    else
        nbnull = N->c;

    int verif;
    verif =  (nz==nbnull) && (n-rkflint == nz);

    if (nz !=0) {

        nmod_poly_mat_t Z;
        nmod_poly_mat_init(Z, m, nbnull, mat->modulus);

        nmod_poly_mat_mul(Z, mat, N);

        verif = verif && nmod_poly_mat_is_zero(Z); 

        nmod_poly_mat_clear(Z);
        nmod_poly_mat_clear(N);
    }

    // printf("m %ld   n %ld\n",m,n);
    // printf("nz %ld\n",nz);
    // printf("nbnull %ld\n",nbnull);
    // printf("flint null %ld\n",n-rkflint);

    return verif; 

}

TEST_FUNCTION_START(nmod_poly_mat_kernel, state)
{
    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    int i,result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {

        ulong nbits = 2 + n_randint(state, 30);
        ulong rdim = 1 + n_randint(state, 20);
        ulong cdim = rdim + 1 + n_randint(state, 40);
        ulong deg = n_randint(state, 20);

        ulong prime = n_randprime(state, nbits, 1);


        nmod_poly_mat_t A;

        if (i < 4) {
            nmod_poly_mat_init(A, rdim, rdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 1.0);
        }
        else if (i < 8) {
            nmod_poly_mat_init(A, rdim, rdim+1, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.8);
        }
        else if (i < 50) {
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.2);
        }
        else {
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.84);
        }     

        result = core_test_kernel(A);

        nmod_poly_mat_clear(A);

        if (!result) {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

TEST_FUNCTION_END(state);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
