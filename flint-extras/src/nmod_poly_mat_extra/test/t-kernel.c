/*
    Copyright (C) 2025 Gilles Villard 

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/
#include <math.h> 

#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_extra.h"  // det_iter currently in _extra.h

#include "nmod_poly_mat_description.h"


// test one given input
int core_test_kernel(const nmod_poly_mat_t A, slong delta)
{

    nmod_poly_mat_t N;

    slong nbnull;

    nbnull = nmod_poly_mat_kernel(N, A, delta);

    if (nbnull == 0)
        return 0;

    nmod_poly_mat_t T;
    nmod_poly_mat_init(T, A->r, nbnull, A->modulus);

    nmod_poly_mat_mul(T,A,N);

    // printf("--------- %ld \n",nbnull);
    // nmod_poly_mat_print_pretty(N, "x");
    // printf("\n");

    if (nmod_poly_mat_is_zero(T) == 0)
        {
            
            nbnull = 0; 
        }

    nmod_poly_mat_clear(N); 
    nmod_poly_mat_clear(T);

    return nbnull;
}



TEST_FUNCTION_START(nmod_poly_mat_kernel, state)
{
    ulong nbits; 
    ulong rdim;
    ulong cdim;
    ulong deg; 
    ulong delta; 
    ulong prime;

    int i, res;

    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    if (res == 0)
            TEST_FUNCTION_FAIL("");

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nbits = 2 + n_randint(state, 62);
        rdim = 1 + n_randint(state, 20);
        deg = 1+ n_randint(state, 40);

        prime = n_randprime(state, nbits, 1);

        cdim = rdim + 1 + n_randint(state, 20);

        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, rdim, cdim, prime);

        flint_printf("-- rdim %ld  cdim %ld   deg %ld\n", rdim, cdim, deg); 


        //nmod_poly_mat_rand(A, state, deg+1);

        if (i < 30)
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.9);
        else if (i < 60)
            nmod_poly_mat_rand(A, state, deg+1);
        else
            nmod_poly_mat_randtest(A, state, deg+1);

        //++++++++++
            nmod_poly_mat_t T;
            nmod_poly_mat_init(T, A->r, A->r, A->modulus);

            for (int ii = 0; ii < rdim; ii++)
            {
                for (int j = 0; j < rdim; j++)
                nmod_poly_set(nmod_poly_mat_entry(T, ii, j), nmod_poly_mat_entry(A, ii, j));
            }


            nmod_poly_mat_truncate(T,2);

            nmod_poly_t det;
            nmod_poly_init(det, A->modulus);
            nmod_poly_mat_det(det, T);

            printf("DET \n");
            nmod_poly_print_pretty(det,"x");
            printf("\n");

            nmod_poly_mat_clear(T);
        //+++++++++++++=

        delta = floor((double) (rdim*deg)/(cdim-rdim));
        flint_printf("delta %ld\n", delta); 
        flint_printf("prime %ld\n", prime); 
        res = core_test_kernel(A, delta);
        flint_printf("res %ld  diff %ld\n", res, cdim-rdim); 
        nmod_poly_mat_clear(A);

        if (res ==0)
            TEST_FUNCTION_FAIL("");
    }

    TEST_FUNCTION_END(state);
}
