#include <stdlib.h>
//#include <flint/profiler.h>
#include <time.h>
#include <flint/ulong_extras.h>
#include <flint/test_helpers.h>

#include <nmod_poly_mat_utils.h>
#include <nmod_poly_mat_io.h>

//#include "../../nmod_poly_mat_extras/test/testing_collection.h"

#include "nmod_poly_mat_description.h"




// Random poly_mat which is nonsingular for x=0
void nmod_poly_mat_rand_origin(nmod_poly_mat_t A, flint_rand_t state, slong order) 
{

    nmod_mat_t T1,T2;
    nmod_mat_init(T1, A->r, A->r, A->modulus);
    nmod_mat_init(T2, A->r, A->r, A->modulus);

    nmod_mat_randtril(T1,state,0);
    nmod_mat_randtriu(T2,state,0);

    nmod_mat_mul(T1,T1,T2);

    nmod_poly_mat_t B;
    nmod_poly_mat_init(B, A->r, A->r, A->modulus);

    // nmod_poly_mat_set_nmod_mat(B,T1); // GV does not work, bug ? Sam 31 mai 2025 14:29:48 CEST

    for (slong i = 0; i < A->r; i++)
    {
        for (slong j = 0; j < A->r; j++)
        {
            if (nmod_mat_entry(T1, i, j) == 0)
                nmod_poly_zero(nmod_poly_mat_entry(B, i, j));
            else
            {
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, i, j),0,nmod_mat_entry(T1, i, j));
            }
        }
    }

    nmod_poly_mat_rand(A, state, order-1);
    nmod_poly_mat_shift_left(A,A,1);

    nmod_poly_mat_add(A,A,B);

    nmod_mat_clear(T1);
    nmod_mat_clear(T2); 
    nmod_poly_mat_clear(B); 
}



// test one given input
int core_test_description(slong prime, slong rdim, slong order, slong Bcdim, slong delta, flint_rand_t state)
{
    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, rdim, prime);
    nmod_poly_mat_rand_origin(A, state, order+1);

    nmod_poly_mat_t B;
    nmod_poly_mat_init(B, rdim, Bcdim, A->modulus);
    nmod_poly_mat_rand(B, state, order+1);

    nmod_poly_mat_t X;
    nmod_poly_mat_init(X, rdim, B->c, A->modulus);

    //nmod_poly_mat_inv_trunc(S,A,order-1);

    slong sigma;
    sigma = ceil((double) (rdim + Bcdim)*delta/rdim +1);   // sigma >= ceil((double) (rdim + Bcdim)*delta/rdim +1);


    nmod_poly_mat_dixon(X, A, B, order, sigma);    

    // printf("\n");
    // nmod_poly_mat_print_pretty(A, "x");
    // printf("\n");
    // nmod_poly_mat_print_pretty(B, "x");
    // printf("sigma %ld\n",sigma); 
    // printf("\n");


    nmod_poly_mat_t N;
    nmod_poly_mat_init(N, rdim, Bcdim, A->modulus);

    nmod_poly_mat_t D;
    nmod_poly_mat_init(D, Bcdim, Bcdim, A->modulus);

    int res=0;
    res=nmod_poly_mat_description(N, D, X, delta);

    if (res != 0)
    {
        nmod_poly_mat_t T1;
        nmod_poly_mat_init(T1, rdim, Bcdim, A->modulus);

        nmod_poly_mat_t T2;
        nmod_poly_mat_init(T2, rdim, Bcdim, A->modulus);


        nmod_poly_mat_mul(T1,A,N);
        nmod_poly_mat_mul(T2,B,D);

        nmod_poly_mat_sub(T1,T1,T2);

        printf("\n");
        nmod_poly_mat_print_pretty(T1, "x");
        printf("\n");

        if (nmod_poly_mat_is_zero(T1) !=0) 
        {
            nmod_poly_mat_clear(T1); 
            return 1;
        }
        
    }
    
    return 0;
   
}


TEST_FUNCTION_START(nmod_poly_mat_descriptions, state)
{

    int res=0;  

    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());


    core_test_description(5, 8, 1, 1, 8, state);


    res=collection_test_dixon(state);

    if (res == 0)
    {
       TEST_FUNCTION_FAIL("");
    }
    else
    {
    TEST_FUNCTION_END(state);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
