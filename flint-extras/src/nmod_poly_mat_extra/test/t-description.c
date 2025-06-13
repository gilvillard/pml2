#include <stdlib.h>
#include <math.h> 

#include <time.h>
#include <flint/ulong_extras.h>
#include <flint/test_helpers.h>

#include <nmod_poly_mat_utils.h>
#include <nmod_poly_mat_io.h>

#include "nmod_poly_mat_extra.h"



// test one given input: here, description of a polynomial matrix ! 
int poly_test_description(slong prime, slong rdim, slong cdim, slong order, slong delta, flint_rand_t state)
{
    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, cdim, prime);
    nmod_poly_mat_rand(A, state, order+1);

    slong sigma;
    sigma = ceil((double) (rdim + cdim)*(delta+1)/rdim +1);  


    nmod_poly_mat_t N;
    nmod_poly_mat_init(N, rdim, cdim, A->modulus);

    nmod_poly_mat_t D;
    nmod_poly_mat_init(D, cdim, cdim, A->modulus);

    int res=0;
    res=nmod_poly_mat_description(N, D, A, delta);

    if (res != 0)
    {
        nmod_poly_mat_t T1;
        nmod_poly_mat_init(T1, rdim, cdim, A->modulus);

        nmod_poly_mat_mul(T1,A,D);
        nmod_poly_mat_sub(T1,T1,N);

        if (nmod_poly_mat_is_zero(T1) !=0) 
            res=1;
        else 
            res=0;

        nmod_poly_mat_clear(T1); 
        }
    
    nmod_poly_mat_clear(A); 
    nmod_poly_mat_clear(N); 
    nmod_poly_mat_clear(D); 
    
    return res;
}

// test one given input: here, description of a known fraction 
int core_test_description(slong prime, slong rdim, slong order, slong Bcdim, slong delta, flint_rand_t state)
{
    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, rdim, prime);
    nmod_poly_mat_rand_at_zero(A, state, order+1);

    nmod_poly_mat_t B;
    nmod_poly_mat_init(B, rdim, Bcdim, A->modulus);
    nmod_poly_mat_rand(B, state, order+1);

    nmod_poly_mat_t X;
    nmod_poly_mat_init(X, rdim, B->c, A->modulus);

    slong sigma;
    sigma = ceil((double) (rdim + Bcdim)*(delta+1)/rdim +1);   


    nmod_poly_mat_dixon(X, A, B, order, sigma);    

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

        nmod_poly_mat_clear(T2); 

        if (nmod_poly_mat_is_zero(T1) !=0) 
            res=1;
        else 
            res=0;

        nmod_poly_mat_clear(T1); 
        }
    
    nmod_poly_mat_clear(A); 
    nmod_poly_mat_clear(B); 
    nmod_poly_mat_clear(X); 
    nmod_poly_mat_clear(N); 
    nmod_poly_mat_clear(D); 

    return res;
}

// TODO, matrix allocation once 
TEST_FUNCTION_START(nmod_poly_mat_description, state)
{
    int i;

    int res=0;  

    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    slong order;
    slong Bcdim;
    slong delta;
    slong prime;
    slong rdim;
    slong cdim;
    slong nbits;

 /** ------------------------- */
    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {

        nbits = 2 + n_randint(state, 62);
        rdim = 1 + n_randint(state, 20);
        cdim = 1 + n_randint(state, 10);
        order = 1+ n_randint(state, 20);
        prime = n_randprime(state, nbits, 1);

        delta=floor((double) rdim*(order)/cdim);
        res=poly_test_description(prime, rdim, cdim, order, delta, state);
        
        if (res == 0)
            TEST_FUNCTION_FAIL("");
    }

    /** ------------------------- */
    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        prime=2; 
        rdim = 1 + n_randint(state, 80);
        order=1 + n_randint(state, 4);
        Bcdim=1 + n_randint(state, 4);
        delta=floor((double) rdim*(order)/Bcdim);
        res=core_test_description(prime, rdim, order, Bcdim, delta, state);

        if (res == 0)
            TEST_FUNCTION_FAIL("");
    }
    
    /** ------------------------- */
    nbits=10;
    prime=n_randprime(state, nbits, 1);
    rdim = 1 + n_randint(state, 40);
    order=1 + n_randint(state, 4);
    Bcdim=1 + n_randint(state, 4);
    delta=floor((double) rdim*(order)/Bcdim);
    res=core_test_description(prime, rdim, order, Bcdim, delta, state);
    
    if (res == 0)
       TEST_FUNCTION_FAIL("");

    /** ------------------------- */
   for (i = 0; i < 4 * flint_test_multiplier(); i++)
   {
    nbits=40;
    prime=n_randprime(state, nbits, 1);
    rdim = 1 + n_randint(state, 40);
    order=1 + n_randint(state, 20);
    Bcdim=1 + n_randint(state, 6);
    delta=floor((double) rdim*(order)/Bcdim);
    res=core_test_description(prime, rdim, order, Bcdim, delta, state);
    
    if (res == 0)
        TEST_FUNCTION_FAIL("");
}
    
    /** ------------------------- */
   for (i = 0; i < 4 * flint_test_multiplier(); i++)
   {
        nbits=80;
        prime=n_randprime(state, nbits, 1);
        rdim = 1 + n_randint(state, 10);
        order=1 + n_randint(state, 10);
        Bcdim=rdim + n_randint(state, 20);
        delta=floor((double) rdim*(order)/Bcdim);
        res=core_test_description(prime, rdim, order, Bcdim, delta, state);
    
        if (res == 0)
        TEST_FUNCTION_FAIL("");
    }

#if FLINT64

    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        nbits=60;
        prime=n_randprime(state, nbits, 1);
        rdim = 1 + n_randint(state, 20);
        order=1 + n_randint(state, 12);
        Bcdim=1 + n_randint(state, 8);
        delta=floor((double) rdim*(order)/Bcdim);
        res=core_test_description(prime, rdim, order, Bcdim, delta, state);
    
        if (res == 0)
        TEST_FUNCTION_FAIL("");
    }

   /** ------------------------- */
    nbits=60;
    prime=n_randprime(state, nbits, 1);
    rdim = 1 + n_randint(state, 20);
    order=1 + n_randint(state, 12);
    Bcdim=rdim + n_randint(state, 20);
    delta=floor((double) rdim*(order)/Bcdim);
    res=core_test_description(prime, rdim, order, Bcdim, delta, state);
    
    if (res == 0)
       TEST_FUNCTION_FAIL("");

#endif  

TEST_FUNCTION_END(state);
    
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
