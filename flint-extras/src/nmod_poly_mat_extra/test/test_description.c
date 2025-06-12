#include <stdlib.h>
#include <math.h> 
#include <time.h>
#include <flint/ulong_extras.h>


#include <pml/nmod_poly_mat_utils.h>
#include <pml/nmod_poly_mat_io.h>
//#include "../../nmod_poly_mat_extras/test/testing_collection.h"

#include "nmod_poly_mat_dixon.h"
#include "nmod_poly_mat_description.h"


// test one given input
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

    printf("\n");
    flint_printf("DIM: %ld",res);
    printf("\n");


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


int main(int argc, char ** argv)
{
    printf("Usage: %s OR %s [nbits] [rdim] [order] [Bcdim] [delta] \n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    slong nbits = atoi(argv[1]);
    slong rdim = atoi(argv[2]);
    slong order = atoi(argv[3]);
    slong Bcdim = atoi(argv[4]);
    slong delta = atoi(argv[5]);

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    slong prime = n_randprime(state, nbits, 1);
    printf("Launching  with\n\tprime = %ld,\n\trdim = %ld,\n\torder = %ld,\
            \n\tBcdim = %ld, \n\tdelta = %ld...\n",prime,rdim,order,Bcdim,delta);


    core_test_description(prime, rdim, order, Bcdim, delta, state);


    flint_rand_clear(state);

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
