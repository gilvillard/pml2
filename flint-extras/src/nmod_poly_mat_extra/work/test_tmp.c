#include <stdlib.h>
#include <math.h> 
//#include <flint/profiler.h>
#include <time.h>
#include <flint/ulong_extras.h>

#include <pml/nmod_poly_mat_utils.h>
#include <pml/nmod_poly_mat_io.h>
//#include "../../nmod_poly_mat_extras/test/testing_collection.h"

#include "nmod_poly_mat_dixon.h"
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



int main(int argc, char ** argv)
{
    printf("Usage: %s OR %s [nbits] [rdim] [cdim] [order] [delta] \n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    slong nbits = atoi(argv[1]);
    slong rdim = atoi(argv[2]);
    slong cdim = atoi(argv[3]);
    slong order = atoi(argv[4]);
    slong delta = atoi(argv[5]);

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    slong prime = n_randprime(state, nbits, 1);
    printf("Launching  with\n\tprime = %ld,\n\trdim = %ld,\n\tcdim = %ld,\
            \n\torder = %ld, \n\tdelta = %ld...\n",prime,rdim,cdim,order,delta);

    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, cdim, prime);
    nmod_poly_mat_rand(A, state, order+1);

    nmod_poly_mat_t N; // Not initialized 

    slong nbnull;

    nbnull = nmod_poly_mat_kernel(N,A,delta);

    if (nbnull > 0) 
    {
       printf("--------- %ld \n",nbnull);
       //nmod_poly_mat_print_pretty(N, "x");
       printf("\n");
   }
   else 
        flint_printf("\nNo kernel vector found of degree %ld or less\n\n", delta);

    nmod_poly_mat_clear(A);

    flint_rand_clear(state);

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
