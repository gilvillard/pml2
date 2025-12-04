#include <stdlib.h>
#include <math.h> 
//#include <flint/profiler.h>
#include <time.h>
#include <flint/ulong_extras.h>

#include <pml/nmod_poly_mat_utils.h>
#include <pml/nmod_poly_mat_io.h>
//#include "../../nmod_poly_mat_extras/test/testing_collection.h"



int main(int argc, char ** argv)
{
    printf("Usage: %s OR %s [nbits] [rdim] [cdim] [degree] \n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    slong nbits = atoi(argv[1]);
    slong rdim = atoi(argv[2]);
    slong cdim = atoi(argv[3]);
    slong order = atoi(argv[4]);

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());


    slong prime = n_randprime(state, nbits, 1);

    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, cdim, prime);

    nmod_poly_mat_randtest_sparse(mat, state, order+1, 0.05);


    nmod_poly_mat_t N; // Not initialized 

    printf("Launching  with\n\tprime = %ld,\n\trdim = %ld,\n\tcdim = %ld,\
            \n\tdegree = %ld ...\n",prime,rdim,cdim,order);



static inline void _test_collection_mat_sparse(nmod_poly_mat_t mat, slong deg, flint_rand_t state)
{
   
}


    

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
