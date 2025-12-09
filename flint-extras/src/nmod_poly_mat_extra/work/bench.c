#include <stdlib.h>
#include <math.h> 
//#include <flint/profiler.h>
#include <time.h>
#include <flint/ulong_extras.h>

#include <pml/nmod_poly_mat_utils.h>
#include <pml/nmod_poly_mat_io.h>
//#include "../../nmod_poly_mat_extras/test/testing_collection.h"

#include "nmod_poly_mat_kernel.h"


int main(int argc, char ** argv)
{
    printf("Usage: %s OR %s [nbits] [rdim] [cdim] [degree] [kappa] \n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    slong nbits = atoi(argv[1]);
    slong rdim = atoi(argv[2]);
    slong cdim = atoi(argv[3]);
    slong order = atoi(argv[4]);
    slong kappa = atoi(argv[5]);

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());


    slong prime = n_randprime(state, nbits, 1);

    prime=2;


    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, cdim, prime);

    nmod_poly_mat_randtest_sparse(A, state, order+1, 0.48);


    nmod_poly_mat_t N; // Not initialized 

    printf("Launching  with\n\tprime = %ld,\n\trdim = %ld,\n\tcdim = %ld,\
        \n\tdegree = %ld, \n\tkappa = %ld ...\n",prime,rdim,cdim,order,kappa);


    slong i;
    slong iz[rdim];

    
    slong tshift[cdim];

    for (i = 0; i < rdim; i++) 
        iz[i]=0; 

    // Initial modification of A

    //sortM(A,ishift,perm,iz);   


    // printf("Permutation \n [ ");
    // for (int j=0; j<cdim; j++) 
    //     printf(" %ld, ",perm[j]);
    // printf(" %ld ]\n",perm[cdim-1]);


    printf("A input \n");
    nmod_poly_mat_print_pretty(A, "x");
    printf("\n");


    double t = 0.0;
    clock_t tt;


    printf("~~~~WARMUP~~~~\n");
    nmod_poly_mat_zls(N, tshift, A, iz, kappa);
    printf("~~~~WARMUP DONE~~~~\n");

    tt = clock();
    nmod_poly_mat_zls(N, tshift, A, iz, kappa);
    t += (double)(clock()-tt) / CLOCKS_PER_SEC;



    // printf("N output \n");
    // nmod_poly_mat_print_pretty(N, "x");
    // printf("\n");


    slong nbnull;

    nbnull=N->c;

    nmod_poly_mat_t Z;
    nmod_poly_mat_init(Z, rdim, nbnull, A->modulus);

    nmod_poly_mat_mul(Z, A, N);

    // printf("\n");
    // nmod_poly_mat_print_pretty(Z, "x");
    // printf("\n");

    nmod_poly_mat_t Nflint;
    nmod_poly_mat_init(Nflint, cdim, cdim, A->modulus);
    
    nmod_poly_mat_nullspace(Nflint, A);
    
    double t2 = 0.0;
    tt = clock();
    nmod_poly_mat_nullspace(Nflint, A);
    t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;

    // printf("Nflint output \n");
    // nmod_poly_mat_print_pretty(Nflint, "x");
    // printf("\n");


    if (nmod_poly_mat_is_zero(Z) !=0) 
    {
        printf("Nullspace of dimension %ld, time: %f , time: %f \n",nbnull,t,t2);
        return 1;
    }



    nmod_poly_mat_clear(A);

    flint_rand_clear(state);

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
