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
    printf("Usage: %s OR %s [nbits] [rdim] [cdim] [degree] [kappa]\n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    slong nbits = atoi(argv[1]);
    slong rdim = atoi(argv[2]);
    slong cdim = atoi(argv[3]);
    slong order = atoi(argv[4]);
    double kappa = atof(argv[5]);

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());


    slong prime = n_randprime(state, nbits, 1);


    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, cdim, prime);

    nmod_poly_mat_randtest_sparse(A, state, order+1, 0.9);
    //nmod_poly_mat_randtest(A, state, order+1);

    nmod_poly_mat_t N; 
    nmod_poly_mat_init(N, cdim, cdim, A->modulus);

    printf("Launching  with\n\tprime = %ld,\n\trdim = %ld,\n\tcdim = %ld,\
        \n\tdegree = %ld, \n\tkappa = %.1f ...\n",prime,rdim,cdim,order,kappa);


    slong i,j;
    slong iz[cdim];

    
    slong degN[cdim];

    double t = 0.0;
    clock_t tt;


    slong nz;
 

    //printf("~~~~WARMUP~~~~\n");
    //nz=nmod_poly_mat_kernel(N, degN, A, NULL, kappa);
    //printf("~~~~WARMUP DONE~~~~\n");


    for (i = 0; i < 2; i++) 
        iz[i]=0; 
    for (i = 2; i < cdim; i++) 
        iz[i]=0;

    tt = clock();
    nz=nmod_poly_mat_kernel(N, degN, A, iz, kappa);
    //nz=nmod_poly_mat_approximant_kernel(N, degN, A, iz);
    t += (double)(clock()-tt) / CLOCKS_PER_SEC;


    nmod_poly_mat_t Nflint;
    nmod_poly_mat_init(Nflint, cdim, cdim, A->modulus);
    
    //nmod_poly_mat_nullspace(Nflint,A);
    
    slong nullflint;

    double t2 = 0.0;
    tt = clock();
    nullflint=nmod_poly_mat_nullspace(Nflint,A);
    //nmod_poly_mat_approximant_kernel(N, degN, A, iz);
    t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;

    
    if (nz != nullflint){

        printf("\n Differ from flint\n");
            return 0;

    }

    if (nz !=0) {

        nmod_poly_mat_t NN; 
        nmod_poly_mat_init(NN, cdim, nz, A->modulus);


        for (i = 0; i < cdim; i++)
            for (j = 0; j < nz; j++) {
                nmod_poly_set(nmod_poly_mat_entry(NN, i, j), nmod_poly_mat_entry(N, i, j));
            }

        // ++++++++++++++++++++
            // printf("Kernel -- \n");
            // nmod_poly_mat_print_pretty(NN, "x");
            // printf("\n");

// printf("NN output 0\n");
//     nmod_poly_mat_print_pretty(NN, "x");
//     printf("\n");

            slong pivind[nz];

            nmod_poly_mat_pivot_index(pivind,NN,iz,COL_UPPER);


// printf("NN output \n");
//     nmod_poly_mat_print_pretty(NN, "x");
//     printf("\n");


            // printf("Pivind \n [ ");
            // for (int j=0; j<nz-1; j++) 
            //     printf(" %ld, ",pivind[j]);
            // printf(" %ld ]\n",pivind[nz-1]);

            // printf("degN \n [ ");
            // for (int j=0; j<nz-1; j++) 
            //     printf(" %ld, ",degN[j]);
            // printf(" %ld ]\n",degN[nz-1]);


            int is_weak_popov;

            is_weak_popov=nmod_poly_mat_is_ordered_weak_popov(NN, iz, COL_UPPER);
                            
            printf("\n is ordered weak Popov: %d \n ",is_weak_popov);

        // +++++++++++++++++++


            nmod_poly_mat_t Z;
            nmod_poly_mat_init(Z, rdim, nz, A->modulus);

            nmod_poly_mat_mul(Z, A, NN);

        if (nmod_poly_mat_is_zero(Z) !=0) 
        {
            printf("\n Nullspace of dimension %ld, time: %f, time flint: %f, ratio: %f\n\n",nz,t,t2,(double) t2/t);
            return 1;
        }

    }
    else {

        printf("\n Nullspace of dimension %ld, time: %f, time flint: %f, ratio: %f\n\n",nz,t,t2,t2/t);

    }


   

    // printf("Nflint output \n");
    // nmod_poly_mat_print_pretty(Nflint, "x");
    // printf("\n")
    

    nmod_poly_mat_clear(A);

    flint_rand_clear(state);

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
