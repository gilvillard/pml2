#include <stdlib.h>
#include <math.h> 

#include <time.h>
#include <flint/ulong_extras.h>
#include <flint/profiler.h>

#include <pml/nmod_poly_mat_utils.h>
#include <pml/nmod_poly_mat_io.h>

#include "nmod_poly_mat_kernel.h"


// !!!!! cdim > rdim 

void benchmark_zls(ulong prime, slong rdim, slong cdim, ulong deg, double sparse, slong kappa, \
                    slong threshold, flint_rand_t state)
{
    

    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, rdim, cdim, prime);

    nmod_poly_mat_randtest_sparse(A, state, deg+1, sparse);


    nmod_poly_mat_t N; // Not initialized for the moment 

    slong i;
    slong iz[rdim];

    
    slong degN[cdim];

    for (i = 0; i < rdim; i++) 
        iz[i]=0; 


    nmod_poly_mat_t Nflint;
    nmod_poly_mat_init(Nflint, cdim, cdim, A->modulus);

    //slong nz;
    //slong rkflint;

    // for timing
    //long thres = 100000; // 1000ms = 1s
    long nb_iter;
    long t_pml = 0;  //time, in ms
    long t_flint = 0; // time, in ms
    timeit_t timer1,timer2;

    flint_printf("%ld\t%ld\t%ld\t|\t", rdim, cdim, deg);

    nb_iter = 4; t_pml = 0; t_flint = 0;

    for (int i=0; i<nb_iter; i++)
    {
        timeit_start(timer1);
        nmod_poly_mat_zls(N, degN, A, iz, kappa,threshold);
        timeit_stop(timer1);
        t_pml += timer1->wall;

        timeit_start(timer2);
        nmod_poly_mat_nullspace(Nflint,A),
        //nmod_poly_mat_rank(A);

        timeit_stop(timer2);
        t_flint += timer2->wall;
    }

    double pml = (double)t_pml/nb_iter;
    double flint = (double)t_flint/nb_iter;

    //flint_printf("%0.1e\t%0.1e\t\t%f\t",
    flint_printf("%.1f\t%.1f\t|\t%.1f\t",
       (double)pml,
       (double)flint,
       (flint/pml));


    flint_printf("\n");
}



int main(int argc, char ** argv)
{

    setlinebuf(stdout);
    srand(time(NULL));
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, rand(), rand());


    //printf("Usage: %s OR %s [nbits] [rdim] [cdim] [degree] [sparsity] [kappa] \n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    // slong nbits = atoi(argv[1]);
    // slong rdim = atoi(argv[2]);
    // slong cdim = atoi(argv[3]);
    // slong order = atoi(argv[4]);
    // slong kappa = atoi(argv[5]);

    // flint_rand_t state;
    // flint_rand_init(state);
    // srand(time(NULL));
    // flint_rand_set_seed(state, rand(), rand());

    int m;

    for (m=60; m<80; m++) {
        benchmark_zls(104729, m, m+10, 4, 1.0, 2, 2, state);
    }



    flint_rand_clear(state);

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
