#include <stdlib.h>
#include <math.h> 
#include <flint/nmod.h> 
#include <flint/nmod_poly.h> 


#include "nmod_poly_mat_extra.h"
#include "nmod_poly_mat_io.h"

#include "nmod_mat_extra.h"
#include "nmod_poly_mat_kernel.h"


/**
 *  ZLS 
 * 
 *  N initialized inside ?
 * 
 */


int nmod_poly_mat_kernel(nmod_poly_mat_t N, const nmod_poly_mat_t M)
{
    
    printf("\n");
    nmod_poly_mat_print_pretty(M, "x");
    printf("\n");



    // slong i,j;

    // slong n = M->r;
    // slong m = M-> c;


    // nmod_poly_mat_t A;
    // nmod_poly_mat_init(A, n, n, M->modulus);

    // nmod_poly_mat_t B;
    // nmod_poly_mat_init(B, n, m-n, M->modulus);


    // for (i = 0; i < n; i++)
    // {
    //     for (j = 0; j < n; j++)
    //         nmod_poly_set(nmod_poly_mat_entry(A, i, j), nmod_poly_mat_entry(M, i, j));
    // }

    // for (i = 0; i < n; i++)
    // {
    //     for (j = 0; j < m-n; j++)
    //         nmod_poly_set(nmod_poly_mat_entry(B, i, j), nmod_poly_mat_entry(M, i, n+j));
    // }       

    // slong sigma;
    // sigma = ceil((double) m*delta/n +1);   // sigma >= ceil((double) (rdim + Bcdim)*delta/rdim +1);

    // nmod_poly_mat_t X;
    // nmod_poly_mat_init(X, n, m-n, M->modulus);

    // slong order = nmod_poly_mat_degree(A);

    // nmod_poly_mat_dixon(X, A, B, order, sigma);    

    // // printf("\n");
    // // nmod_poly_mat_print_pretty(A, "x");
    // // printf("\n");
    // // nmod_poly_mat_print_pretty(B, "x");
    // // printf("sigma %ld\n",sigma); 
    // // printf("\n");


    // nmod_poly_mat_t P;
    // nmod_poly_mat_init(P, n, m-n, A->modulus);

    // nmod_poly_mat_t Q;
    // nmod_poly_mat_init(Q, m-n, m-n, A->modulus);

    // slong nbnull;
    // nbnull=nmod_poly_mat_right_description(P, Q, X, delta);


    // if (nbnull != 0) 
    // {
    //     nmod_poly_mat_init(N, m, nbnull, A->modulus);


    //     for (i = 0; i < n; i++)
    //     {
    //         for (j = 0; j < nbnull; j++)
    //             nmod_poly_set(nmod_poly_mat_entry(N, i, j), nmod_poly_mat_entry(P, i, j));
    //     }


    //     for (i = 0; i < m-n; i++)
    //     {
    //         for (j = 0; j < nbnull; j++)
    //         {
    //             nmod_poly_set(nmod_poly_mat_entry(N, n+i, j), nmod_poly_mat_entry(Q, i, j));
    //             nmod_poly_neg(nmod_poly_mat_entry(N, n+i, j), nmod_poly_mat_entry(N, n+i, j));
    //         }
    //     }      

    //     nmod_poly_mat_t T;
    //     nmod_poly_mat_init(T, n, nbnull, M->modulus);

    //     nmod_poly_mat_mul(T,M,N);

    //     // printf("--------- %ld \n",nbnull);
    //     // nmod_poly_mat_print_pretty(N, "x");
    //     // printf("\n");

    //     if (nmod_poly_mat_is_zero(T) == 0)
    //     {
    //         nmod_poly_mat_clear(N); 
    //         nbnull = 0; 
    //     }

    //     nmod_poly_mat_clear(T);
    // }
       

    // nmod_poly_mat_clear(A);
    // nmod_poly_mat_clear(B);
    // nmod_poly_mat_clear(X);
    // nmod_poly_mat_clear(P);
    // nmod_poly_mat_clear(Q);

    return 0;
}







/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
