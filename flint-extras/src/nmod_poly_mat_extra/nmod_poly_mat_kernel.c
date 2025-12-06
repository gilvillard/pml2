#include <stdlib.h>
#include <math.h> 
#include <flint/nmod.h> 
#include <flint/nmod_poly.h> 


#include "nmod_poly_mat_extra.h"
#include "nmod_poly_mat_io.h"

#include "nmod_mat_extra.h"
#include "nmod_poly_mat_kernel.h"




/**
 * 
 * In place for the moment 
 * 
 * Couplé : calcule le shifted degree et tri : returns sdeg, initialized outside
 * 
 */
 

void sortM(nmod_poly_mat_t M, slong *sdeg, const slong *ishift)
{

    slong n = M->c;
    slong j;

    nmod_poly_mat_column_degree(sdeg, M, ishift);

    // Only required from the input matrix, not for the subsequent calls 
    for (j=0; j<n; j++) 
        if (sdeg[j] < 0) sdeg[j]=0;



    // printf("\n [ ");
    // for (int j=0; j<n-1; j++) 
    //     printf(" %ld, ",sdeg[j]);
    // printf(" %ld ]\n",sdeg[n-1]);


    slong * perm = flint_malloc(n * sizeof(slong));

    _nmod_poly_mat_permute_columns_by_sorting_vec(M, n, sdeg, perm);


    // printf("\n [ ");
    // for (int j=0; j<n-1; j++) 
    //     printf(" %ld, ",sdeg[j]);
    // printf(" %ld ]\n",sdeg[n-1]);

}


    //nmod_poly_mat_permute_columns();

/**
 *  ZLS 
 * 
 *  !!!  n >= m  (does not work otherwise, pb with degs/shift lengths)
 * 
 *  N initialized inside ?
 * 
 *  input ishift 
 *  output tshift, initialized outside
 * 
 *  returns the number of nullspace vectors
 * 
 */


int nmod_poly_mat_zls(nmod_poly_mat_t N, slong *tshift, const nmod_poly_mat_t A, const slong *ishift)
{

    slong i,j,k;

    slong m = A->r;
    slong n = A->c;


    slong kappa=2;
    slong min_mn;
    slong rho=0;

    if (m <= n) 
    {
        min_mn=m;
        for (i=n-m; i<n; i++)   // Pb null shift 
            rho+=ishift[i];   
    }
    else 
    {
        min_mn=n;

        for (i=0; i<n; i++)   // Pb null shift 
            rho+=ishift[i];  
        printf("\n m > n  %ld  %ld",m,n);   // To see, rectangular case 
    }  

    slong s;
    s = ceil((double) rho/min_mn); 

    printf("\n rho: %ld    s: %ld \n",rho, s);

    nmod_poly_mat_t AT;
    nmod_poly_mat_init(AT, n, m, A->modulus);
    nmod_poly_mat_transpose(AT,A);

    nmod_poly_mat_t PT;
    nmod_poly_mat_init(PT, n, n, A->modulus);

   // shift is modified in place 
    slong shift[n];
    for (i=0; i<n; i++)
    {
        shift[i]=ishift[i];
    }
    
    nmod_poly_mat_pmbasis(PT, shift, AT, kappa*s+1);

    nmod_poly_mat_t P;
    nmod_poly_mat_init(P, n, n, A->modulus);
    nmod_poly_mat_transpose(P,PT);


    // printf("\n");
    // nmod_poly_mat_print_pretty(P, "x");
    // printf("\n");


    // printf("\n [ ");
    // for (i=0; i<n-1; i++) 
    //     printf(" %ld, ",shift[i]);
    // printf(" %ld ]\n",shift[i]);


    nmod_poly_mat_t R;
    nmod_poly_mat_init(R, m, n, A->modulus);

    nmod_poly_mat_mul(R,A,P);

    // printf("\n");
    // nmod_poly_mat_print_pretty(R, "x");
    // printf("\n");

    slong cdeg[n];

    slong zshift[m];
    for (i=0; i<m; i++)
    {
        zshift[i]=0;
    }

    nmod_poly_mat_column_degree(cdeg, R, zshift);


    // printf("\n [ ");
    // for (i=0; i<n-1; i++) 
    //     printf(" %ld, ",cdeg[i]);
    // printf(" %ld ]\n",cdeg[i]);

    slong n1=0;
    slong n2;

    for (j=0; j<n; j++) {
        if (cdeg[j]<0) 
            n1+=1;
    }


    nmod_poly_mat_t P1;

    if (n1>0) {
        nmod_poly_mat_init(P1, n, n1, A->modulus);

        k=0;
        for (j = 0; j < n; j++)
        {
            if (cdeg[j]<0) {

                for (i = 0; i < n; i++)
                    nmod_poly_set(nmod_poly_mat_entry(P1, i, k), nmod_poly_mat_entry(P, i, j));
                k+=1;
            }

        }
    }
    

    n2=n-n1;

    if (n2==0) {

        nmod_poly_mat_init_set(N,P1);
        nmod_poly_mat_column_degree(tshift, P1, ishift);

        return n1;
    }

    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, n, n2, A->modulus);

    k=0;
    for (j = 0; j < n; j++)
    {
        if (cdeg[j]>=0) {

            for (i = 0; i < n; i++)
                nmod_poly_set(nmod_poly_mat_entry(P2, i, k), nmod_poly_mat_entry(P, i, j));
            k+=1;
        }

    }


    if (m==1){   // Then n2=0 ? 

        if (n1==0) 
            return 0; 
        else {
            nmod_poly_mat_init_set(N,P1);
            nmod_poly_mat_column_degree(tshift, P1, ishift);

            return n1;
        }

    }


    // n2 <> 0 and m> 1
    // ----------------

    sortM(P2,tshift,ishift);

    for (i = 0; i < n2; i++) {
        tshift[i]=tshift[i]-kappa*s;
    }
    
    printf("\n [ ");
    for (i=0; i<n2-1; i++) 
        printf(" %ld, ",tshift[i]);
    printf(" %ld ]\n",tshift[i]);




///+++++++++++++++++++++++


   //P:=res[1];
   //shift:=res[2];


//slong shift[2*n];




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
