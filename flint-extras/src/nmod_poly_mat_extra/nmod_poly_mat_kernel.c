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
 

void sortM(nmod_poly_mat_t M, slong *sdeg, slong *perm, const slong *ishift)
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


    //slong * perm = flint_malloc(n * sizeof(slong));

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


int nmod_poly_mat_zls_sorted(nmod_poly_mat_t N, slong *tshift, const nmod_poly_mat_t A, const slong *ishift)
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

        // for (i=0; i<n; i++)   // Pb null shift 
        //     rho+=ishift[i];  
        // printf("\n m > n  %ld  %ld",m,n);   // To see, rectangular case 
    }  

    slong s;
    s = ceil((double) rho/min_mn); 

    //printf("\n rho: %ld    s: %ld \n",rho, s);

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


    // !!!!!  CHECK
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

    slong * perm = flint_malloc(n2 * sizeof(slong));

    sortM(P2,tshift,perm,ishift);

    for (i = 0; i < n2; i++) {
        tshift[i]=tshift[i]-kappa*s;
    }
    
    // printf("\n [ ");
    // for (i=0; i<n2-1; i++) 
    //     printf(" %ld, ",tshift[i]);
    // printf(" %ld ]\n",tshift[n2-1]);


    nmod_poly_mat_t G;
    nmod_poly_mat_init(G, m, n2, A->modulus);

    nmod_poly_mat_mul(G, A, P2);

    nmod_poly_mat_t TT;
    nmod_poly_mat_init(TT, m, n2, A->modulus);


    nmod_poly_mat_shift_right(TT,G,kappa*s);

    slong new_m=floor((double) m/2);

    nmod_poly_mat_t G1;
    nmod_poly_mat_init(G1, new_m, n2, A->modulus);

    for (i = 0; i < new_m; i++){
        for (j = 0; j < n2; j++) {
            nmod_poly_set(nmod_poly_mat_entry(G1, i, j), nmod_poly_mat_entry(TT, i, j));
        }
    }

    nmod_poly_mat_t G2;
    nmod_poly_mat_init(G2, m-new_m, n2, A->modulus);

    for (i = 0; i < m-new_m; i++) {
        for (j = 0; j < n2; j++) {
            nmod_poly_set(nmod_poly_mat_entry(G2, i, j), nmod_poly_mat_entry(TT, i+new_m, j));
        }
    }

    // printf("--- G1 \n");
    // nmod_poly_mat_print_pretty(G1, "x");
    // printf("\n");

    // printf("--- G2 \n");
    // nmod_poly_mat_print_pretty(G2, "x");
    // printf("\n");

    // Change the dimension of the shift? tshift can be reused in output

    //+++++++++++++++++++++++++++++++++++++++++

    for (i=0; i<n2; i++) {
        shift[i]=tshift[i];
    }

    slong res;

    nmod_poly_mat_t N1;
    nmod_poly_mat_t N2;

    res=nmod_poly_mat_zls_sorted(N1, tshift, G1, shift); 

    slong c1 = N1->c;

    if (res != 0) {
        
        nmod_poly_mat_t G3;
        nmod_poly_mat_init(G3, m-new_m, c1, A->modulus);

        nmod_poly_mat_mul(G3, G2, N1);

        for (i=0; i<c1; i++) {
            shift[i]=tshift[i];
        }

        res=nmod_poly_mat_zls_sorted(N2, tshift, G3, shift); 

    }

    if (res==0) {

        if (n1==0) {
            return 0; 
        }
        else {
            nmod_poly_mat_init_set(N,P1);
            nmod_poly_mat_column_degree(tshift, P1, ishift);

            return n1;
        }

    }

    slong c2 = N2->c;

    nmod_poly_mat_t Q1;
    nmod_poly_mat_init(Q1, n, c1, A->modulus);

    nmod_poly_mat_mul(Q1, P2, N1);

    nmod_poly_mat_t Q;
    nmod_poly_mat_init(Q, n, c2, A->modulus);

    nmod_poly_mat_mul(Q, Q1, N2);

    if (n1 ==0) {

        nmod_poly_mat_init_set(N,Q);
        nmod_poly_mat_column_degree(tshift, Q, ishift);

        return c2;

    }
    else {

        nmod_poly_mat_init(N, n, n1+c2, A->modulus);

        for (i = 0; i < n; i++) {
            for (j = 0; j < n1; j++) {
                nmod_poly_set(nmod_poly_mat_entry(N, i, j), nmod_poly_mat_entry(P1,i,j));
            }
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < c2; j++) {
                nmod_poly_set(nmod_poly_mat_entry(N, i, j+n1), nmod_poly_mat_entry(Q, i, j));
            }
        }


        slong oshift[m]; // Too big, to see

        nmod_poly_mat_column_degree(oshift, P1, ishift);

        for (i=0; i<n1; i++) {
            tshift[i]=oshift[i];
        }

        nmod_poly_mat_column_degree(oshift, Q, ishift);

        for (i=0; i<c2; i++) {
            tshift[i+n1]=oshift[i];
        }

        return n1+c2;

    }

    // nmod_poly_mat_clear(A);
    // nmod_poly_mat_clear(B);
    // nmod_poly_mat_clear(X);
    // nmod_poly_mat_clear(P);
    // nmod_poly_mat_clear(Q);

    return 0; // ??? 
}


int nmod_poly_mat_zls(nmod_poly_mat_t N, slong *tshift, const nmod_poly_mat_t iA, const slong *ishift)
{

    slong i,j,k;

    slong m = iA->r;
    slong n = iA->c;


    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, m, n, iA->modulus);


    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            nmod_poly_set(nmod_poly_mat_entry(A, i, j), nmod_poly_mat_entry(iA,i,j));
        }
    }


    slong * perm = flint_malloc(n * sizeof(slong));

    slong shift[n];

    sortM(A,shift,perm,ishift);   

    slong res;

    nmod_poly_mat_t NT;
    res=nmod_poly_mat_zls_sorted(NT, tshift, A, shift);


    // printf("--- NT \n");
    // nmod_poly_mat_print_pretty(NT, "x");
    // printf("\n");

    // Back to A 
    if (res !=0) {

        slong nz = NT->c;

        nmod_poly_mat_init(N, n, nz, A->modulus);

        for (k = 0; k < n; k++) {
            for (j = 0; j < nz; j++){

                nmod_poly_set(nmod_poly_mat_entry(N, perm[k], j), nmod_poly_mat_entry(NT,k,j));
            }
        }

        return nz; 
    }
    
    return 0; 
    
}





/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
