#include <math.h> 
#include <stdlib.h>

#include "nmod_poly_mat_extra.h"
#include "nmod_poly_mat_kernel.h"

/**
 * 
 * The content of the file has not yet been tested very extensively
 * 
 */ 

/**
 * 
 * In place, M is m x n, shifted column degrees 
 *   permute the columns in nondecreasing order 
 * 
 * Todo/to see: zero is set to degree 0 (not -1), for specific use 
 *     in ZLS algorithm for the kernel
 * 
 * Input:
 *  M  m x n
 *  ishift[m]
 * 
 * Output:
 *  M is modified in place 
 *  perm[n],initialized outside, the carried out permutation 
 *  sdeg[n] initialized outside
 * 
 */
 

// Dummy function for qsort 
 int cmp(const void *a, const void *b)
{
    const slong *x = a;
    const slong *y = b;

    if (*x < *y) return -1;
    if (*x > *y) return 1;
    return 0;
}


/**
 *  
 *  Internal function 
 * 
 *  Right shifted kernel of a polynomial matrix, assuming that the columns has been 
 *     sorted by shifted degree 
 * 
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 * 
 *  TODO/TO SEE: 
 *    
 *  Input:
 *    A in m x n 
 *    ishift[n], NULL (the degrees are computed) or initialized outside, 
 *      the shift for the kernel   
 *      values should be at least 0 (even for zero columns in A
 *      "with entries arranged in non-decreasing order and bounding the 
 *       corresponding column degrees of A." 
 *    kappa, a double >= 2, for the order of the order bases 
 *              kappa * s instead of 3 *s in ZLS  
 * 
 *  Output:
 *    returns the dimension w of the kernel, which may be zero 
 *    N, n x w, is initialized here only when w > 0, 
 *           should be freed outside in that case only,
 *            gives a minimal basis of the right kernel   
 *    degN[n], initialized outside, its first w entries are concerned,
 *        they are the ishift shift degrees of the kernel basis 
 * 
 */


/**
 *  
 *  Right shifted kernel of a polynomial matrix, assuming that the columns has been 
 *     sorted by shifted degree 
 * 
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 * 
 *  Calls nmod_poly_mat_zls_sorted after an initial sorting 
 * 
 *  TODO/TO SEE: 
 *    
 * Input: 
 *    iA in m x n 
 *     ishift[n], NULL (the degrees are computed) or initialized outside, 
 *      the shift for the kernel    
 *      values should be at least 0 (even for zero columns in A
 *      "with entries arranged in non-decreasing order and bounding the 
 *       corresponding column degrees of A." 
 *    kappa, a double >= 2, for the order of the order bases 
 *              kappa * s instead of 3 *s in ZLS  
 *
 *  Output:
 *    returns the dimension w of the kernel, which may be zero 
 *    N, is initialized  n x n outside  
 *       its first w columns give a minimal basis of the kernel   
 *    degN[n], initialized outside, its first w entries are concerned,
 *        they are the ishift shifted degrees of the kernel basis 
 * 
 */

int nmod_poly_mat_kernel(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                                 const slong *input_shift, const double kappa)
{

    slong i,j,k;

    slong m = A->r;
    slong n = A->c;


    //  Input shift modified to ishift, to meet ZLS requirements,  
    //    internal ishift bounding the column degrees 
    //    we add delta to input_shift 
    // ----------------------------------------------------------

    slong sdeg[n];

    nmod_poly_mat_column_degree(sdeg, A, NULL);
    
    for (j=0; j<n; j++) {  // to prevent values -1 with zero columns in the shift 
        if (sdeg[j] < 0) sdeg[j]=0;
    }

    slong ishift[n];
    slong delta;

    if (input_shift == NULL) {

        delta=sdeg[0]; // will be the maximum degree,    
        for (j=1; j<n; j++) {
            if (sdeg[j] > delta) delta=sdeg[j];
        }

        for (j=0; j<n; j++) {
            ishift[j]=delta;
        }
    }
    else {

        delta=sdeg[0]-input_shift[0];

        for (j=1; j<n; j++) {
            if ((sdeg[j]-input_shift[j]) > delta) delta=sdeg[j]-input_shift[j];
        }
  
        for (j=0; j<n; j++) {
            ishift[j]=input_shift[j]+delta;
        }
    }

    // Tuning the order of the subsequent approximant 
    // ----------------------------------------------

    slong min_mn;
    slong rho=0;

    long tmp_shift[n]; // temporary variable, since sometimes the shift is modified inplace 

    // In order to sort for computing the order of the approximant, for the test m=1 to be ok 
    for (i=0; i<n; i++)
    {
        tmp_shift[i]=ishift[i];
    }

    qsort(tmp_shift, n, sizeof(slong), cmp);

    if (m <= n) 
    {
        min_mn=m;
        for (i=n-m; i<n; i++)   
            rho+=tmp_shift[i];   
    }
    else 
    {   // No sort needed, we consider all degrees 
        min_mn=n;
        for (i=0; i<n; i++)   
            rho+=ishift[i];  
    }  


    slong s;
    s = ceil((double) rho/min_mn); 

    // Transposed A for the approximant PT on the left 
    // PT will then be use without transposing 

    nmod_poly_mat_t AT;
    nmod_poly_mat_init(AT, n, m, A->modulus);
    nmod_poly_mat_transpose(AT,A);

    nmod_poly_mat_t PT;
    nmod_poly_mat_init(PT, n, n, A->modulus);

    // Approximant PT 
    // tmp_shift is modified in place 
    // --------------------------

    for (i=0; i<n; i++)
    {
        tmp_shift[i]=ishift[i];
    }
    
    slong ks;
    ks=ceil((double) kappa*s); 
    nmod_poly_mat_pmbasis(PT, tmp_shift, AT, ks+1);


    // Looking for zero residues and non zero residues
    //  the global residue is reused later 
    // -----------------------------------------------

    nmod_poly_mat_t RT;
    nmod_poly_mat_init(RT, n, m, A->modulus);
    nmod_poly_mat_mul(RT,PT,AT);
    
    slong cdeg[n];
    nmod_poly_mat_row_degree(cdeg, RT, NULL);

    nmod_poly_mat_clear(AT);

    // Zero residue matrix P1  n x n1 
    // nonzero residues matrix P2  n x n2
    //-----------------------------------

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
                    nmod_poly_set(nmod_poly_mat_entry(P1, i, k), nmod_poly_mat_entry(PT, j, i));
                k+=1;
            }
        }
    }
    
    //++++++++
    slong nb=0;

    for (int j=0; j<n; j++) {

        if (tmp_shift[j] < (ks+1))
            nb+=1;
    }

    if (nb != n1) {
        printf("*******   nb: %ld  n1: %ld \n\n",nb,n1);

        printf("sigma: %ld  \n [ ",ks+1);
        for (j=0; j<n-1; j++) 
            printf(" %ld, ",tmp_shift[j]);
        printf(" %ld ]\n\n",tmp_shift[n-1]);
    }


    n2=n-n1;

    // the kernel is found 
    if (n2==0) {
        nmod_poly_mat_init_set(N,P1);
        nmod_poly_mat_column_degree(degN, P1, input_shift); ///

        nmod_poly_mat_clear(P1);
        return n1;
    }

    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, n, n2, A->modulus);

    k=0;
    for (j = 0; j < n; j++)
    {
        if (cdeg[j]>=0) {

            for (i = 0; i < n; i++)
                nmod_poly_set(nmod_poly_mat_entry(P2, i, k), nmod_poly_mat_entry(PT, j, i));
            k+=1;
        }

    }

    //nmod_poly_mat_clear(PT);

    //  Special case m=1
    //  before the divide and conquer on m
    //  the kernel is found 
    // -----------------------------------

    if (m==1){  

        if (n1==0) 
            return 0; 
        else {

            nmod_poly_mat_init_set(N,P1);
            nmod_poly_mat_column_degree(degN, P1, input_shift); ///

            nmod_poly_mat_clear(P1);
            return n1;
        }
    }

    //  Now n2 <> 0 and m > 1
    //    one can proceed to the divide and conquer from the rows 
    //    of the nonzero residue P2, and of A.P2
    // --------------------------------------------------------

    slong degP2[n2];

    nmod_poly_mat_column_degree(degP2, P2, ishift);

    for (i = 0; i < n2; i++) {
        degP2[i]=degP2[i]-ks;  // used below for the recursive calls 
    }

    nmod_poly_mat_t TT;
    nmod_poly_mat_init(TT, m, n2, A->modulus);

    // We extract G (temporary TT) from RT as we have been extracting P2 from PT
    //   and apply above ordering of P2 
    k=0;
    for (j = 0; j < n; j++)
    {
        if (cdeg[j]>=0) {

            for (i = 0; i < m; i++)
                nmod_poly_set(nmod_poly_mat_entry(TT, i, k), nmod_poly_mat_entry(RT, j, i));
            k+=1;
        }

    }

    nmod_poly_mat_clear(RT);

    nmod_poly_mat_t G;
    nmod_poly_mat_init(G, m, n2, A->modulus);

    nmod_poly_mat_shift_right(G,TT,ks);
    nmod_poly_mat_clear(TT);

    // We split G for the recursive call
    // ---------------------------------

    slong new_m=floor((double) m/2);

    nmod_poly_mat_t G1;
    nmod_poly_mat_init(G1, new_m, n2, A->modulus);

    for (i = 0; i < new_m; i++){
        for (j = 0; j < n2; j++) {
            nmod_poly_set(nmod_poly_mat_entry(G1, i, j), nmod_poly_mat_entry(G, i, j));
        }
    }

    nmod_poly_mat_t G2;
    nmod_poly_mat_init(G2, m-new_m, n2, A->modulus);

    for (i = 0; i < m-new_m; i++) {
        for (j = 0; j < n2; j++) {
            nmod_poly_set(nmod_poly_mat_entry(G2, i, j), nmod_poly_mat_entry(G, i+new_m, j));
        }
    }

    nmod_poly_mat_clear(G);
   

    // Recursive calls 
    // ---------------

    nmod_poly_mat_t N1;
    nmod_poly_mat_t N2;

    slong c1=0;
    slong c2=0;


    c1=nmod_poly_mat_kernel(N1, degN, G1, degP2, kappa); 


    if (c1 != 0) {
        
        nmod_poly_mat_t G3;
        nmod_poly_mat_init(G3, m-new_m, c1, A->modulus);

        nmod_poly_mat_mul(G3, G2, N1);

        for (i=0; i<c1; i++) {
            tmp_shift[i]=degN[i];
        }

        c2=nmod_poly_mat_kernel(N2, degN, G3, tmp_shift, kappa); 
        nmod_poly_mat_clear(G3);
    }


    nmod_poly_mat_clear(G1);
    nmod_poly_mat_clear(G2);

    // the recursive calls did not provide anything more 

    if ((c1==0) || (c2==0)){

        if (n1==0) {
            return 0; 
        }
        else {
            nmod_poly_mat_init_set(N,P1);
            nmod_poly_mat_column_degree(degN, P1, input_shift); ///

            nmod_poly_mat_clear(P1);
            nmod_poly_mat_clear(P2);
            return n1;
        }
    }

    // A new part of the kernel has been found by the recursive calls, Q
    // -----------------------------------------------------------------


    nmod_poly_mat_t Q1;
    nmod_poly_mat_init(Q1, n, c1, A->modulus);

    nmod_poly_mat_mul(Q1, P2, N1);

    nmod_poly_mat_t Q;
    nmod_poly_mat_init(Q, n, c2, A->modulus);

    nmod_poly_mat_mul(Q, Q1, N2);

    nmod_poly_mat_clear(N1);
    nmod_poly_mat_clear(N2);
    nmod_poly_mat_clear(Q1);    
    nmod_poly_mat_clear(P2);


    if (n1 ==0) {

        nmod_poly_mat_init_set(N,Q); // We should not need to copy 
        nmod_poly_mat_column_degree(degN, Q, input_shift); /// 

        nmod_poly_mat_clear(Q); 
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

        slong pivind_N[n1+c2];

        nmod_poly_mat_pivot_index(pivind_N,N,ishift,COL_UPPER);

        slong * perm = flint_malloc((n1+c2) * sizeof(slong));
        _nmod_poly_mat_permute_columns_by_sorting_vec(N, n1+c2, pivind_N, perm);

        nmod_poly_mat_column_degree(degN, N, input_shift); 

        nmod_poly_mat_clear(P1); 
        nmod_poly_mat_clear(Q); 
        
        return n1+c2;

    }

    return 0; 
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
