#include <stdlib.h>
#include <math.h> 
#include <time.h>

#include <flint/nmod_poly.h> 
#include <flint/nmod_mat.h> 

#include <nmod_poly_mat_io.h>

#include "nmod_poly_mat_dixon.h"
#include "nmod_poly_mat_description.h"


/**
 *  Left description computation for H(x) n x m in K(x) (power series)
 *    with target degree delta
 * 
 *  requires enough precision in input: ie at least (m+n)*(delta+1)/m +1
 * 
 *  returns nbrows and a partial (or full) description when 0 < nbrows <= n rows, 
 *    or zero if no candidates 
 * 
 *   N  n x m and and D n x n  are initialized outside 
 *    hence the result is part of N and D  
 * 
 *  Should work for constant matrices 
 * 
 */

/**
 * n x m 
 */

slong guess_left_sigma(slong n, slong m, slong delta)
{
    return ceil((double) (m+n)*(delta+1)/m +1);
}


slong nmod_poly_mat_left_description(nmod_poly_mat_t D, nmod_poly_mat_t N,
                            const nmod_poly_mat_t H, 
                            slong delta)

{

    slong i,j;

    slong n = H->r;
    slong m = H->c;

    slong sigma;
    sigma = guess_left_sigma(n, m, delta);

    nmod_poly_mat_t M;
    nmod_poly_mat_init(M, n+m, m, H->modulus);

    nmod_poly_mat_set(M, H); // GV, correct even the dimension of M is larger ? 

    nmod_poly_t mone;
    nmod_poly_init(mone, H->modulus);
    nmod_poly_set_coeff_ui(mone, 0, H->modulus -1); // 1

    for (i = 0; i < m; i++)
        nmod_poly_set(nmod_poly_mat_entry(M, n+i, i), mone);

    nmod_poly_mat_t B;
    nmod_poly_mat_init(B, n+m, n+m, H->modulus);

    // Appropriate shift
    slong shift[n+m];

    for (i = 0; i < n+m; i++) 
        shift[i]=0; 

    // Shifted approximant computation 

    nmod_poly_mat_pmbasis(B, shift, M, sigma); 

    slong rows[n+m];
    slong nbrows=0;

    for (i = 0; i < n+m; i++) 
    {
        if (shift[i] <= delta) 
        {
            rows[nbrows]=i;
            nbrows +=1;           
        }
    }

    // See how to do depending on the level of informations wanted 
    // if (nbrows != n) 
    // {
    //     flint_printf("A complete description of degree at most %{slong} probably doesn't exist\n", delta);
    // }

    if (nbrows == 0) 
        return 0; 
    else 
    {
        for (i = 0; i < nbrows; i++)
        {
            for (j = 0; j < m; j++)
                nmod_poly_set(nmod_poly_mat_entry(N, i, j), nmod_poly_mat_entry(B, rows[i], n+j));
        }


        for (i = 0; i < nbrows; i++)
        {
            for (j = 0; j < n; j++)
                nmod_poly_set(nmod_poly_mat_entry(D, i, j), nmod_poly_mat_entry(B, rows[i], j));
        }
    }

    nmod_poly_mat_clear(M); 
    nmod_poly_mat_clear(B);

    nmod_poly_clear(mone); 

    return nbrows;
}


/**
 *  Right description computation for H(x) m x n in K(x) (power series)
 *    with target degree delta
 * 
 *  requires enough precision in input: ie at least (m+n)*(delta+1)/m +1
 * 
 *  returns nbcols and a partial (or full) description when 0 < nbcols <= n rows, 
 *    or zero if no candidates 
 * 
 *    N m x n and and D n x n  are initialized outside  
 *       hence the result is part of N and D
 * 
 *  Should work for constant matrices 
 * 
 */

/**
 * n x m 
 */

slong guess_right_sigma(slong n, slong m, slong delta)
{
    return ceil((double) (m+n)*(delta+1)/n +1);
}


slong nmod_poly_mat_description(nmod_poly_mat_t N, nmod_poly_mat_t D,
                            const nmod_poly_mat_t H, 
                            slong delta)

{
    slong nbcols;

    nmod_poly_mat_t NT;
    nmod_poly_mat_init(NT, H->c, H->r, H->modulus);

    nmod_poly_mat_t DT;
    nmod_poly_mat_init(DT, H->c, H->c, H->modulus);

    nmod_poly_mat_t HT;
    nmod_poly_mat_init(HT, H->c, H->r, H->modulus);


    nmod_poly_mat_transpose(HT,H);

    nbcols = nmod_poly_mat_left_description(DT,NT,HT,delta);

    if (nbcols == 0) 
        return nbcols;

    nmod_poly_mat_transpose(N,NT);
    nmod_poly_mat_transpose(D,DT);

    nmod_poly_mat_clear(NT); 
    nmod_poly_mat_clear(DT);
    nmod_poly_mat_clear(HT);

    return nbcols;

}


/**
 * 
 *  more columns m > n 
 * 
 *  Right kernel computation for M(x) n x m in K[x]
 *   with target degree delta
 * 
 *  N is intitialized inside for correct dimensions, n x nbnull
 *  returns the number nbnull and N, or zero (N not initialized)   
 *  
 *  Should work for constant matrices 
 * 
*/

int nmod_poly_mat_kernel(nmod_poly_mat_t N, const nmod_poly_mat_t M, slong delta)
{
    
    // printf("\n");
    // nmod_poly_mat_print_pretty(M, "x");
    // printf("\n");

    slong i,j;

    slong n = M->r;
    slong m = M-> c;

    if (m <= n) 
    {
            printf("error in nmod_poly_mat_kernel, the column dimension must be greater than the row one\n");
            exit(-1);
    }

    nmod_poly_mat_t A;
    nmod_poly_mat_init(A, n, n, M->modulus);

    nmod_poly_mat_t B;
    nmod_poly_mat_init(B, n, m-n, M->modulus);


    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(A, i, j), nmod_poly_mat_entry(M, i, j));
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m-n; j++)
            nmod_poly_set(nmod_poly_mat_entry(B, i, j), nmod_poly_mat_entry(M, i, n+j));
    }       


    nmod_poly_mat_t X;
    nmod_poly_mat_init(X, n, m-n, M->modulus);


    slong order=1;
    
    // TODO, which best order ?
    if (nmod_poly_mat_degree(A) !=0)
        order  = nmod_poly_mat_degree(A);
    

    slong sigma;
    sigma = guess_right_sigma(n,m,delta);  

    nmod_poly_mat_dixon(X, A, B, order, sigma);   

    // printf("\n");
    // nmod_poly_mat_print_pretty(A, "x");
    // printf("\n");
    // nmod_poly_mat_print_pretty(B, "x");
    // printf("sigma %ld\n",sigma); 
    // printf("\n");


    nmod_poly_mat_t P;
    nmod_poly_mat_init(P, n, m-n, A->modulus);

    nmod_poly_mat_t Q;
    nmod_poly_mat_init(Q, m-n, m-n, A->modulus);

    slong nbnull;
    nbnull=nmod_poly_mat_description(P, Q, X, delta);


    if (nbnull != 0) 
    {
        nmod_poly_mat_init(N, m, nbnull, A->modulus);


        for (i = 0; i < n; i++)
        {
            for (j = 0; j < nbnull; j++)
                nmod_poly_set(nmod_poly_mat_entry(N, i, j), nmod_poly_mat_entry(P, i, j));
        }


        for (i = 0; i < m-n; i++)
        {
            for (j = 0; j < nbnull; j++)
            {
                nmod_poly_set(nmod_poly_mat_entry(N, n+i, j), nmod_poly_mat_entry(Q, i, j));
                nmod_poly_neg(nmod_poly_mat_entry(N, n+i, j), nmod_poly_mat_entry(N, n+i, j));
            }
        }      

        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, n, nbnull, M->modulus);

        nmod_poly_mat_mul(T,M,N);

        // printf("--------- %ld \n",nbnull);
        // nmod_poly_mat_print_pretty(N, "x");
        // printf("\n");

        if (nmod_poly_mat_is_zero(T) == 0)
        {
            nmod_poly_mat_clear(N); 
            nbnull = 0; 
        }

        nmod_poly_mat_clear(T);
    }
       

    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(X);
    nmod_poly_mat_clear(P);
    nmod_poly_mat_clear(Q);

    return nbnull;
}


int nmod_poly_mat_kernel2(nmod_poly_mat_t N, const nmod_poly_mat_t M_input, slong delta)
{
   

    int i,j;

    slong n = M_input->r;
    slong m = M_input-> c;


    nmod_mat_t U0;
    nmod_mat_init(U0, m, n, M_input->modulus);

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    nmod_mat_randtril(U0,state,2);

    nmod_poly_mat_t U;
    nmod_poly_mat_init(U, m, n, M_input->modulus);

    //nmod_poly_mat_set_nmod_mat(U,U0); // GV does not work, bug ? Sam 14 jui 2025 18:01:54 CEST

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (nmod_mat_entry(U0, i, j) == 0)
                nmod_poly_zero(nmod_poly_mat_entry(U, i, j));
            else
            {
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(U, i, j),0,nmod_mat_entry(U0, i, j));
            }
        }
    }

    nmod_poly_mat_t M;
    nmod_poly_mat_init(M, n, m, M_input->modulus);

    nmod_poly_mat_mul(M,M_input,U);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m-n; j++)
        {
            nmod_poly_set(nmod_poly_mat_entry(M, i, n+j),nmod_poly_mat_entry(M_input, i, n+j));
            
        }
    }

    //++++++++++
            nmod_poly_mat_t T;
            nmod_poly_mat_init(T, M->r, M->r, M->modulus);

            for (int ii = 0; ii < n; ii++)
            {
                for (int j = 0; j < n; j++)
                nmod_poly_set(nmod_poly_mat_entry(T, ii, j), nmod_poly_mat_entry(M, ii, j));
            }


            nmod_poly_mat_truncate(T,2);

            nmod_poly_t det;
            nmod_poly_init(det, M->modulus);
            nmod_poly_mat_det(det, T);

            printf("\n");
            printf("DET :");
            nmod_poly_print_pretty(det,"x");
            printf("\n");
            printf("\n");
            
            nmod_poly_mat_clear(T);
    //+++++++++++++=

    slong nbnull;

    nbnull = nmod_poly_mat_kernel(N, M, delta);

    // nbnull not zero           

    nmod_poly_mat_t UH;
    nmod_poly_mat_init(UH, n, n, M_input->modulus);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(UH, i, j),nmod_poly_mat_entry(U, i, j));
            
    nmod_poly_mat_t NH;
    nmod_poly_mat_init(NH, n, nbnull, M_input->modulus);
    for (i = 0; i < n; i++)
        for (j = 0; j < nbnull; j++)
            nmod_poly_set(nmod_poly_mat_entry(NH, i, j),nmod_poly_mat_entry(N, i, j));
                    
    nmod_poly_mat_t TMP1;
    nmod_poly_mat_init(TMP1, n, nbnull, M_input->modulus);
    nmod_poly_mat_mul(TMP1, UH, NH);


    nmod_poly_mat_t UL;
    nmod_poly_mat_init(UL, m-n, n, M_input->modulus);
    for (i = 0; i < m-n; i++)
        for (j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(UL, i, j),nmod_poly_mat_entry(U, n+i, j));
            
    nmod_poly_mat_t NL;
    nmod_poly_mat_init(NL, m-n, nbnull, M_input->modulus);
    for (i = 0; i < m-n; i++)
        for (j = 0; j < nbnull; j++)
            nmod_poly_set(nmod_poly_mat_entry(NL, i, j),nmod_poly_mat_entry(N, n+i, j));

    nmod_poly_mat_t TMP2;
    nmod_poly_mat_init(TMP2, m-n, nbnull, M_input->modulus);
    nmod_poly_mat_mul(TMP2, UL, NH);
    nmod_poly_mat_add(TMP2, TMP2, NL);


    for (i = 0; i < n; i++)
        for (j = 0; j < nbnull; j++)
            nmod_poly_set(nmod_poly_mat_entry(N, i, j),nmod_poly_mat_entry(TMP1, i, j));

    for (i = 0; i < m-n; i++)
        for (j = 0; j < nbnull; j++)
            nmod_poly_set(nmod_poly_mat_entry(N, n+i, j),nmod_poly_mat_entry(TMP2, i, j));


 nmod_poly_mat_t Z;
    nmod_poly_mat_init(Z, n, nbnull, M_input->modulus);

nmod_poly_mat_mul(Z, M_input, N);

    printf("\n");
    nmod_poly_mat_print_pretty(Z, "x");
    printf("\n");




    return nbnull;
    // return, changer pour N à l'extérieur 

}





/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
