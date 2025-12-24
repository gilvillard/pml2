
#include <stdio.h>   // pour FILE, fopen, fprintf, fclose
#include <string.h>
#include <stdlib.h>


#include <math.h> 
//#include <flint/profiler.h>
#include <time.h>
#include <flint/ulong_extras.h>

#include <pml/nmod_poly_mat_utils.h>
#include <pml/nmod_poly_mat_io.h>
//#include "../../nmod_poly_mat_extras/test/testing_collection.h"

#include "nmod_poly_mat_kernel.h"


void nmod_poly_mat_fprint_pretty(FILE *file, const nmod_poly_mat_t mat, const char * var)
{
    slong rdim = mat->r, cdim = mat->c;

    flint_fprintf(file, "<%wd x %wd matrix over Z/nZ[%s]>\n", mat->r, mat->c, var);
    flint_fprintf(file, "[");
    for (slong i = 0; i < rdim; i++)
    {
        flint_fprintf(file, "[");
        for (slong j = 0; j < cdim; j++)
        {
            nmod_poly_fprint_pretty(file,nmod_poly_mat_entry(mat, i, j), var);
            if (j+1 < cdim)
                flint_fprintf(file,", ");
        }
        if (i != rdim -1)
            flint_fprintf(file,"],\n");
        else
            flint_fprintf(file,"]");
    }
    flint_fprintf(file,"]\n");
}


int main(int argc, char ** argv)
{

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    int i,j;

    slong prime = 2; 


    nmod_poly_mat_t A;

    slong m=90;
    slong n=120;
    slong deg=12;


    slong rkflint;
    nmod_poly_mat_t N; 
    nmod_poly_mat_init(N, n, n, prime);

    slong nz;
    slong degN[n];

    slong verif;

    char namef[50];

    FILE* file;

    int nbfile = 0;

    
    printf("Launching  with\n\tprime = %ld,\n\tm = %ld,\n\tm = %ld,\
        \n\tdegree = %ld ...\n",prime,m, n, deg);

    // NB runs
    slong K=1000;

    for (int k=0; k<K; k++) {


    printf("\n=== %d ===\n",k);


        nmod_poly_mat_init(A, m,n , prime);

        nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.3);

        sprintf(namef, "resultat_%d.txt", nbfile);

        file = fopen(namef, "w");

        flint_fprintf(file,"A input \n");
        nmod_poly_mat_fprint_pretty(file, A, "x");

        fprintf(file,"\n");
        fclose(file);


        rkflint=nmod_poly_mat_rank(A);
    
        nz=nmod_poly_mat_kernel(N, degN, A, NULL, 2.1);

        verif =  (n-rkflint == nz);

        if (nz !=0) {

            nmod_poly_mat_t NN; 
            nmod_poly_mat_init(NN, n, nz, prime);

            for (i = 0; i < n; i++) {
                for (j = 0; j < nz; j++) {
                    nmod_poly_set(nmod_poly_mat_entry(NN, i, j), nmod_poly_mat_entry(N, i, j));
                }
            }

            nmod_poly_mat_t Z;
            nmod_poly_mat_init(Z, m, nz, prime);

            nmod_poly_mat_mul(Z, A, NN);

            verif = verif && nmod_poly_mat_is_zero(Z); 

            int is_weak_popov;

            is_weak_popov=nmod_poly_mat_is_ordered_weak_popov(NN, NULL, COL_UPPER);

            verif = verif && is_weak_popov;  


            printf("\n %ld  %ld\n",nz,verif);

            nmod_poly_mat_clear(Z);
            nmod_poly_mat_clear(NN);

        }


        if (!verif) {
            printf("Failed with\n\tprime = %ld,\n\tm = %ld,\n\tm = %ld,\
                n\tdegree = %ld ...\n",prime,m, n, deg);

            printf("   nz= %ld, nullity= %ld\n",nz,n-rkflint);

            printf("nbfile %d, k: %d\n\n",nbfile,k);

            nbfile+=1;

            k=K+1;
        }


    } // End of main loop 

    

    return 0;

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
