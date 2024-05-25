#include <flint/nmod_poly.h>

#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_sd_fft_init_set(nmod_sd_fft_t F, ulong w, ulong order, nmod_t mod)
{
    vec1d inv_2, inv, dp;
    ulong k;

    F->mod = mod;
    F->order = order;
    F->w = w;
    F->inv_w = nmod_inv(w, mod);
    F->p = mod.n;
    F->pinv = 1.0/F->p;
    
    F->powers_w = (vec1d **) flint_malloc(sizeof(vec1d *) * (order + 1));
    F->powers_inv_w_t = (ulong **) flint_malloc(sizeof(vec1d *) * (order + 1));

    dp = mod.n;

    for (k = 0; k <= order; k++)
    {
        ulong wr, inv_wr;
        ulong j, K;
        
        wr = w;
        inv_wr = F->inv_w;
        for (j = 0; j < order-k; j++)
        {
            wr = nmod_mul(wr, wr, mod);
            inv_wr = nmod_mul(inv_wr, inv_wr, mod);
        } 
        // wr = w^{2^(order-k)}
        // inv_wr = 1/w^{2^(order-k)}
        
        K = 1L << k;
        if ((sizeof(vec1d) * K) >= 32)
        {
            F->powers_w[k] = (vec1d *) aligned_alloc(32, sizeof(vec1d) * K);
            F->powers_inv_w_t[k] = (ulong *) aligned_alloc(32, sizeof(vec1d) * K);
        }
        else
        {
            F->powers_w[k] = (vec1d *) flint_malloc(sizeof(vec1d) * K);
            F->powers_inv_w_t[k] = (ulong *) flint_malloc(sizeof(vec1d) * K);
        }
        
        if (K == 1)
        {
            F->powers_w[k][0] = 1;
            F->powers_inv_w_t[k][0] = 1;
        }
        else
        {
            ulong nb;
            nb = 0;
            K = K/2;
            
            while (K >= 1)
            {
                ulong wri, inv_wri;
                ulong i;

                wri = 1;     // wr^i
                inv_wri = 1;

                for (i = 0; i < K; i++)
                {
                    F->powers_w[k][nb] = (vec1d) wri;
                    wri = nmod_mul(wri, wr, mod);
                    F->powers_inv_w_t[k][nb] = inv_wri;
                    inv_wri = nmod_mul(inv_wri, inv_wr, mod);
                    nb++;
                }

                wr = nmod_mul(wr, wr, mod);
                inv_wr = nmod_mul(inv_wr, inv_wr, mod);
                K /= 2;
            }
        }
    }


    F->powers_inv_2 = (nn_ptr) flint_malloc(sizeof(vec1d) * (order + 1));
    inv_2 = nmod_inv(2, mod);
    inv = 1;
    for (k = 0; k <= order; k++)
    {
        F->powers_inv_2[k] = inv;
        inv = nmod_mul(inv, inv_2, mod);
    }

    if (order > 0)
    {
        F->powers_inv_w_over_2 = (vec1d **) flint_malloc(sizeof(vec1d *) * order);
        for (k = 0; k < order; k++)
        {
            ulong i, K;
            ulong *src;
            
            K = 1L << k;
            
            if ((sizeof(vec1d) * K) >= 32)
            {
                F->powers_inv_w_over_2[k] = (vec1d *) aligned_alloc(32, sizeof(vec1d) * K);
            }
            else
            {
                F->powers_inv_w_over_2[k] = (vec1d *) flint_malloc(sizeof(vec1d) * K);
            }
            src = F->powers_inv_w_t[k+1];
            
            for (i = 0; i < K; i++)
            {
                F->powers_inv_w_over_2[k][i] = vec1d_reduce_0n_to_pmhn((vec1d) nmod_mul(src[i], F->powers_inv_2[k], mod), dp);
            }
        }
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
