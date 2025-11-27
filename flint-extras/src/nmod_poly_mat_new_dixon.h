
#ifndef NMOD_POLY_MAT_DIXON_H
#define NMOD_POLY_MAT_DIXON_H

#include <flint/nmod_types.h>
#include <flint/nmod_poly_mat.h>



void nmod_poly_mat_inv_trunc(nmod_poly_mat_t S, 
                            const nmod_poly_mat_t A, 
                            ulong order);


void nmod_poly_mat_dixon(nmod_poly_mat_t X, 
                            const nmod_poly_mat_t A, 
                            const nmod_poly_mat_t B, 
                            ulong order,
                            ulong sigma);

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


