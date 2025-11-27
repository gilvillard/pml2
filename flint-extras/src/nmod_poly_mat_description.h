
#ifndef NMOD_POLY_MAT_DESCRIPTION_H
#define NMOD_POLY_MAT_DESCRIPTION_H

#include <flint/nmod_types.h>
#include <flint/nmod_poly_mat.h>



slong nmod_poly_mat_left_description(nmod_poly_mat_t N,  nmod_poly_mat_t D,
                            const nmod_poly_mat_t H, 
                            slong delta);


slong nmod_poly_mat_description(nmod_poly_mat_t N,  nmod_poly_mat_t D,
                            const nmod_poly_mat_t H, 
                            slong delta);


int nmod_poly_mat_kernel(nmod_poly_mat_t N, const nmod_poly_mat_t M, 
                            slong delta);


#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


