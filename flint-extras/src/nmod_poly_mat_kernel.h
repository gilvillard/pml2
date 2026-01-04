/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_KERNEL_H
#define NMOD_POLY_MAT_KERNEL_H

#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  
 *  Right kernel of a polynomial matrix, in shift-ordered weak Popov form
 *   
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 * 
 *  TODO/TO SEE: ....
 *    
 * Input: 
 *     - A in m x n 
 *     - input_shift[n], NULL (the degrees are computed) or initialized outside, 
 *         the shift for the output kernel    
 *     - kappa, a double >= 2, for the order of the approximant bases that are used 
 *         i.e. we use kappa * s instead of 3s in ZLS  
 *
 *  Output:
 *    - returns the dimension w of the kernel, which may be zero 
 *    - N is initialized n x w by the procedure if w >0,  
 *       its w columns give a minimal basis of the kernel   
 *    - degN[n], initialized outside, its first w entries are concerned,
 *        they are the shifted degrees of the kernel basis 
 * 
 */

int nmod_poly_mat_kernel(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                         const slong *ishift, const double kappa); 



/**
 * Experimental, should not be really considered  
 *
 */

int nmod_poly_mat_approximant_kernel(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                                 const slong *ishift);


#ifdef __cplusplus
}
#endif

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


