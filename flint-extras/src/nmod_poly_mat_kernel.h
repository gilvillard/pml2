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

#include <flint/nmod_types.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif


void sortM(nmod_poly_mat_t M, slong *sdeg, const slong *ishift);

int nmod_poly_mat_zls(nmod_poly_mat_t N, slong *tshift, const nmod_poly_mat_t M, const slong *ishift); 


#ifdef __cplusplus
}
#endif

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


