/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_vec.h>  // for _nmod_vec_set

#include "nmod_mat_extra.h"


/**********************************************************************
*                    ROW ROTATION DOWNWARD/UPWARD                    *
**********************************************************************/

void _nmod_mat_rotate_rows_downward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        ulong * tmp_mat = mat->rows[j];
        for (slong ii = j; ii > i; ii--)
            mat->rows[ii] = mat->rows[ii-1];
        mat->rows[i] = tmp_mat;
#else
        ulong * tmp_row = _nmod_vec_init(mat->c);
        _nmod_vec_set(tmp_row, nmod_mat_entry_ptr(mat, j, 0), mat->c);
        for (slong ii = j; ii > i; ii--)
            _nmod_vec_set(nmod_mat_entry_ptr(mat, ii, 0),
                          nmod_mat_entry_ptr(mat, ii-1, 0),
                          mat->c);
        _nmod_vec_set(nmod_mat_entry_ptr(mat, i, 0), tmp_row, mat->c);
        _nmod_vec_clear(tmp_row);
#endif
    }
}

void _nmod_mat_rotate_rows_upward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[i];
            for (slong ii = i; ii < j; ii++)
                vec[ii] = vec[ii+1];
            vec[j] = tmp_vec;
        }

#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        ulong * tmp_mat = mat->rows[i];
        for (slong ii = i; ii < j; ii++)
            mat->rows[ii] = mat->rows[ii+1];
        mat->rows[j] = tmp_mat;
#else
        ulong * tmp_row = _nmod_vec_init(mat->c);
        _nmod_vec_set(tmp_row, nmod_mat_entry_ptr(mat, i, 0), mat->c);
        for (slong ii = i; ii < j; ii++)
            _nmod_vec_set(nmod_mat_entry_ptr(mat, ii, 0),
                          nmod_mat_entry_ptr(mat, ii+1, 0),
                          mat->c);
        _nmod_vec_set(nmod_mat_entry_ptr(mat, j, 0), tmp_row, mat->c);
        _nmod_vec_clear(tmp_row);
#endif
    }
}


/**********************************************************************
*                COLUMN ROTATION RIGHTWARD/LEFTWARD                  *
**********************************************************************/

void _nmod_mat_rotate_columns_rightward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

        for (slong ii = 0; ii < mat->r; ii++)
        {
            ulong tmp_mat = nmod_mat_entry(mat, ii, j);
            for (slong jj = j; jj > i; jj--)
                nmod_mat_entry(mat, ii, jj) = nmod_mat_entry(mat, ii, jj-1);
            nmod_mat_entry(mat, ii, i) = tmp_mat;
        }
    }
}

void _nmod_mat_rotate_columns_leftward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[i];
            for (slong ii = i; ii < j; ii++)
                vec[ii] = vec[ii+1];
            vec[j] = tmp_vec;
        }

        for (slong ii = 0; ii < mat->r; ii++)
        {
            ulong tmp_mat = nmod_mat_entry(mat, ii, i);
            for (slong jj = i; jj < j; jj++)
                nmod_mat_entry(mat, ii, jj) = nmod_mat_entry(mat, ii, jj+1);
            nmod_mat_entry(mat, ii, i) = tmp_mat;
        }
    }
}

/**********************************************************************
*                         COLUMN PERMUTATION                         *
**********************************************************************/

/** Permute columns of a matrix `mat` according to `perm_act`, and propagate
 * the action on `perm_store`.
 * That is, performs for each appropriate index `j`, the operations
 * `perm_store[j] <- perm_store[perm_act[j]]`
 * `mat[i][j] <- mat[i][perm_act[j]] for all row indices i`
 **/
void nmod_mat_permute_columns(nmod_mat_t mat, const slong * perm_act, slong * perm_store)
{
    ulong * row_buffer = (ulong *) flint_malloc(mat->c * sizeof(ulong));

    /* perm_store[j] <- perm_store[perm_act[j]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    /* rows[i] <- rows[perm_act[i]]  */
    for (slong i = 0; i < mat->r; i++)
    {
        // copy row i into buffer
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        _nmod_vec_set(row_buffer, mat->rows[i], mat->c);
#else
        _nmod_vec_set(row_buffer, nmod_mat_entry_ptr(mat, i, 0), mat->c);
#endif
        // permute row i
        for (slong j = 0; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = row_buffer[perm_act[j]];
    }

    flint_free(row_buffer);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
