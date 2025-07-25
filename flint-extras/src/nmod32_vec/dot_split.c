/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "flint/nmod.h"  // for NMOD_RED
#include "flint/nmod_vec.h"  // for DOT_SPLIT_MASK

#include "nmod32_vec.h"

uint _nmod32_vec_dot_split(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
{
    ulong dp_lo = 0;
    ulong dp_hi = 0;

    slong i = 0;

    for ( ; i+7 < len; i+=8)
    {
        dp_lo += (ulong)vec1[i+0] * vec2[i+0]
               + (ulong)vec1[i+1] * vec2[i+1]
               + (ulong)vec1[i+2] * vec2[i+2]
               + (ulong)vec1[i+3] * vec2[i+3]
               + (ulong)vec1[i+4] * vec2[i+4]
               + (ulong)vec1[i+5] * vec2[i+5]
               + (ulong)vec1[i+6] * vec2[i+6]
               + (ulong)vec1[i+7] * vec2[i+7];
        dp_hi += (dp_lo >> DOT_SPLIT_BITS);
        dp_lo &= DOT_SPLIT_MASK;
    }

    // less than 8 terms remaining, can be accumulated
    for ( ; i < len; i++)
        dp_lo += (ulong)vec1[i] * vec2[i];
    dp_hi += (dp_lo >> DOT_SPLIT_BITS);
    dp_lo &= DOT_SPLIT_MASK;

    ulong res;
    NMOD_RED(res, pow2_precomp * dp_hi + dp_lo, mod);
    return (uint)res;
}

#if PML_HAVE_AVX2

uint _nmod32_vec_dot_split_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);
    __m256i dp_lo0 = _mm256_setzero_si256();
    __m256i dp_hi0 = _mm256_setzero_si256();

    slong i = 0;

    for ( ; i+31 < len; i+=32)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 0));
        __m256i v1_1 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 8));
        __m256i v1_2 = _mm256_loadu_si256((const __m256i *) (vec1+i+16));
        __m256i v1_3 = _mm256_loadu_si256((const __m256i *) (vec1+i+24));
        __m256i v2_0 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 0));
        __m256i v2_1 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 8));
        __m256i v2_2 = _mm256_loadu_si256((const __m256i *) (vec2+i+16));
        __m256i v2_3 = _mm256_loadu_si256((const __m256i *) (vec2+i+24));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
        __m256i dp_lo1 = _mm256_mul_epu32(v1_1, v2_1);
        __m256i dp_lo2 = _mm256_mul_epu32(v1_2, v2_2);
        __m256i dp_lo3 = _mm256_mul_epu32(v1_3, v2_3);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v1_2 = _mm256_shuffle_epi32(v1_2, 0xB1);
        v1_3 = _mm256_shuffle_epi32(v1_3, 0xB1);
        v2_0 = _mm256_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm256_shuffle_epi32(v2_1, 0xB1);
        v2_2 = _mm256_shuffle_epi32(v2_2, 0xB1);
        v2_3 = _mm256_shuffle_epi32(v2_3, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_1, v2_1));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_2, v2_2));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_3, v2_3));

        // gather results in dp_lo0
        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo1);
        dp_lo2 = _mm256_add_epi64(dp_lo2, dp_lo3);
        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo2);

        // split
        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
    }

    // the following loop iterates < 4 times,
    // each iteration accumulates 2 terms in dp_lo0
    for ( ; i+7 < len; i+=8)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i));
        __m256i v2_0 = _mm256_loadu_si256((const __m256i *) (vec2+i));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v2_0 = _mm256_shuffle_epi32(v2_0, 0xB1);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
    }

    dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);

    ulong hsum_lo = _mm256_hsum(dp_lo0);
    ulong hsum_hi = _mm256_hsum(dp_hi0) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    // less than 8 terms remaining, can accumulate
    for (; i < len; i++)
        hsum_lo += (ulong)vec1[i] * vec2[i];

    hsum_hi += (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    // the requirement on "len <= DOT2_ACC8_MAX_LEN"
    // ensures pow2_precomp * hsum_hi + hsum_lo fits in 64 bits
    ulong res;
    NMOD_RED(res, pow2_precomp * hsum_hi + hsum_lo, mod);
    return (uint)res;
}

#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512

uint _nmod32_vec_dot_split_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64(DOT_SPLIT_MASK);
    __m512i dp_lo0 = _mm512_setzero_si512();
    __m512i dp_hi0 = _mm512_setzero_si512();

    slong i = 0;

    for ( ; i+63 < len; i+=64)
    {
        __m512i v1_0 = _mm512_loadu_si512((const __m512i *) (vec1+i+ 0));
        __m512i v1_1 = _mm512_loadu_si512((const __m512i *) (vec1+i+16));
        __m512i v1_2 = _mm512_loadu_si512((const __m512i *) (vec1+i+32));
        __m512i v1_3 = _mm512_loadu_si512((const __m512i *) (vec1+i+48));
        __m512i v2_0 = _mm512_loadu_si512((const __m512i *) (vec2+i+ 0));
        __m512i v2_1 = _mm512_loadu_si512((const __m512i *) (vec2+i+16));
        __m512i v2_2 = _mm512_loadu_si512((const __m512i *) (vec2+i+32));
        __m512i v2_3 = _mm512_loadu_si512((const __m512i *) (vec2+i+48));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        __m512i dp_lo1 = _mm512_mul_epu32(v1_1, v2_1);
        __m512i dp_lo2 = _mm512_mul_epu32(v1_2, v2_2);
        __m512i dp_lo3 = _mm512_mul_epu32(v1_3, v2_3);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm512_shuffle_epi32(v1_1, 0xB1);
        v1_2 = _mm512_shuffle_epi32(v1_2, 0xB1);
        v1_3 = _mm512_shuffle_epi32(v1_3, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm512_shuffle_epi32(v2_1, 0xB1);
        v2_2 = _mm512_shuffle_epi32(v2_2, 0xB1);
        v2_3 = _mm512_shuffle_epi32(v2_3, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm512_srli_epi64(v1_0, 32)

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_1, v2_1));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_2, v2_2));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_3, v2_3));

        // gather results in dp_lo0
        dp_lo0 = _mm512_add_epi64(dp_lo0, dp_lo1);
        dp_lo2 = _mm512_add_epi64(dp_lo2, dp_lo3);
        dp_lo0 = _mm512_add_epi64(dp_lo0, dp_lo2);

        // split
        dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo0
    for ( ; i+15 < len; i+=16)
    {
        __m512i v1_0 = _mm512_loadu_si512((const __m512i *) (vec1+i));
        __m512i v2_0 = _mm512_loadu_si512((const __m512i *) (vec2+i));

        // 1st term
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        // 2nd term
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
    }

    // finally, do a last iteration which may be "incomplete": use mask load
    // (at each position we accumulate 0 or 2 terms, so including the above
    // loop this is a total of <= 8 terms)
    if (i < len)
    {
        // mask == 0b0...01...1 with number of 1's == number of remaining terms
        __mmask16 mask = (1 << (len-i)) - 1;
        __m512i v1_0 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec1+i));
        __m512i v2_0 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2+i));

        // 1st term
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        // 2nd term
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
    }

    // split
    dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);

    // gather 8 terms in single ulong
    ulong hsum_lo = _mm512_hsum(dp_lo0);
    ulong hsum_hi = _mm512_hsum(dp_hi0) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    // the requirement "len <= DOT2_ACC8_MAX_LEN"
    // ensures pow2_precomp * hsum_hi + hsum_lo fits in 64 bits
    ulong res;
    NMOD_RED(res, pow2_precomp * hsum_hi + hsum_lo, mod);
    return (uint)res;
}


/*------------------*/
/* TODO IN PROGRESS */
/*------------------*/

uint _nmod32_vec_dot_ifma_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x((1L<<32) - 1);
    __m256i dp_lo0 = _mm256_setzero_si256();
    __m256i dp_hi0 = _mm256_setzero_si256();
    __m256i dp_lo1 = _mm256_setzero_si256();
    __m256i dp_hi1 = _mm256_setzero_si256();
    __m256i dp_lo2 = _mm256_setzero_si256();
    __m256i dp_hi2 = _mm256_setzero_si256();
    __m256i dp_lo3 = _mm256_setzero_si256();
    __m256i dp_hi3 = _mm256_setzero_si256();

    slong i = 0;

    for ( ; i+31 < len; i += 32)
    {
        __m256i vec1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 0));
        __m256i vec2_0 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 0));
        __m256i vec1_1 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 8));
        __m256i vec2_1 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 8));
        __m256i vec1_2 = _mm256_loadu_si256((const __m256i *) (vec1+i+16));
        __m256i vec2_2 = _mm256_loadu_si256((const __m256i *) (vec2+i+16));
        __m256i vec1_3 = _mm256_loadu_si256((const __m256i *) (vec1+i+24));
        __m256i vec2_3 = _mm256_loadu_si256((const __m256i *) (vec2+i+24));

        // handle low 32 bit word of each 64 bit word
        __m256i v1_0 = _mm256_and_si256(vec1_0, low_bits);
        __m256i v1_1 = _mm256_and_si256(vec1_1, low_bits);
        __m256i v2_0 = _mm256_and_si256(vec2_0, low_bits);
        __m256i v2_1 = _mm256_and_si256(vec2_1, low_bits);
        __m256i v1_2 = _mm256_and_si256(vec1_2, low_bits);
        __m256i v1_3 = _mm256_and_si256(vec1_3, low_bits);
        __m256i v2_2 = _mm256_and_si256(vec2_2, low_bits);
        __m256i v2_3 = _mm256_and_si256(vec2_3, low_bits);
        dp_lo0 = _mm256_madd52lo_epu64(dp_lo0, v1_0, v2_0);
        dp_hi0 = _mm256_madd52hi_epu64(dp_hi0, v1_0, v2_0);
        dp_lo1 = _mm256_madd52lo_epu64(dp_lo1, v1_1, v2_1);
        dp_hi1 = _mm256_madd52hi_epu64(dp_hi1, v1_1, v2_1);
        dp_lo2 = _mm256_madd52lo_epu64(dp_lo2, v1_2, v2_2);
        dp_hi2 = _mm256_madd52hi_epu64(dp_hi2, v1_2, v2_2);
        dp_lo3 = _mm256_madd52lo_epu64(dp_lo3, v1_3, v2_3);
        dp_hi3 = _mm256_madd52hi_epu64(dp_hi3, v1_3, v2_3);

        // handle high 32 bit word of each 64 bit word
        v1_0 = _mm256_srli_epi64(vec1_0, 32);
        v2_0 = _mm256_srli_epi64(vec2_0, 32);
        v1_1 = _mm256_srli_epi64(vec1_1, 32);
        v2_1 = _mm256_srli_epi64(vec2_1, 32);
        v1_2 = _mm256_srli_epi64(vec1_2, 32);
        v2_2 = _mm256_srli_epi64(vec2_2, 32);
        v1_3 = _mm256_srli_epi64(vec1_3, 32);
        v2_3 = _mm256_srli_epi64(vec2_3, 32);
        dp_lo0 = _mm256_madd52lo_epu64(dp_lo0, v1_0, v2_0);
        dp_hi0 = _mm256_madd52hi_epu64(dp_hi0, v1_0, v2_0);
        dp_lo1 = _mm256_madd52lo_epu64(dp_lo1, v1_1, v2_1);
        dp_hi1 = _mm256_madd52hi_epu64(dp_hi1, v1_1, v2_1);
        dp_lo2 = _mm256_madd52lo_epu64(dp_lo2, v1_2, v2_2);
        dp_hi2 = _mm256_madd52hi_epu64(dp_hi2, v1_2, v2_2);
        dp_lo3 = _mm256_madd52lo_epu64(dp_lo3, v1_3, v2_3);
        dp_hi3 = _mm256_madd52hi_epu64(dp_hi3, v1_3, v2_3);
    }

    ulong hsum_lo = _mm256_hsum(dp_lo0) + _mm256_hsum(dp_lo1) + _mm256_hsum(dp_lo2) + _mm256_hsum(dp_lo3);
    ulong hsum_hi = _mm256_hsum(dp_hi0) + _mm256_hsum(dp_hi1) + _mm256_hsum(dp_hi2) + _mm256_hsum(dp_hi3);

    // TODO handle remaining terms
    //for (; i < len; i++)
    //{
    //    hsum_lo += (ulong)vec1[i] * (ulong)vec2[i];
    //    hsum_hi += (hsum_lo >> DOT_SPLIT_BITS);
    //    hsum_lo &= DOT_SPLIT_MASK;
    //}

    NMOD_RED(hsum_lo, hsum_lo, mod);
    NMOD_RED(hsum_hi, hsum_hi, mod);
    ulong res = nmod_add(hsum_lo, nmod_mul(pow2_precomp, hsum_hi, mod), mod);
    return (uint)res;
}

uint _nmod32_vec_dot_ifma_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64((1L<<32) - 1);
    __m512i dp_lo0 = _mm512_setzero_si512();
    __m512i dp_hi0 = _mm512_setzero_si512();
    __m512i dp_lo1 = _mm512_setzero_si512();
    __m512i dp_hi1 = _mm512_setzero_si512();

    slong i = 0;

    for ( ; i+31 < len; i += 32)
    {
        __m512i vec1_0 = _mm512_loadu_si512((const __m512i *) (vec1+i+ 0));
        __m512i vec2_0 = _mm512_loadu_si512((const __m512i *) (vec2+i+ 0));
        __m512i vec1_1 = _mm512_loadu_si512((const __m512i *) (vec1+i+16));
        __m512i vec2_1 = _mm512_loadu_si512((const __m512i *) (vec2+i+16));

        // handle low 32 bit word of each 64 bit word
        __m512i v1_0 = _mm512_and_si512(vec1_0, low_bits);
        __m512i v1_1 = _mm512_and_si512(vec1_1, low_bits);
        __m512i v2_0 = _mm512_and_si512(vec2_0, low_bits);
        __m512i v2_1 = _mm512_and_si512(vec2_1, low_bits);
        dp_lo0 = _mm512_madd52lo_epu64(dp_lo0, v1_0, v2_0);
        dp_hi0 = _mm512_madd52hi_epu64(dp_hi0, v1_0, v2_0);
        dp_lo1 = _mm512_madd52lo_epu64(dp_lo1, v1_1, v2_1);
        dp_hi1 = _mm512_madd52hi_epu64(dp_hi1, v1_1, v2_1);

        // handle high 32 bit word of each 64 bit word
        v1_0 = _mm512_srli_epi64(vec1_0, 32);
        v2_0 = _mm512_srli_epi64(vec2_0, 32);
        v1_1 = _mm512_srli_epi64(vec1_1, 32);
        v2_1 = _mm512_srli_epi64(vec2_1, 32);
        dp_lo0 = _mm512_madd52lo_epu64(dp_lo0, v1_0, v2_0);
        dp_hi0 = _mm512_madd52hi_epu64(dp_hi0, v1_0, v2_0);
        dp_lo1 = _mm512_madd52lo_epu64(dp_lo1, v1_1, v2_1);
        dp_hi1 = _mm512_madd52hi_epu64(dp_hi1, v1_1, v2_1);
    }

    ulong hsum_lo = _mm512_hsum(dp_lo0) + _mm512_hsum(dp_lo1);
    ulong hsum_hi = _mm512_hsum(dp_hi0) + _mm512_hsum(dp_hi1);

    // TODO handle remaining terms
    //for (; i < len; i++)
    //{
    //    hsum_lo += (ulong)vec1[i] * (ulong)vec2[i];
    //    hsum_hi += (hsum_lo >> DOT_SPLIT_BITS);
    //    hsum_lo &= DOT_SPLIT_MASK;
    //}

    NMOD_RED(hsum_lo, hsum_lo, mod);
    NMOD_RED(hsum_hi, hsum_hi, mod);
    ulong res = nmod_add(hsum_lo, nmod_mul(pow2_precomp, hsum_hi, mod), mod);
    return (uint)res;
}

#endif  /* PML_HAVE_AVX512 */
