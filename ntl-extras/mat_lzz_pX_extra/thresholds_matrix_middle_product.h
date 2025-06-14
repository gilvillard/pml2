#ifndef __THRESHOLDS_MATRIX_MIDDLE_PRODUCT__H
#define __THRESHOLDS_MATRIX_MIDDLE_PRODUCT__H

#include "lzz_p_extra.h"
#include "thresholds_mp_naive_evaluate.h"
#include "util.h"

PML_OPEN_NNS
NTL_USE_NNS

/*------------------------------------------------------------*/
/* max degree for naive                                       */
/*------------------------------------------------------------*/
inline long max_degree_mp_naive(long sz)
{
    long t = type_of_prime();
    long i;
    for (i = 0; i < MATRIX_MP_THRESHOLDS_LEN - 1; )
        if (sz > MATRIX_MP_THRESHOLDS_SIZES[i])
            i++;
        else 
            break;

    if (t == TYPE_FFT_PRIME)
        return MATRIX_MP_NAIVE_THRESHOLDS_FFT[i];
    else 
        if (t == TYPE_SMALL_PRIME)
            return MATRIX_MP_NAIVE_THRESHOLDS_SMALL[i];
        else
            return MATRIX_MP_NAIVE_THRESHOLDS_LARGE[i];
}

PML_CLOSE_NNS

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
