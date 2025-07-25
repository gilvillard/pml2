#include <algorithm> // for minmax
#include "mat_lzz_pX_multiply.h"

PML_START_IMPL

/*------------------------------------------------------------*/
/* in-place reduction modulo the current prime                */
/*------------------------------------------------------------*/
void reduce_mod_p(Mat<zz_pX> & a)
{
    const long p = zz_p::modulus();
    sp_reduce_struct red_struct = zz_pInfo->red_struct;
    for (long i = 0; i < a.NumRows(); ++i) 
        for (long j = 0; j < a.NumCols(); ++j)
        {
            for (long k = 0; k <= deg(a[i][j]); ++k)
                a[i][j][k].LoopHole() = rem(a[i][j][k]._zz_p__rep, p, red_struct);
            a[i][j].normalize();
        }
}

/*------------------------------------------------------------*/
/* copy and reduce modulo the current prime                   */
/*------------------------------------------------------------*/
void reduce_mod_p(Mat<zz_pX> & ap, const Mat<zz_pX> & a)
{
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = zz_p::modulus();
    sp_reduce_struct red_struct = zz_pInfo->red_struct;
    ap.SetDims(m,n);
    for (long i = 0; i < m; ++i) 
        for (long j = 0; j < n; ++j)
        {
            const long d = deg(a[i][j]);
            ap[i][j].SetLength(d+1);
            for (long k = 0; k <= d; ++k)
                ap[i][j][k].LoopHole() = rem(a[i][j][k]._zz_p__rep, p, red_struct);
            ap[i][j].normalize();
        }
}

/*------------------------------------------------------------*/
/* matrix CRT modulo 2 primes, result reduced modulo p        */
/* c cannot alias c0 or c1; c does not have to be zero        */
/*------------------------------------------------------------*/
static void reconstruct_2CRT(Mat<zz_pX> & c, const Mat<zz_pX> & c0, long p0, const Mat<zz_pX> & c1, long p1)
{
    if (p0 > p1) // ensures that p0 <= p1
    {
        reconstruct_2CRT(c, c1, p1, c0, p0);
        return;
    }

    const long p0_inv = InvMod(p0, p1);
    mulmod_precon_t p0_inv_prec = PrepMulModPrecon(p0_inv, p1, PrepMulMod(p1)); 

    const long p = zz_p::modulus();
    sp_reduce_struct red_struct = zz_pInfo->red_struct;
    const long p0_red = rem(p0, p, red_struct);
    mulmod_precon_t p0_prec = PrepMulModPrecon(p0_red, p, PrepMulMod(p)); 

    const long r = c0.NumRows();
    const long s = c0.NumCols();
    c.SetDims(r, s);

    for (long i = 0; i < r; ++i)
    {
        for (long j = 0; j < s; ++j)
        {
            const std::pair<long,long> degs = std::minmax(deg(c0[i][j]), deg(c1[i][j]));
            c[i][j].SetLength(degs.second+1);
            for (long k = 0; k <= degs.first; ++k)
            {
                long m0 = c0[i][j][k]._zz_p__rep;
                long alpha = SubMod(c1[i][j][k]._zz_p__rep, m0, p1);
                alpha = MulModPrecon(alpha, p0_inv, p1, p0_inv_prec);
                alpha = rem(alpha, p, red_struct);
                m0 = rem(m0, p, red_struct);
                c[i][j][k] = AddMod(m0, MulModPrecon(alpha, p0_red, p, p0_prec), p);
            }
            for (long k = degs.first+1; k <= degs.second; ++k)
            {
                long m0 = coeff(c0[i][j], k)._zz_p__rep;
                long alpha = SubMod(coeff(c1[i][j], k)._zz_p__rep, m0, p1);
                alpha = MulModPrecon(alpha, p0_inv, p1, p0_inv_prec);
                alpha = rem(alpha, p, red_struct);
                m0 = rem(m0, p, red_struct);
                c[i][j][k] = AddMod(m0, MulModPrecon(alpha, p0_red, p, p0_prec), p);
            }
            c[i][j].normalize();
        }
    }
}

/*------------------------------------------------------------*/
/* matrix CRT modulo 3 primes, result reduced modulo p        */
/* c cannot alias c0, c1 or c2; c does not have to be zero    */
/*------------------------------------------------------------*/
static void reconstruct_3CRT(Mat<zz_pX> & c, const Mat<zz_pX> & c0, long p0, const Mat<zz_pX> & c1, long p1, const Mat<zz_pX> & c2, long p2)
{

    // ensure that p0 <= p1 <= p2
    if (p0 > p1) // ensures that p0 <= p1
    {
        reconstruct_3CRT(c, c1, p1, c0, p0, c2, p2);
        return;
    }

    if (p1 > p2) // ensures that p1 <= p2
    {
        reconstruct_3CRT(c, c0, p0, c2, p2, c1, p1);
        return;
    }	

    long p = zz_p::modulus();
    sp_reduce_struct red_struct = zz_pInfo->red_struct;
    // sp_reduce_struct red_struct = zz_p::red_struct();

    // p0 mod p
    long p0_p = rem(p0, p, red_struct);
    mulmod_precon_t p0_p_prec = PrepMulModPrecon(p0_p, p, PrepMulMod(p)); 
    // p1 mod p
    long p1_p = rem(p1, p, red_struct);
    // p0p1 mod p
    long p0p1_p = MulMod(p0_p, p1_p, p);
    mulmod_precon_t p0p1_p_prec = PrepMulModPrecon(p0p1_p, p, PrepMulMod(p)); 

    // 1/p0 mod p1
    long p0_inv_p1  = InvMod(p0, p1);
    mulmod_precon_t p0_inv_p1_prec = PrepMulModPrecon(p0_inv_p1, p1, PrepMulMod(p1)); 

    // p0 mod p2 -- already reduced
    mulmod_precon_t p0_p2_prec = PrepMulModPrecon(p0, p2, PrepMulMod(p2)); 
    // p0p1 and 1/p0p1 mod p2;
    long p0p1_p2 = MulMod(p0, p1, p2);
    long p0p1_inv_p2 = InvMod(p0p1_p2, p2);
    mulmod_precon_t p0p1_inv_p2_prec = PrepMulModPrecon(p0p1_inv_p2, p2, PrepMulMod(p2)); 

    long r = c0.NumRows();
    long s = c0.NumCols();
    c.SetDims(r, s);

    for (long i = 0; i < r; ++i)
    {
        for (long j = 0; j < s; ++j)
        {
            const std::pair<long,long> degs = std::minmax({deg(c0[i][j]), deg(c1[i][j]), deg(c2[i][j])});
            c[i][j].SetLength(degs.second+1);
            for (long k = 0; k <= degs.first; ++k)
            {
                // find coefficients alpha0, alpha1, alpha2 s.t. c = alpha0 + p0 alpha1 + p0 p1 alpha2
                long alpha0 = c0[i][j][k]._zz_p__rep;
                long m1 = c1[i][j][k]._zz_p__rep;
                long m2 = c2[i][j][k]._zz_p__rep;

                long alpha1 = SubMod(m1, alpha0, p1);
                alpha1 = MulModPrecon(alpha1, p0_inv_p1, p1, p0_inv_p1_prec);

                // p0 < p1 < p2 so alpha0, alpha1 reduced mod p2
                long alpha1_p0_p2 = MulModPrecon(alpha1, p0, p2, p0_p2_prec);
                alpha1_p0_p2 = AddMod(alpha0, alpha1_p0_p2, p2);
                long alpha2 = SubMod(m2, alpha1_p0_p2, p2);
                alpha2 = MulModPrecon(alpha2, p0p1_inv_p2, p2, p0p1_inv_p2_prec);

                // reduce everything mod p
                alpha0 = rem(alpha0, p, red_struct);
                alpha1 = rem(alpha1, p, red_struct);
                alpha2 = rem(alpha2, p, red_struct);

                c[i][j][k] = AddMod(alpha0, 
                                    AddMod(MulModPrecon(alpha1, p0_p, p, p0_p_prec), 
                                           MulModPrecon(alpha2, p0p1_p, p, p0p1_p_prec), p), p);
            }
            for (long k = degs.first+1; k <= degs.second; ++k)
            {
                // find coefficients alpha0, alpha1, alpha2 s.t. c = alpha0 + p0 alpha1 + p0 p1 alpha2
                long alpha0 = coeff(c0[i][j], k)._zz_p__rep;
                long m1 = coeff(c1[i][j], k)._zz_p__rep;
                long m2 = coeff(c2[i][j], k)._zz_p__rep;

                long alpha1 = SubMod(m1, alpha0, p1);
                alpha1 = MulModPrecon(alpha1, p0_inv_p1, p1, p0_inv_p1_prec);

                // p0 < p1 < p2 so alpha0, alpha1 reduced mod p2
                long alpha1_p0_p2 = MulModPrecon(alpha1, p0, p2, p0_p2_prec);
                alpha1_p0_p2 = AddMod(alpha0, alpha1_p0_p2, p2);
                long alpha2 = SubMod(m2, alpha1_p0_p2, p2);
                alpha2 = MulModPrecon(alpha2, p0p1_inv_p2, p2, p0p1_inv_p2_prec);

                // reduce everything mod p
                alpha0 = rem(alpha0, p, red_struct);
                alpha1 = rem(alpha1, p, red_struct);
                alpha2 = rem(alpha2, p, red_struct);

                c[i][j][k] = AddMod(alpha0, 
                                    AddMod(MulModPrecon(alpha1, p0_p, p, p0_p_prec), 
                                           MulModPrecon(alpha2, p0p1_p, p, p0p1_p_prec), p), p);
            }
            c[i][j].normalize();
        }
    }
}

/*------------------------------------------------------------*/
/* FFT multiplication done modulo FFTprime[idx]               */
/* c can alias a or b; c does not have to be zero             */
/*------------------------------------------------------------*/
static void multiply_modulo_FFT_prime(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long idx)
{
    long p = zz_p::modulus();
    zz_pPush push;
    zz_p::FFTInit(idx);
    long fft_p = zz_p::modulus();

    if (fft_p < p) // entries may not be reduced mod fft_p
    {
        Mat<zz_pX> ap, bp;
        reduce_mod_p(ap, a);
        reduce_mod_p(bp, b);
        multiply_evaluate_FFT(c, ap, bp);
    }
    else
        multiply_evaluate_FFT(c, a, b);
}


/*------------------------------------------------------------*/
/* middle product done modulo FFTprime[idx]                   */
/* b can alias a or c; b does not have to be zero             */
/*------------------------------------------------------------*/
static void middle_product_modulo_FFT_prime(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB, long idx)
{
    long p = zz_p::modulus();
    zz_pPush push;
    zz_p::FFTInit(idx);
    long fft_p = zz_p::modulus();

    if (fft_p < p) // entries may not be reduced mod fft_p
    {
        Mat<zz_pX> ap, cp;
        reduce_mod_p(ap, a);
        reduce_mod_p(cp, c);
        middle_product_evaluate_FFT(b, ap, cp, dA, dB);
    }
    else
        middle_product_evaluate_FFT(b, a, c, dA, dB);
}

/*------------------------------------------------------------*/
/* constructor of lzz_p_3_primes                              */
/*------------------------------------------------------------*/
lzz_pX_3_primes::lzz_pX_3_primes(long ncols, long dA, long dB)
{
    long p = zz_p::modulus();
    ZZ max_coeff = ZZ(ncols) * ZZ(min(dA + 1, dB + 1)) * ZZ(p - 1) * ZZ(p - 1);
    {
        zz_pPush push;
        zz_p::FFTInit(0);
        fft_p0 = zz_p::modulus();
        zz_p::FFTInit(1);
        fft_p1 = zz_p::modulus();
        zz_p::FFTInit(2);
        fft_p2 = zz_p::modulus();
    }

    nb_primes = 0;
    if (ZZ(fft_p0) > max_coeff)
        nb_primes = 1;
    else if (ZZ(fft_p0)*ZZ(fft_p1) > 2*max_coeff)
        nb_primes = 2;
    else if (ZZ(fft_p0)*ZZ(fft_p1)*ZZ(fft_p2) > 2*max_coeff)   // x2 so that normalized remainders are in [0..p0p1p2/2]
        nb_primes = 3;
    else
        LogicError("size too large for 3 primes FFT");
}

/*------------------------------------------------------------*/
/* reconstructs c from its images                             */
/*------------------------------------------------------------*/
void lzz_pX_3_primes::reconstruct(Mat<zz_pX>& c, const Vec<Mat<zz_pX>>& cs) const
{
    switch(nb_primes)
    {
    case 1:
        c = cs[0];
        reduce_mod_p(c);
        break;
    case 2:
        reconstruct_2CRT(c, cs[0], fft_p0, cs[1], fft_p1);
        break;
    case 3:
        reconstruct_3CRT(c, cs[0], fft_p0, cs[1], fft_p1, cs[2], fft_p2);
        break;
    default:
        LogicError("impossible branch for 3 primes");
    }
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* chooses either 1, 2 or 3 FFT primes                        */
/* c can alias a or b; c does not have to be zero             */
/*------------------------------------------------------------*/
void multiply_3_primes(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    lzz_pX_3_primes primes(a.NumCols(), deg(a), deg(b));
    const long nb = primes.nb();

    if (nb == 1)
    {
        multiply_modulo_FFT_prime(c, a, b, 0); // no need to copy stuff around
        reduce_mod_p(c);
        return;
    }
    else
    {
        Vec<Mat<zz_pX>> cs;
        cs.SetLength(nb);
        for (long i = 0; i < nb; ++i)
            multiply_modulo_FFT_prime(cs[i], a, b, i);
        primes.reconstruct(c, cs);
    }
}


/*------------------------------------------------------------*/
/* middle product                                             */
/* chooses either 1, 2 or 3 FFT primes                        */
/* b can alias a or c; b does not have to be zero             */
/*------------------------------------------------------------*/
void middle_product_3_primes(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    lzz_pX_3_primes primes(a.NumCols(), dA, dB);
    const long nb = primes.nb();

    if (nb == 1)
    {
        middle_product_modulo_FFT_prime(b, a, c, dA, dB, 0);
        reduce_mod_p(b);
        return;
    }
    else
    {
        Vec<Mat<zz_pX>> bs;
        bs.SetLength(nb);
        for (long i = 0; i < nb; ++i)
            middle_product_modulo_FFT_prime(bs[i], a, c, dA, dB, i);
        primes.reconstruct(b, bs);
    }
}

PML_END_IMPL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
