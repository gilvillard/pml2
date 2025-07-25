#include <assert.h>
#include <NTL/vec_lzz_p.h>

#include "util.h"
#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

PML_CLIENT

long next(long i, long step, long q)
{
    return (i * step) % q;
}

/*------------------------------------------------------------*/
/* solves some systems                                        */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    Vec<long> rows;
    Vec<long> cols;
  
    long idx = 1;
    long step = 19;
    long per = 101;

    for (long rep = 0; rep < 100; rep++)
        for (long i = 1; i < 5; i++)
            for (long j = 1; j < 5; j++)
            {
                rows.SetLength(i);
                cols.SetLength(j);
                
                for (long k = 0; k < i; k++)
                {
                    idx = next(idx, step, per);
                    rows[k] = idx;
                }
                
                for (long k = 0; k < j; k++)
                {
                    idx = next(idx, step, per);
                    cols[k] = idx;
            }
                
                Vec<Vec<toeplitz_lzz_p>> mat;
                mat.SetLength(i);
                for (long k = 0; k < i; k++)
                {
                    Vec<toeplitz_lzz_p> row;
                    row.SetLength(j);
                    for (long ell = 0; ell < j; ell++)
                        row[ell] = toeplitz_lzz_p(random_vec_zz_p(rows[k] + cols[ell] - 1), rows[k], cols[ell]);
                    mat[k] = row;
                }
                mosaic_toeplitz_lzz_p MH = mosaic_toeplitz_lzz_p(mat);
                long n = MH.NumCols();
                Vec<zz_p> u, v, w;
                u = random_vec_zz_p(n);
                v = MH.mul_right(u);
                long r = MH.solve(w, v);
                assert (r == 1);
                assert (MH.mul_right(w) == v);
            }
           
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(0);
    check(786433);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
