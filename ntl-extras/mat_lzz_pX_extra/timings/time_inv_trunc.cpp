#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

PML_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    const double thres = 0.5;
    double t_middle, t_geometric;
    long nb;
    Mat<zz_pX> a, x;
    Mat<zz_p> a0;

    do
    {
        random(a, sz, sz, deg);
        GetCoeff(a0, a, 0);
    }
    while (determinant(a0) == 0);

    t_middle = GetWallTime();
    nb = 0;
    do
    {
        newton_inv_trunc_middle_product(x, a, deg);
        nb++;
    }
    while ((GetWallTime()-t_middle) <= thres);
    t_middle = (GetWallTime()-t_middle) / nb;

    t_geometric = GetWallTime();
    nb = 0;
    do
    {
        newton_inv_trunc_geometric(x, a, deg);
        nb++;
    }
    while ((GetWallTime()-t_geometric) <= thres);
    t_geometric = (GetWallTime()-t_geometric) / nb;

    cout << sz << " " << deg << " " << t_middle << " " << t_geometric << " ";

    if (is_FFT_prime())
    {
        double t_FFT;
        t_FFT = GetWallTime();
        nb = 0;
        do
        {
            newton_inv_trunc_FFT(x, a, deg);
            nb++;
        }
        while ((GetWallTime()-t_FFT) <= thres);
        t_FFT = (GetWallTime()-t_FFT) / nb;
        cout << t_FFT;
    }
    cout << endl;

}

/*------------------------------------------------------------*/
/* for a given prime, checks some (size, degree)              */
/*------------------------------------------------------------*/
void all_checks()
{

    VecLong szs =
    {
        1, 2, 3, 5, 10, 20, 30, 50, 100
    };

    VecLong degs =
    {
        20, 50, 75, 99, 150, 200, 500, 1000
    };

    cout << "p=" << zz_p::modulus() << "\nFFT=" << is_FFT_prime() << endl;
    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check(szs[si], degs[di]);

}



/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);
    all_checks();
    //zz_p::UserFFTInit(786433);
    //all_checks();
    //zz_p::init(288230376151711813);
    //all_checks();
    //zz_p::init(786433);
    //all_checks();
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup();
    if (argc == 1)
        check();
    else if (argc == 3) // sz, deg
    {
        zz_p::FFTInit(0);
        one_check(atoi(argv[1]), atoi(argv[2]));
    }
    else
    {
        std::cout << "Usage: ./" << argv[0] << " OR ./" << argv[0] << "sz deg" << std::endl;
    }
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
