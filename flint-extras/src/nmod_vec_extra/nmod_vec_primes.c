#include <flint/ulong_extras.h>  // for n_nextprime

#include "nmod_vec_extra.h"

/*--------------------------------------------------------------*/
/* vector of n consecutive primes of exactly s bits             */
/*--------------------------------------------------------------*/
void nmod_vec_primes(nn_ptr v, slong n, flint_bitcnt_t s)
{
    slong i;
    v[0] = n_nextprime(UWORD(1) << (s-1), 0);
    for (i = 1; i < n; i++)
    {
        v[i] = n_nextprime(v[i-1], 0);
    }

    if (FLINT_BIT_COUNT(v[n-1]) > s)
    {
	flint_throw(FLINT_ERROR, "too few primes of given bitlength\n");
    }
}
