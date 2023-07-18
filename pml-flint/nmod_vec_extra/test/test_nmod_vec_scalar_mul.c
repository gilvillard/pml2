#include <assert.h>
#include <flint/flint.h>

#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/* multiplies a vector by a scalar                            */
/*------------------------------------------------------------*/
void check()
{
    slong i;
    flint_rand_t state;
    mp_limb_t nn;
    nmod_t mod;
    
    flint_randinit(state);
    nn = 65537;
    nmod_init(&mod, nn);

    for (i = 1; i < 100; i+=1)
    {
        slong j;
        mp_limb_t r;
        mp_ptr val, val2;
        
        val = _nmod_vec_init(i);
        val2 = _nmod_vec_init(i);
        _nmod_vec_randtest(val, state, i, mod);
        
        for (j = 0; j < i; j++)
        {
            val2[j] = val[j];
        }

        r = nn/2;
        _nmod_vec_scalar_mul(val, val, i, r, mod);

        for (j = 0; j < i; j++)
        {
            assert (val[j] == nmod_mul(val2[j], r, mod));
        }
        
        _nmod_vec_clear(val);
        _nmod_vec_clear(val2);
    }
    flint_randclear(state);
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
