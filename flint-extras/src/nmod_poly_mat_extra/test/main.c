/*
    This file is part of PML

    Using FLINT functionalities 
*/

/* Include functions *********************************************************/

#include "t-dixon.c"
#include "t-description.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_poly_mat_dixon),
    TEST_FUNCTION(nmod_poly_mat_dixon);
};

/* main function *************************************************************/

TEST_MAIN(tests)
