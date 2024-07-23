/* mth 06/2010 */

this is a test-case for testing WENO and Lagrane interpolation

the test function is
    f = sqrt(0.25-x)    | 0.   < x < 0.25
    f = sin(2.*Pi*x)    | 0.25 < x < 0.5
    f = 1-sin(2.*Pi*x)  | 0.5  < x < 0.75
    f = 1               | 0.75 < x < 1.0
with 100 points

this is interpolated to a finer grid (half the gridspace) 
like bam uses meshrefinement

x   x   x   x       / coarser grid
 o o o o o o        / finer grid

several testcases:
    Lagrange  2,4,6,8
    WENO      4,6

