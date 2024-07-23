
Nina and Wolfgang made the suggestion that the Mathematica file should be
split in several pieces. Here is a starting point.

There is a script "TensorEquationsToC.m", which is the only script
that the user loads, but which in turn loads other .m files in the
math directory.

It would be neat to have not only a testsuite for numerical outputs,
but also a test suite which tests whether the C code is still
identical!  I.e. whether some change in a shift derivative introduced
an unexpected bug somewhere else.

