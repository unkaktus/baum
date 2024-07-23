%MATLAB SCRIPT to compile mex-files

fprintf(1,'===> compiling EoSColdAnalFitsMex.c ...')
mex EoSColdAnalFitsMex.c -v CFLAGS="\$CFLAGS -Wall" -o EoSColdAnalFits
fprintf(1,' done\n')

