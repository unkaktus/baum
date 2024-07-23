This directory contains 1D (Cold) EoS tables
Taken mostly from Lorene lib
They can be directly read into bam passing the file name

Table format is the following:

	#.... comments
	#.... 
	#
	71  <-- Number of lines
	#
	#        n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]
	#
	    1    0.979240E-14    0.131971E+02    0.302909E+15
	    ...	 ...		 ..		 ...

Remark:

    n_B = baryon number density,
    e   = total energy density,
    p   = pressure

The subdir matlab/ contains mex/matlab routines to produce analytic
fits of SLy, APR and FPS tables.
