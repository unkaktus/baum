 
 This BAM project implements the "matter" part of Einstein equations.
 Here "matter" can be general, and refers to a RHS of an evolution
 system in flux-balance form (conservation law) on 3+1 dynamical
 spacetimes, i.e.

 (1)    d q /dt = d F /dx + s
 
 where

         q conservative vars
         F fluxes
         s sources

 and additional primitives variables, w, as well as algebraic
 relations E(w)=0 (the equation of state) are usually necessary to
 close the system.

 The system (1) is solved with an high-resolution-shock-capturing
 (HRSC) scheme which relies on 
 - method-of-line time integration (shared with the spacetime)
 - un-split numerical fluxes constructed as successive 1d operations 
 Within this framework different schemes for the numerical fluxes can
 be implemented. 

 The project interfaces with the spacetime evolution via ADM
 variables, that couple matter with spacetime.

 The project usually needs another one which implements the routines
 for 1d fluxes, conservative<=>primitives, source terms, etc, specific
 for given conservation law system. (see e.g. grhd/ )
     
 For more see documentation ( doc/ ) 

 sbernuz, mth 03/2012

   
