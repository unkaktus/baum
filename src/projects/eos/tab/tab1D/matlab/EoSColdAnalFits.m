%EoSColdAnalFits analytic fits of neutron-star cold EoS SLy,FPS and APR  
%
%   This is a MATLAB interface for the routines in 
%                   HPFits.c
%                   ShibFits.c
%   which implement analytical representation of neutron-star cold
%   EoS SLy, FPS and APR following:
%
%   [1] Haensel, P.; Potekhin, A. Y. 
%       Astronomy and Astrophysics 428, p.191-197 (2004) 
%       arXiv:astro-ph/0408324 
%   [2] Shibata, M. et al., 
%       Phys. Rev. D71, 084021 (2005) 
%       arXiv:gr-qc/0503119
%       Phys. Rev. D73, 064027 (2006)
%       arXiv:astro-ph/0603145 
%
%   Options for the fitting formulas P(rho):
%
%   'hp'   : epsl(rho) Eq.(15) of [1] + P(ene) Eq.(14) of [1] 
%
%   'hps'  : epsl(rho) Eq.(15) of [1] + P from thermodynamics [2]
%
%   'shib' : epsl(rho) Eq.(25) of [2] + P from thermodynamics [2]
%
%   Usage:
%
%   [p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits(eos,fit,rho);
%   [p,DpDrho,DpDepsl,cs2,ene,epsl,DepslDrho,D2epslD2rho]=EoSColdAnalFits(eos,fit,rho);
%
%   eos = eos name <string> {'sly' , 'fps' , 'apr' }
%   fit = fit name <string> {'hp' , 'hps' , 'shib' } 
%   rho = rest-mass density <vector>
%   
%   Units: Msun=c=G=1 
%
%   Author: SBernuz 01/11
%

