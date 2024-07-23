function [ rho, epsl, pre ] = writeTab1DPoly( fname, n, K, gamma )
%WRITETAB1DPOLY writes 1D tables with polytropic EoS
%
%                         P(rho) = K rho^gamma
%
%         [ rho, epsl, pre ] =  writeTab1DPoly( fname, n, K, gamma )
%
%         The script computes the EoS in dimensionless units 
%                     c = G = Msun = 1
%         using n points logaritmically spaced: rho in [1e-16 5e-1]
%
%         It returns rho (rest-mass density), epsl (specific
%         energy) and p (pressure) in dimensionless units and a file
%         'fname' with the line  
%           n  <-- Number of lines
%         and 4 cols
%           index  n_B [fm^{-3}]  e [g/cm^3]   p [dyn/cm^2]
%
%         n_B = baryon number density, 
%         e   = total energy density,
%         p   = pressure.
%
%         plus comments lines starting with # 
         
% units 
utsMdens_cgs = 6.1764e+17; % g cm^-3
utsPress_cgs = 5.5511e+38; % dyn cm^-2
m0           = 1.66e-24;   % g
utsN_cgs     = 1e+39;      % cm^-3

% dimensionless EoS
rho  = logspace(-16,-1,n);
pre  = K*rho.^gamma;
epsl = pre./((gamma-1)*rho);
ene  = rho.*(epsl+1);

% change units
tab_nB  = rho *utsMdens_cgs/(utsN_cgs*m0);
tab_ene = ene *utsMdens_cgs;
tab_pre = pre *utsPress_cgs;

% write data
fprintf(1,'file : %s\n', fname)

fid = fopen( fname, 'wt' );
fprintf(fid,'#  Date: %s\n',date);
fprintf(fid,'#  Created with: writeTab1DPoly.m (Author: S.Bernuzzi)\n');
fprintf(fid,'#  Polytropic EoS: K=%g Gamma=%g (c=G=Msunc=1)\n',K,gamma);
fprintf(fid,'%d  <-- Number of lines\n',n);
fprintf(fid,'#\n');
fprintf(fid,'#            n_B [fm^{-3}]       rho [g/cm^3]        p [dyn/cm^2]\n');
fprintf(fid,'#\n');
for k=1:n
  fprintf(fid,'%d %2.16e %2.16e %2.16e\n',...
          k,tab_nB(k),tab_ene(k),tab_pre(k));
end
fclose(fid);


