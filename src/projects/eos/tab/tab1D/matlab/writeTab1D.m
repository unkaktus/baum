function varargout = writeTab1D( fname, rho, epsl, pre )
%WRITETAB1D writes 1D EoS tables given input quanities
%
%         Usage:
%
%         writeTab1D(fname,rho,epsl,pre);
%         [tab_rho,tab_epsl,tab_pre] =  writeTab1D(fname,rho,epsl,pre);
%
%         The script requires rest-mass energy denisty, specific
%         energy and pressure in dimensionless units 
%                     c = G = Msun = 1
%         and output them in an ASCII table constructed as follows
%
%              # ... comments ...
%              #
%              n  <-- Number of lines
%              # ... comments ...
%              # ... comments ...
%                  index  n_B [fm^{-3}]  e [g/cm^3]   p [dyn/cm^2]
%                   ...     ...             ...          ...
%                  
%        with in the
%
%         n_B = baryon number density [fm^{-3}]
%         e   = total energy density  [g/cm^3]
%         p   = pressure              [dyn/cm^2]
%
%        Author: SBenruz 01/11
%
         
% units 
utsMdens_cgs = 6.1764e+17; % g cm^-3
utsPress_cgs = 5.5511e+38; % dyn cm^-2
m0           = 1.66e-24;   % g
utsN_cgs     = 1e+39;      % cm^-3

% get size
n = length(rho);

% construct energy
ene  = rho.*(epsl+1);

% change units
tab_nb  = rho *utsMdens_cgs/(utsN_cgs*m0);
tab_ene = ene *utsMdens_cgs;
tab_pre = pre *utsPress_cgs;

% write data
fprintf(1,'===> writing table %s ...',fname)

fid = fopen(fname, 'wt');
fprintf(fid,'#  %s\n',fname);
fprintf(fid,'#  Date: %s\n',date);
fprintf(fid,'#  Created with: writeTab1D.m (Author: S.Bernuzzi)\n');
fprintf(fid,'%d  <-- Number of lines\n',n);
fprintf(fid,'#\n');
fprintf(fid,'#            n_B [fm^{-3}]       rho [g/cm^3]        p [dyn/cm^2]\n');
fprintf(fid,'#\n');
for k=1:n
  fprintf(fid,'%d %2.16e %2.16e %2.16e\n',...
          k,tab_nb(k),tab_ene(k),tab_pre(k));
end
fprintf(1,' done\n')

% output
if nargout==3
  varargout = {tab_nb, tab_ene, tab_pre};
elseif nargout~=0
  error('3 o 0 output arguments required')
end
