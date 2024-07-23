function [ rho, epsl, pre ] = readTab1D( fname )
%READTAB1D imports data from 1D tables to matlab
%
%         [ rho, epsl, pre ] = readTab1D( fname )
%
%         Tables must contain the line 
%           172  <-- Number of lines
%         and 4 cols
%           index  n_B [fm^{-3}]  e [g/cm^3]   p [dyn/cm^2]
%
%         n_B = baryon number density, 
%         e   = total energy density,
%         p   = pressure.
%
%         Lines starting with # are ignored
%
%         The script returns rho (rest-mass density), epsl
%         (specific energy) and p (pressure) in dimensionless unit:
%                     c = G = Msun = 1
%

% units 
utsMdens_cgs = 6.1764e+17; % g cm^-3
utsPress_cgs = 5.5511e+38; % dyn cm^-2
m0           = 1.66e-24;   % g
utsN_cgs     = 1e+39;      % cm^-3

% import data
fprintf(1,'file : %s\n', fname)

fid = fopen( fname );

n=1;
nlines=2;
while n<=nlines
  sline = fgetl(fid);
  if ~ischar(sline); break; end;
  if( strcmp(sline(1),'#') )
    fprintf(1,'infos: %s\n',sline);
    continue; 
  end
  tmpi = strfind(sline,'<-- Number of lines');
  if ~isempty(tmpi) 
    nlines = str2num(sline(1:tmpi-1));
    continue;
  end
  tmpd = str2num(sline);
  tab_nB(n) = tmpd(2); tab_ene(n) = tmpd(3); tab_pre(n) = tmpd(4);
  n = n+1;
end

if( (n-1)~=nlines )
  fprintf(1,'error [readTab1D.m]: too many lines read');
  return
end
fprintf(1,'lines: %d\n', nlines);

% change units
rho = tab_nB*utsN_cgs*m0/utsMdens_cgs;
ene = tab_ene/utsMdens_cgs;
pre = tab_pre/utsPress_cgs;
epsl = ene./rho-1.;

