% script to generate tabs from fits

% -----------------------------------------------------------------

% the following tables sample poorly data at low densities used only for ATM
% and has 2000 pts for rho > 1e8 g / cm^3

%

% apr shib
% rem rho<1e-6 inaccurate but you need something for atm 
%     rho>1.5e-2 unphysical values p<0
rho=logspace(-16,-9.1,48); 
rho=[rho logspace(-9,-2,2000)];
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('apr','shib',rho);
writeTab1D('eos_apr_shibfit_2048.d',rho,ene./rho-1,p);

% rem the following energy density range apply to all sly and fps fits
rho=logspace(-16,-9.1,48);
rho=[rho logspace(-9,-1,2000)];

% sly shib
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('sly','shib',rho);
writeTab1D('eos_sly_shibfit_2048.d',rho,ene./rho-1,p);

% sly hps
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('sly','hps',rho);
writeTab1D('eos_sly_hpsfit_2048.d',rho,ene./rho-1,p);

% fps hps
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('fps','hps',rho);
writeTab1D('eos_fps_hpsfit_2048.d',rho,ene./rho-1,p);

% fps shib
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('fps','shib',rho);
writeTab1D('eos_fps_shibfit_2048.d',rho,ene./rho-1,p);

%


% -----------------------------------------------------------------


% the following table sample poorly data at low densities used only for ATM
% and has 500 pts for rho > 1e8 g / cm^3

%{

% apr shib
% rem rho<1e-6 inaccurate but you need something for atm 
%     rho>1.5e-2 unphysical values p<0
rho=logspace(-16,-9.1,12); 
rho=[rho logspace(-9,-2,500)];
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('apr','shib',rho);
writeTab1D('eos_apr_shibfit.d',rho,ene./rho-1,p);

%}


% -----------------------------------------------------------------


% the following are uniformaly spaced tables
% most of the pts are probabily not used during simulations 
% since they refer to too low densities

%{

% apr shib
% rem rho<1-6 inaccurate but you need something for atm 
%     rho>1.5e-2 p<0
rho=logspace(-15,-2,2048);
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('apr','shib',rho);
writeTab1D('eos_apr_shibfit.d',rho,ene./rho-1,p);

% sly shib
rho=logspace(-16,-1,2048);
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('sly','shib',rho);
writeTab1D('eos_sly_shibfit.d',rho,ene./rho-1,p);

% sly hps
rho=logspace(-16,-1,2048);
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('sly','hps',rho);
writeTab1D('eos_sly_hpsfit.d',rho,ene./rho-1,p);

% fps hps
rho=logspace(-16,-1,2048);
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('fps','hps',rho);
writeTab1D('eos_fps_hpsfit.d',rho,ene./rho-1,p);

% fps shib
rho=logspace(-16,-1,2048);
[p,DpDrho,DpDepsl,cs2,ene]=EoSColdAnalFits('fps','shib',rho);
writeTab1D('eos_fps_shibfit.d',rho,ene./rho-1,p);

%}


