% Script to test fits against tables 

% choose table 
% 'eos_akmalpr.d' 'eos_fps.d' 'eos_sly4.d'  'eos_sly4_hp.d'
tabname = '../eos_akmalpr.d';

% choose eos 
% 'fps'  'sly'  'apr'
eosname = 'apr';

% choose fit
% 'hp' 'hps' 'shib'
fitname = 'shib';

% set limits for rho in figs
rholim=[1.6e-8 5e-1];

% save figs?
sv='no';

% do not change below here

if strcmp(fitname,'shib')
  warning('shib fits not accurate for rho<10^{10} g/cm^3 (1.6 10^{-8})');
end

% load table (return dimensionless units!)
[trho,tepsl,tpre] = readTab1D(tabname);
tene=trho.*(tepsl+1);

% compute fit
[pre,DpDrho,DpDepsl,cs2,ene,epsl,DepslDrho,D2epslD2rho]=...
    EoSColdAnalFits(eosname,fitname,trho);

% compare epsl Vs rho
fh(1)=figure; fname{1}='epsl_rho';
subplot(3,1,1:2)
loglog( trho,tepsl,'b-', ...
        trho,epsl,'rx-' )
xlim(rholim); 
ylim([1e-4 10])
%xlabel('\rho')
ylabel('\epsilon') 
title(sprintf('EoS: %s  Fit: %s',eosname,fitname))
legend('tab','fit','Location','SouthEast')
subplot(3,1,3)
semilogx( trho, 100*abs(tepsl-epsl)./tepsl, 'r-' )
xlim(rholim)
ylabel('Diff.[%]')
xlabel('\rho')

% compare pre Vs rho
fh(2)=figure; fname{2}='p_rho';
subplot(3,1,1:2)
loglog( trho,tpre,'b-', ...
        trho,pre,'rx-' )
xlim(rholim); 
ylim([1e-12 1])
%xlabel('\rho')
ylabel('P') 
title(sprintf('EoS: %s  Fit: %s',eosname,fitname))
legend('tab','fit','Location','SouthEast')
subplot(3,1,3)
semilogx( trho, 100*abs(tpre-pre)./tpre, 'r-' )
xlim(rholim)
ylabel('Diff.[%]')
xlabel('\rho')

% compare pres Vs ene
fh(3)=figure; fname{3}='p_e';
subplot(3,1,1:2)
loglog( tene,tpre,'b-', ...
        ene,pre,'rx-' )
%xlabel('e')
ylabel('P') 
title(sprintf('EoS: %s  Fit: %s',eosname,fitname))
legend('tab','fit','Location','SouthEast')
subplot(3,1,3)
semilogx( tene, 100*abs(tpre-pre)./tpre, 'r-' )
ylabel('Diff.[%]')
xlabel('e')

% save figs
if strcmp(lower(sv),'yes')
  for f=1:length(fh)
    saveas(fh(f),[fname{f},'_',eosname,'_',fitname,'.eps'], 'psc2')
  end
end