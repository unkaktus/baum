## units
hplanck_cgs = 6.62606876e-27
clight      = 2.99792458e+10
Msun        = 1.98892e+33
Length_cgs  = 0.5*2.953250077*1e5
Length_km   =  0.5*2.953250077

Mass_cgs    = Msun
Time_cgs    = Length_cgs/clight*1e3
Volume_cgs  = Length_cgs**3
Energy_cgs  = Msun*clight*clight
Edens_cgs   = Energy_cgs/Volume_cgs;
Mdens_cgs   = Msun/Volume_cgs
rho_nuc     = 2.3e14

Press_cgs   = Msun/(Time_cgs*Time_cgs*Length_cgs)*1e6;
Energy_MeV   = Energy_cgs / 1.602176565e-6

mbar_cgs = 1.67262192e-24
mbar_BAM = 1.67262192e-24/Msun
