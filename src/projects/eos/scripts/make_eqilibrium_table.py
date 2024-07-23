##   Script for computing the beta equilibrum configuration in every point of the EoS table as in Perego et. at. 2018 (https://arxiv.org/abs/1903.07898)
##   Needed for computing the emissivities in the M1 neutrino transport
##   Give as an input the path to the EoS file in .d format

import numpy as np
import units
import h5py
from argparse import ArgumentParser
import scipy.optimize

parser = ArgumentParser()
parser.add_argument("-eos", dest="eosfile",  required=True)
args = parser.parse_args()

eosfile = args.eosfile
input_name = eosfile.partition('.')[0]
out_fname_txt = input_name + "_equilibrium.d"
out_fname_h5  = input_name + "_equilibrium.h5"

m_b = 1.6588056935297748e-24 / units.Mass_cgs  ## baryon mass in solar masses
m_n = 939.565413                               ## neutron mass in MeV, c = 1

h = units.hplanck_cgs*units.Energy_MeV/units.Energy_cgs;  ## Planck const in MeV.s
c = units.clight;                                         ## speed of light in cm/s


########################
## table interpolator ##
########################

print("Loading EoS table from", args.eosfile)
eos = np.loadtxt(args.eosfile)

eos_rho = np.unique(eos[:,1])*m_b*1e+39/units.Mdens_cgs*units.Mass_cgs
eos_T   = np.unique(eos[:,0])
eos_Y   = np.unique(eos[:,2])

nrho = len(eos_rho)
nT   = len(eos_T)
nYe  = len(eos_Y)

eos_epsl = eos[:,8]/units.Energy_MeV/m_b - 1
eos_mu_e = eos[:,6] - eos[:,5]
eos_mu_n = eos[:,4] + m_n
eos_mu_p = eos[:,5] + eos_mu_n

eos_epsl = np.reshape(eos_epsl, (nT, nrho, nYe))
eos_mu_e = np.reshape(eos_mu_e, (nT, nrho, nYe))
eos_mu_n = np.reshape(eos_mu_n, (nT, nrho, nYe))
eos_mu_p = np.reshape(eos_mu_p, (nT, nrho, nYe))

eos_mu_nue = eos_mu_e + eos_mu_p - eos_mu_n
eos_mu_nua = -eos_mu_nue

eos_logrho = np.log10(eos_rho)
eos_logT   = np.log10(eos_T)

eos_rho_max = eos_rho[-1]
eos_rho_min = eos_rho[0]

eos_T_max = eos_T[-1]
eos_T_min = eos_T[0]

eos_Ye_max = eos_Y[-1]
eos_Ye_min = eos_Y[0]

dxi = 1/(eos_logT[1]-eos_logT[0])
dyi = 1/(eos_logrho[1]-eos_logrho[0])
dzi = 1/(eos_Y[1]-eos_Y[0])

dxyi = dxi*dyi
dxzi = dxi*dzi
dyzi = dyi*dzi
dxyzi = dxyi*dzi


def interp_3d(x, y, z, data):  # taken from BAM

    ix = 1+np.floor((x-eos_logT[0]   -1e-10)*dxi)
    iy = 1+np.floor((y-eos_logrho[0] -1e-10)*dyi)
    iz = 1+np.floor((z-eos_Y[0]      -1e-10)*dzi)

    ix = ix.astype(int)
    iy = iy.astype(int)
    iz = iz.astype(int)

    delx = eos_logT[ix]-x
    dely = eos_logrho[iy]-y
    delz = eos_Y[iz]-z

    fh0 = data[ix,     iy,     iz  ]
    fh1 = data[ix-1,   iy,     iz  ]
    fh2 = data[ix,     iy-1,   iz  ]
    fh3 = data[ix,     iy,     iz-1]
    fh4 = data[ix-1,   iy-1,   iz  ]
    fh5 = data[ix-1,   iy,     iz-1]
    fh6 = data[ix,     iy-1,   iz-1]
    fh7 = data[ix-1,   iy-1,   iz-1]

    a1 = fh0;
    a2 = (fh1 - fh0)*dxi;
    a3 = (fh2 - fh0)*dyi;
    a4 = (fh3 - fh0)*dzi;
    a5 = (fh4 - fh1 - fh2 + fh0)*dxyi;
    a6 = (fh5 - fh1 - fh3 + fh0)*dxzi;
    a7 = (fh6 - fh2 - fh3 + fh0)*dyzi;
    a8 = (fh7 - fh0 + fh1 + fh2 + fh3 - fh4 - fh5 - fh6)*dxyzi;

    f = a1 + a2*delx + a3*dely + a4*delz + a5*delx*dely + a6*delx*delz + a7*dely*delz + a8*delx*dely*delz

    return f


#####################
## Fermi integrals ##
#####################

def F2(eta):

    if (eta>1e-3):
        return 	(eta*eta*eta/3+3.2899*eta)/(1.0-np.exp(-1.8246*eta))
    else:
        return 2*np.exp(eta)/(1+0.1092*np.exp(0.8908*eta))

def F3(eta):

    if (eta>1e-3):
        return (eta*eta*eta*eta/4.0+4.9348*eta**2+11.3644)/(1.0+np.exp(-1.9039*eta))
    else:
        return 6.0*np.exp(eta)/(1.0+0.0559*np.exp(0.9069*eta))



#################################
## Number and energy densities ##
#################################

rho_lim_nue = 5e5 / units.Mdens_cgs
rho_lim_nux = 5e5 / units.Mdens_cgs

def Y_nue(rho,T,Ye):
    mu_nu = interp_3d(np.log10(T), np.log10(rho), Ye, eos_mu_nue)
    return 4*np.pi*m_b/(rho*h**3*c**3) * units.Volume_cgs * T*T*T * F2(mu_nu/T) * np.exp(-rho_lim_nue/rho)

def Y_nua(rho,T,Ye):
    mu_nu = interp_3d(np.log10(T), np.log10(rho), Ye, eos_mu_nua)
    return 4*np.pi*m_b/(rho*h**3*c**3) * units.Volume_cgs * T*T*T  * F2(mu_nu/T) * np.exp(-rho_lim_nue/rho)

def Z_nue(rho,T,Ye):
    mu_nu = interp_3d(np.log10(T), np.log10(rho), Ye, eos_mu_nue)
    return 4*np.pi/(h**3*c**3) * units.Volume_cgs/units.Energy_MeV * T*T*T*T  * F3(mu_nu/T) * np.exp(-rho_lim_nue/rho)

def Z_nua(rho,T,Ye):
    mu_nu = interp_3d(np.log10(T), np.log10(rho), Ye, eos_mu_nua)
    return 4*np.pi/(h**3*c**3) * units.Volume_cgs/units.Energy_MeV * T*T*T*T * F3(mu_nu/T) * np.exp(-rho_lim_nue/rho)

def Z_nux(rho,T,Ye):
    return 4*np.pi/(h**3*c**3) * units.Volume_cgs/units.Energy_MeV * T*T*T*T  * F3(0)       * np.exp(-rho_lim_nux/rho)

def energy(rho,T,Ye):
    return  rho*interp_3d(np.log10(T), np.log10(rho), Ye, eos_epsl)


#######################
## Function to solve ##
#######################

def bound_Y(x):
    if (x>eos_Ye_max):
        return eos_Ye_max
    elif (x<eos_Ye_min):
        return eos_Ye_min
    else:
        return x

def bound_T(x):
    if (x>eos_T_max):
        return eos_T_max
    elif (x<eos_T_min):
        return eos_T_min
    else:
        return x

def f(x, args):  ## args[0]=rho  args[1]=e  args[2]=Ye 

    res=np.zeros(2)

    res[0] = x[0]-args[2] + Y_nue(args[0],x[1],x[0]) - Y_nua(args[0],x[1],x[0])

    res[1] = energy(args[0],x[1],x[0]) - args[1] + Z_nue(args[0],x[1],x[0]) + Z_nua(args[0],x[1],x[0])  + 4*Z_nux(args[0],x[1],x[0])

    if (res[0]==1e50 or res[1]>1e50):
        print (x, args)

    return res

def F(x, args): ## args[0]=rho  args[1]=e  args[2]=Ye

    res = f([bound_Y(x[0]),bound_T(x[1])],args)

    if (x[0]>eos_Ye_max):
        res += (f([eos_Y[-1],bound_T(x[1])],args) - f([eos_Y[-2],bound_T(x[1])],args)) / (eos_Y[-1] - eos_Y[-2]) * (x[0] - eos_Y[-1])
   
    if (x[0]<eos_Ye_min):
        res += (f([eos_Y[1],bound_T(x[1])],args) - f([eos_Y[0],bound_T(x[1])],args)) / (eos_Y[1] - eos_Y[0]) * (x[0] - eos_Y[0])

    if (x[1]>eos_T_max):
        res += (f([bound_Y(x[0]),eos_T[-1]],args) - f([bound_Y(x[0]),eos_T[-2]],args)) / (eos_T[-1] - eos_T[-2]) * (x[1] - eos_T[-1])
    
    if (x[1]<eos_T_min):
        res += (f([bound_Y(x[0]),eos_T[1]],args) - f([bound_Y(x[0]),eos_T[0]],args)) / (eos_T[1] - eos_T[0]) * (x[1] - eos_T[0])

    res[0] /= args[2]
    res[1] /= args[1]

    return res


##############################
## solve for Ye_eq and T_eq ##
##############################
Yeq     = np.zeros((nYe,nT,nrho))
Teq     = np.zeros((nYe,nT,nrho))
success = np.zeros((nYe,nT,nrho))
err_Ye  = np.zeros((nYe,nT,nrho))
err_T   = np.zeros((nYe,nT,nrho))

out = h5py.File(out_fname_h5, "w")

out.create_dataset("rho", (nrho,), dtype=float)
out['rho'][:] = eos_rho

out.create_dataset("T", (nT,), dtype=float)
out['T'][:] = eos_T

out.create_dataset("Ye", (nYe,), dtype=float)
out['Ye'][:] = eos_Y

out.create_dataset("Yeq",     (np.shape(Yeq)),     dtype=float)
out.create_dataset("Teq",     (np.shape(Teq)),     dtype=float)
out.create_dataset("success", (np.shape(success)), dtype=float)
out.create_dataset("err_Ye",  (np.shape(err_Ye)),  dtype=float)
out.create_dataset("err_T",   (np.shape(err_T)),   dtype=float)

out.close()

for i in range(nYe):
    for j in range(nT):
        for k in range(nrho):

            rho = eos_rho[k] 
            T   = eos_T[j]
            Ye  = eos_Y[i]
            
            u = rho*interp_3d(np.log10(T), np.log10(rho), Ye, eos_epsl)
            
            if (k==0):
                x0 = [Ye, T]
            else:
                x0=[Yeq[i,j,k-1],Teq[i,j,k-1]]
            args = [rho, u, Ye]

            eq = scipy.optimize.root(F, x0, args=(args), tol=1e-6)
  
#            if (((eq.fun[0]**2 + eq.fun[1]**2)>1e-6)):
#                eq = scipy.optimize.root(F, x0, args=(args), tol=1e-6, method='hybr')

#            if (eq.success==False):
#                eq = scipy.optimize.root(F, x0, args=(args), tol=1e-6, method='lm')

#            if ((eq.fun[0]**2 + eq.fun[1]**2)>1e-6):
#                eq = scipy.optimize.root(F, x0, args=(args), tol=1e-6, method='krylov')

#            if ((eq.fun[0]**2 + eq.fun[1]**2)>1e-6):
#                eq = scipy.optimize.root(F, x0, args=(args), tol=1e-6, method='anderson')

#            if ((eq.fun[0]**2 + eq.fun[1]**2)>1e-6):
#                eq = scipy.optimize.root(F, x0, args=(args), tol=1e-6, method='df-sane')


            # save result
            Yeq[i,j,k] = eq.x[0]
            Teq[i,j,k] = eq.x[1]
            err_Ye[i,j,k] = eq.fun[0]
            err_T[i,j,k]  = eq.fun[1]

            if ((eq.fun[0]**2 + eq.fun[1]**2)<1e-6):
                success[i,j,k]=True
            else:
                success[i,j,k]=False

            print("Ye =", Ye, "  T =", T, "  rho =", rho, eq.fun[0],  eq.fun[1], Yeq[i,j,k], Teq[i,j,k], success[i,j,k])
#        print("Ye =", Ye, "  T =", T)


    out = h5py.File(out_fname_h5, "r+")

    out['Yeq'][:] = Yeq
    out['Teq'][:] = Teq
    out['success'][:] = success
    out['err_Ye'][:] = err_Ye
    out['err_T'][:] = err_T

    out.close()




##################
## write output ##
##################
Ye, T, rho = np.meshgrid(eos_Y, eos_T, eos_rho, indexing='ij')

## write in txt ##
f = open(out_fname_txt, 'w')

print("Writing table in text file \n")
for y in range(nT):   #loop over T
    for t in range(nrho):  #loop over rho

        to_print = np.array([T[:,y,t], rho[:,y,t], Ye[:,y,t], Yeq[:,y,t], Teq[:,y,t]])

        to_print = np.transpose(to_print)

        for r in range(nYe):  #loop over rho

            line = to_print[r,:]

            line = str(list(line))
            line = line.replace(']', '')
            line = line.replace('[', '')
            line = line.replace(',', '')
            line = line.replace('e', 'E')

            f.write(line + '\n')

f.close()

print("Table has been written in:")
print(out_fname_txt)
print(out_fname_h5)
print("\n")
