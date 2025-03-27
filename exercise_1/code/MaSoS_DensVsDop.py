# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 18:45:46 2017

@author: Robert Fabbro

modifications: Karin Zojer
March 2023
"""

"""
 This Script and all the necessary functions were translated from the 
 corresponding MATLAB script for the first exersice of the lecture 'Modeling 
 and Simulation of Semiconductors'.
"""

import numpy as np
import matplotlib.pylab as plt
import sys # for controlled exit

# Module with all necessary functions for the calculations. 
# The module have to be copied into the individual functions or save folder.  
# A list of all functions and a short description are in the Exercise Skript.

import semiconductor_functions as sf

# To only write the function names without the 'sf.' prefix import the
# semiconductor_functions as followed:
    # from semiconductor_functions import *

##########################################################################
# -----------------------------------------------------------------------#
#      determination of Fermi levels and charge carrier densities        #
# -----------------------------------------------------------------------#
########################################################################## 

# hints for usage
# 
# + all spatial quantities are given in SI
# + temperatures are given in K
# + all energies are given in eV

plt.close('all')

#-------------------------------------------------------------------------
#--- parameter section ---------------------------------------------------
#-------------------------------------------------------------------------
 


k_eV = 1.3806504E-23 / 1.602176487E-19 # Boltzmann constant in eV/K

# external conditions
temperature = 300  # in K

# technical parameters
# number of intervals in DOS and energy
resolution = 1500

# threshold and maximal number of iterations for root finding
tolerance = 1e-4
max_RF_iter = 30

# ----- script control ------------
np.seterr('ignore') # Ignores Over/Underflow warnings/errors

# plot on logarithmic scale
#   useful for quantities that stretch across orders of magnitude, e.g.
#   charge densities and doping densities
plot_xlog = True
plot_ylog = False

# if plots do not show, toogle this switch to "True"
prompt_fig_show = False

# filenames for saving plots
filename_s1 = "chemical_potential_vs_Ndop"
filename_s2 = "relative_elec_density_vs_Ndop"
filename_s3 = "ionized_doping_density_vs_Ndop"

###############################################################################
###############################################################################

# Setting the Properties 

# doping
# ... doping density
dopant_density = 1e22   # Concentration of dopants in m^-3
# range of doping densities
dopmin = 1e15
dopmax = 1e25
N_steps = 500   # Number of steps between dopmin and dopmax 
# ... energy offset to band
E_doffset = 0.1
# ... doping type
doping_type = 'acceptor'

doping_densities_exp = np.linspace(np.log10(dopmin),np.log10(dopmax),N_steps) 
#doping_densities = np.logspace(dopmin,dopmax,N_steps)
doping_densities = np.power(10, doping_densities_exp)


# Initialize Storing Vectors

# Storing vectors for intrinsic and doped electron density
n = np.zeros_like(doping_densities)
n_i = np.zeros_like(doping_densities)
p = np.zeros_like(doping_densities)
trapden = np.zeros_like(doping_densities)

# Storing vectors for ionized dopants
Ndop_ionized = np.zeros_like(doping_densities)

# Storing vectors for charged traps
trap_ionized = np.zeros_like(doping_densities)

# Storing vectros for intrinsic and doped chemical potential
chemical_potential = np.zeros_like(doping_densities)
chemical_potential_i = np.zeros_like(doping_densities)

###############################################################################

# Assign Semiconducter Material. The Semiconductors to be chosen from are Si,
# Ge and GaAs. Anything else will create a artifitial materials. The properties
# can be changed in the routine AssignSemiconductor.

Semiconductor = 'Si'
EC, EV, m_n_eff, m_p_eff = sf.AssignSemiconductor(Semiconductor)

# EC...conduction band minimum; EV...valence band maximum;
# m_n_eff...effective electron mass; m_p_eff...effective hole mass 

# make assignments for dopants
if (doping_type == 'donor'):
    E_dopant = EC - E_doffset
    dopant_label = 'N$_D$'
    dopant_ion = '$^+$'
    dopant_marker = 'N'
elif (doping_type == 'acceptor'):
    E_dopant = EV + E_doffset
    dopant_label = 'N$_A$'
    dopant_ion = '$^-$'
    dopant_marker = 'P'
  
    
# traps
# ... trap density
trap_density = 1E22   # Concentration of traps in m^-3  1e22
# ... energy offset with respect to mid gap
# ... offset > 0 : above mid gap
# ... offset < 0 : below mid gap
E_toffset = -0.2#0.15 -(EC-EV)/2.0
# ... width of associated Gaussian
trap_sigma = 0.1 # in eV
# ... trap type
trap_type = 'acceptor'
    
# make assignments for traps
if (trap_density > 0):
    E_trap = (EC - EV)/2 + E_toffset
    trap_label = 'N$_T$'
    if (trap_type == 'donor'):
        trap_marker = 'N'
        # for a positively charged trap after filling
        trap_ion = '$^+$'
        # for a neutralized trap after filling
        # trap_ion = '$$'
    elif (trap_type == 'acceptor'):
        trap_marker = 'P'
        # for a negatively charged trap after filling
        trap_ion = '$^-$'
        # for a neutralized trap after filling
        # trap_ion = '$$'
    
###############################################################################

# Size of Energy Intervalls

E_min = EV - 0.5    # eV
E_max = EC + 0.5    # eV


                            
###############################################################################

#      evaluate intrinsic Fermi level numerically
#     --------------------------------------------------------------------
#     (a) cast charge neutrality condition into a form F(E,...) = 0 
#     (b) function F has to be provided, here chargeNeutralityIntrinsic()
#     (c) pass F as function handle fh, make sure that E is indicated as the
#         argument to be evaluated.

E_guess = (EC+EV)/2. + 0.2  # Guessing value for the FindRootNestedIntervals routine

    
for k in range(len(doping_densities)):
    
    dopant_density = doping_densities[k]
    
    
    # assemble DOS
    
    # Energy vector. Storing vectors for DOS and occupation
    energies, DOS, occupation = sf.InitializeEnergyAndDOS(E_min, E_max, EV, EC,resolution)
    
    # Initialize DOS_Admin. Capable of holding Information necessary to 
    # describe the nature of the states participating in the overall DOS
    DOS_Admin = sf.InitializeDOSAdministration(energies)
    
    # Adds DOS of the conduction band, adds empty vector for electron 
    # occupation
    DOS_Admin = sf.AddConductionBandToDOS(DOS_Admin,energies,EC,m_n_eff)
    
    # Adds DOS of the valence band, adds empty vector for hole 
    # occupation
    DOS_Admin = sf.AddValenceBandToDOS(DOS_Admin,energies,EV,m_p_eff)
    
    ## Add the Donor (N), Acceptor (P) or Neutral (0) level to DOS
    DOS_Admin  = sf.AddLevelToDOS(DOS_Admin,energies,
                                  dopant_density, E_dopant,dopant_marker)
    
    ## Add trap level if available
    if trap_density >0:
        
            DOS_Admin  = sf.AddGaussToDOS(DOS_Admin,energies,
                                  trap_density, E_trap, trap_sigma, trap_marker)
            
    # define an inline function here using Python lambda command
    # function_name = lambda argument : manipulate(argument)
    
    fh = lambda E: sf.chargeNeutralityIntrinsic(E ,EC,EV,m_n_eff,m_p_eff,temperature)
    
#    (d) employ root-finding algorithm to determine the chemical potential

    chemical_potential_i[k], num_iter, error = sf.FindRootNestedIntervals(fh,energies,E_guess,tolerance, max_RF_iter)


#   Calculation of the effective densities of state NV and NC

    n_t = sf.GetDensityInBand(chemical_potential_i[k],EC,m_n_eff, temperature)
    p_t = sf.GetDensityInBand(chemical_potential_i[k],EV,m_p_eff, temperature)
    n_i[k] = np.sqrt(n_t*p_t) 
    
    
#    ------------------------------------------------------------------
#     evaluate Fermi level numerically for a non-intrinsic system
#    ------------------------------------------------------------------
#
#     (a) cast charge neutrality condition into a form F(E,...) = 0 
#     (b) function F has to be provided, here chargeNeutrality()
#     (c) pass F as function handle fh, make sure that E is indicated as the
#         argument to be evaluated

    fh2 = lambda E: sf.chargeNeutrality(E,DOS_Admin,m_n_eff,m_p_eff,temperature)

#     (d) employ root-finding algorithm to determine the chemical potential

    chemical_potential[k], num_iter, error = sf.FindRootNestedIntervals(fh2,energies,chemical_potential_i[k],tolerance,max_RF_iter)

    n[k] = sf.GetDensityInBand(chemical_potential[k],EC,m_n_eff,temperature)
    p[k] = sf.GetDensityInBand(chemical_potential[k],EV,m_p_eff,temperature)

    if trap_density> 0:
        trapden[k] = sf.GetDensityInGauss(chemical_potential[k], DOS_Admin.Label[3], DOS_Admin.N[3], DOS_Admin.E_ref[3], DOS_Admin.Param[3], temperature)
     
    # this is only valid for acceptor-type dopants
    Ndop_ionized[k] = sf.GetDensityInLevel(chemical_potential[k],
                    DOS_Admin.E_ref[2],DOS_Admin.N[2],temperature)/DOS_Admin.N[2]
    # hence, for donors we need to correct as follows (positively charged)
    if doping_type == 'donor' :
            Ndop_ionized[k] = 1.0 - Ndop_ionized[k]
            
    # and check relative trap filling        
    # this is only valid for acceptor-type traps
    if trap_density> 0:
        trap_ionized[k] = sf.GetDensityInGauss(chemical_potential[k],
                        DOS_Admin.Label[3], 
                        DOS_Admin.N[3],
                        DOS_Admin.E_ref[3],
                        DOS_Admin.Param[3],
                        temperature)/DOS_Admin.N[3]
        # hence, for donors we need to correct as follows (positively charged)
        if trap_type == 'donor' :
            trap_ionized[k] = 1.0 - trap_ionized[k]
    
#################################################################################################################
#################################################################################################################

# Plotting of the Results #############
    
# set colors to distinguish semiconductors    

if Semiconductor == 'Si':
    col = 'red'
elif Semiconductor == 'GaAs':
    col = 'blue'
elif Semiconductor == 'Ge':
    col = 'green'
else:
    col = 'orange'

# Plot Data and set plot properties

#plt.style.use('seaborn-paper')  # Set 'seaborn-poster' for bigger plot and text size.
# type plt.sytle.available into console to see styles available.

#------------------------------------------------------------------------------
# plot chemical potential and intrinsic chemical potential as a function of temperature
    
fig1 = plt.figure(1) 

if plot_xlog:
    line_mu_i = plt.semilogx(doping_densities,chemical_potential_i,
                     linewidth=3,color=col,
                     linestyle='--',label='$\mu_i$ '+Semiconductor)
    line_mu = plt.semilogx(doping_densities,chemical_potential,
                   linewidth = 2,color=col,label='$\mu$ '+Semiconductor)
else:    
    line_mu_i = plt.plot(doping_densities,chemical_potential_i,
                     linewidth=3,color=col,
                     linestyle='--',label='$\mu_i$ '+Semiconductor)
    line_mu = plt.plot(doping_densities,chemical_potential,
                   linewidth = 2,color=col,label='$\mu$ '+Semiconductor)

# Set Labels
plt.xlabel('doping density / m$^{-3}$ ')
plt.ylabel('Energy / eV ')
plt.title('Chemical Potential vs doping density for ' + Semiconductor + ' at $ T= ' + str(temperature) + ' $K')

# Semiconductor Parameter section. If unwanted, comment out
plt.figtext(0.91, 0.5,
         Semiconductor+' Parameter:\n\n $E_{g}$ = '+str(EC)+' eV\n $m_{n}^{*}$= '+\
         str(m_n_eff)+' $m_{e}$\n $m_{p}^{*}$ = '+str(m_p_eff)+' $m_{e}$',
         bbox={'facecolor':'none', 'alpha':0.3, 'pad':5})

plt.legend()
plt.ylim(0,EC)

if prompt_fig_show:
    fig1.show() 

plt.savefig(filename_s1+".svg")
plt.savefig(filename_s1+".pdf")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# plot electron density function of doping densities

fig2 = plt.figure(2) # Electron Density

plot_ylog = False
plot_xlog = True

if plot_xlog:
  

    
    if plot_ylog:
        line_n_i = plt.semilogx(doping_densities,n_i,lw=1,color=col,ls='--',label='$n_i$ '+Semiconductor)
        line_n = plt.semilogx(doping_densities,n,lw=2,color=col,label='n '+Semiconductor)
        line_p = plt.semilogx(doping_densities,p,lw=2,color=[0,0.32,0.7],label='p '+Semiconductor)

        plt.yscale('log')
        plt.ylim(1E6,1E26)
    else:
        line_n_i = plt.semilogx(doping_densities,n_i/doping_densities,lw=1,color=col,ls='--',label='$n_i$ '+Semiconductor)
        line_n = plt.semilogx(doping_densities,n/doping_densities,lw=2,color=col,label='n '+Semiconductor)
        line_p = plt.semilogx(doping_densities,p/doping_densities,lw=2,color=[0,0.32,0.7],label='p '+Semiconductor)    

        plt.ylim(0,1.1)
else:
    line_n_i = plt.plot(doping_densities,n_i/doping_densities,lw=1,color=col,ls='--',label='$n_i$ '+Semiconductor)
    line_n = plt.plot(doping_densities,n/doping_densities,lw=2,color=col,label='n '+Semiconductor)
    line_p = plt.plot(doping_densities,p/doping_densities,lw=2,color=[0,0.32,0.7],label='p '+Semiconductor)

    plt.ylim(0,1)

line_reference = plt.plot()
# end KZ
plt.xlabel('doping density / m$^{-3}$ ')
plt.ylabel('charge density / $m^{-3}$')
plt.title('Charge Density vs doping density for ' + Semiconductor + ' at $ T= ' + str(temperature) + ' $K')

plt.figtext(0.91, 0.5,
         Semiconductor+' Parameter:\n\n $E_{g}$ = '+str(EC)+' eV\n $m_{n}^{*}$= '+\
         str(m_n_eff)+' $m_{e}$\n $m_{p}^{*}$ = '+str(m_p_eff)+' $m_{e}$',
         bbox={'facecolor':'none', 'alpha':0.3, 'pad':5})

plt.legend()


if prompt_fig_show:
    fig2.show() 
    
plt.savefig(filename_s2+".svg")
plt.savefig(filename_s2+".pdf")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# plot density of ionized dopants as function of temperature
fig3 = plt.figure(3) 


label_str = dopant_label+dopant_ion+'/'+dopant_label
trap_str = trap_label+trap_ion+'/'+trap_label

if plot_xlog:
    line_Ndopi = plt.semilogx(doping_densities,Ndop_ionized,lw=2,color=col,
                              label=label_str+' '+Semiconductor)
    
    if trap_density > 0:
        line_trapi = plt.semilogx(doping_densities,trap_ionized,
                                  lw=2,
                                  color='black',
                                  label=trap_str+' '+Semiconductor)
else:
    line_Ndopi = plt.plot(doping_densities,Ndop_ionized,lw=2,color=col,
                          label=label_str+' '+Semiconductor)
    if trap_density > 0:
        line_trapi = plt.semilogx(doping_densities,
                                  trap_ionized,
                                  lw=2,
                                  color='black',
                                  label=trap_str+' '+Semiconductor)

plt.xlabel('doping density / m$^{-3}$ ')
plt.ylabel(label_str)
plt.title('Number of Ionized Dopants vs doping density for ' + Semiconductor + ' at $ T= ' + str(temperature) + ' $K')
plt.legend()

plt.figtext(0.91, 0.5,
         Semiconductor+' Parameter:\n\n $E_{g}$ = '+str(EC)+' eV\n $m_{n}^{*}$= '+\
         str(m_n_eff)+' $m_{e}$\n $m_{p}^{*}$ = '+str(m_p_eff)+' $m_{e}$',
         bbox={'facecolor':'none', 'alpha':0.3, 'pad':5})

plt.ylim(0,1.2)

if prompt_fig_show:
    fig3.show() 

plt.savefig(filename_s3+".svg")
plt.savefig(filename_s3+".pdf")
    
