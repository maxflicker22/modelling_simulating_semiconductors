# -*- coding: utf-8 -*-
"""
Created on Wed May  3 18:45:46 2017

@author: Robert Fabbro

modifications: Karin Zojer
March 2024
"""

"""
 This Script and all the necessary functions were translated from the 
 corresponding MATLAB script for the first exersice of the lecture 'Modeling 
 and Simulation of Semiconductors'.
"""


# %%
import numpy as np
import matplotlib.pylab as plt
from colorsys import rgb_to_hls
import seaborn as sns
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
 
np.seterr('ignore') # Ignores Over/Underflow warnings/errors

k_eV = 1.3806504E-23 / 1.602176487E-19 # Boltzmann constant in eV/K

# external conditions
T = 300  # in K

# technical parameters
# number of intervals in DOS and energy
resolution = 1500

# threshold and maximal number of iterations for root finding
tolerance = 1e-4
max_RF_iter = 30


# ----- script control ------------
# plot on logarithmic scale
plot_xlog = False
plot_ylog = False
# if plots do not show, toogle this switch to "True"
prompt_fig_show = False

# Semiconductor material
# initializue here cause of filename
#  choose from 'Si', 'Ge', 'GaAs'
semiconductor_string = 'MySemi1'


# filenames for saving plots
filename_s1 = "chemical_potential_vs_temp" + semiconductor_string
filename_s2 = "relative_elec_density_vs_temp" + semiconductor_string
filename_s3 = "ionized_doping_density_vs_temp" +semiconductor_string

###############################################################################
###############################################################################

# Setting the Properties 

# doping
# ... doping type
doping_type = 'donor'
# ... doping density
dopant_density = 1e21   # Concentration of dopants in m^-3
# ... energy offset to band
# ... n-doping: requires donors with small offset E_doffset to conduction band
# ... p-doping: requires acceptors with small offset E_doffset to valence band
E_doffset = 0.2


# traps

#   there are no traps in this example

# temperature
Tmin = 10      # Starting temperature in K
Tmax = 1000       # End temperature in K
N_steps = round((Tmax-Tmin)/2.0)   # Number of steps between Tmin and Tmax 

temperature = np.linspace(Tmin,Tmax,N_steps)    # Desired Temperature Range

# Initialize Storing Vectors

# vectors will store intrinsic and doped electron density
n = np.zeros_like(temperature)
p = np.zeros_like(temperature)
n_i = np.zeros_like(temperature)

# vectors will store ionized dopants
ND_ionized = np.zeros_like(temperature)

# vectors will store intrinsic and doped chemical potential
chemical_potential = np.zeros_like(temperature)
chemical_potential_i = np.zeros_like(temperature)

###############################################################################

# Assign Semiconducter Material. The Semiconductors to be chosen from are Si,
# Ge and GaAs. Anything else will create a artifitial materials. The properties
# can be changed in the routine AssignSemiconductor.

Semiconductor = semiconductor_string
EC, EV, m_n_eff, m_p_eff = sf.AssignSemiconductor(Semiconductor)

# EC...conduction band minimum; EV...valence band maximum;
# m_n_eff...effective electron mass; m_p_eff...effective hole mass 


# these lines communicate the position of the doping level
# and prepare some variables to help to generate labels for the graphical displays
if (doping_type == 'donor'):
    E_dopant = EC - E_doffset
    dopant_label = 'N$_D$'
    dopant_ion = '$^{+}$'
    dopant_marker = 'N'
elif (doping_type == 'acceptor'):
    E_dopant = EV + E_doffset
    dopant_label = 'N$_A$'
    dopant_ion = '$^{-}$'
    dopant_marker = 'P'


###############################################################################

#      let us build the density of states
#     --------------------------------------------------------------------
#     (a) define the boundaries of the considered energy interval E_min, E_max
#     (b) initialize the DOS_Admin container in which the nature of the DOS will be stored
#     (c) add DOS of conduction band 
#     (d) add DOS of valence band 
#     (e) add DOS of dopant level     
#

# Size of Energy Intervalls

E_min = EV - 0.5    # eV
E_max = EC + 0.5    # eV

   
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
DOS_Admin  = sf.AddLevelToDOS(DOS_Admin,energies,dopant_density, E_dopant,dopant_marker)


                            
###############################################################################

#      evaluate intrinsic Fermi level numerically
#     --------------------------------------------------------------------
#     (a) cast charge neutrality condition into a form F(E,...) = 0 
#     (b) function F has to be provided; here F is chargeNeutralityIntrinsic()
#     (c) pass F as function handle fh, make sure that E is indicated as the
#         argument to be evaluated.

E_guess = (EC+EV)/2. + 0.2  # Guessing value for the FindRootNestedIntervals routine

    
for k in range(len(temperature)):
    
    Temp = temperature[k]
    
    fh = lambda E: sf.chargeNeutralityIntrinsic(E ,EC, EV, m_n_eff, m_p_eff, Temp)
    
#    (d) employ root-finding algorithm to determine the chemical potential

    chemical_potential_i[k], num_iter, error = sf.FindRootNestedIntervals(fh,energies,E_guess,tolerance, max_RF_iter)


#   Calculation of the effective densities of state NV and NC
#   and the intrinsic electron density n_i    
    
    # check for division by zeros due to either n or p being zero
    # put in limit for T-> 0
    if error == 0:
        print(Temp,'K : prevented DIVISION by zero');
        chemical_potential_i[k] = (EC + EV)/2.;
        n_i[k]= 0;        
    else:
        n_t = sf.GetDensityInBand(chemical_potential_i[k],EC,m_n_eff, Temp)
        p_t = sf.GetDensityInBand(chemical_potential_i[k],EV,m_p_eff,Temp)
        n_i[k] = np.sqrt(n_t*p_t)                     
   
    
###############################################################################    
#    ------------------------------------------------------------------
#     evaluate Fermi level numerically for a non-intrinsic system
#    ------------------------------------------------------------------
#
#     (a) cast charge neutrality condition into a form F(E,...) = 0 
#     (b) function F has to be provided; here F is chargeNeutrality()
#     (c) pass F as function handle fh2, make sure that E is indicated as the
#         argument to be evaluated

    fh2 = lambda E: sf.chargeNeutrality(E,DOS_Admin,m_n_eff,m_p_eff,Temp)

#     (d) employ root-finding algorithm to determine the chemical potential

    chemical_potential[k], num_iter, error = sf.FindRootNestedIntervals(fh2,energies,chemical_potential_i[k],tolerance,max_RF_iter)

    n[k] = sf.GetDensityInBand(chemical_potential[k],EC,m_n_eff,Temp)
    p[k] = sf.GetDensityInBand(chemical_potential[k],EV,m_p_eff,Temp)

#     get  number of ionized dopants:
#     this is only valid for acceptor-type dopants
    ND_ionized[k] = sf.GetDensityInLevel(chemical_potential[k],
                    DOS_Admin.E_ref[2],DOS_Admin.N[2],Temp)/DOS_Admin.N[2]

    # hence, for donors we need to correct as follows (positively charged)
    #     determine how much dopant levels are occupied,  
    #     and substract that from total doping density 
    if doping_type == 'donor' :
            ND_ionized[k] = 1.0 - ND_ionized[k]    

  
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

line_mu_i = plt.plot(temperature,chemical_potential_i,lw=3,color=col,ls='--',label='$\mu_i$ '+Semiconductor)
line_mu = plt.plot(temperature,chemical_potential,lw=2,color=col,label='$\mu$ '+Semiconductor)
line_bandedge = plt.plot(temperature,EC*np.ones(np.shape(temperature)),lw=1,ls='dotted', color=col,label='$E_C$ '+Semiconductor)


# Set Labels
plt.xlabel('Temperature / K')
plt.ylabel('Energy /eV')
plt.title('Chemical Potential vs Temperature for ' + Semiconductor + ' ' + dopant_label + ' ' + str(dopant_density) + ' $m^{-3}$')

# Semiconductor Parameter section. If unwanted, comment out
plt.figtext(0.65, 0.25,
         Semiconductor+' Parameter:\n\n $E_{g}$ = '+str(EC)+' eV\n $m_{n}^{*}$= '+\
         str(m_n_eff)+' $m_{e}$\n $m_{p}^{*}$ = '+str(m_p_eff)+' $m_{e}$',
         bbox={'facecolor':'none', 'alpha':0.3, 'pad':5})

plt.legend()
plt.ylim(0,EC+0.1)

if prompt_fig_show:
    fig1.show()  
#plt.savefig(filename_s1+".svg")
plt.savefig(filename_s1+".pdf")

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# plot electron density function of temperature

fig2 = plt.figure(2) # Electron Density

line_n_i = plt.plot(temperature,n_i/dopant_density,lw=1,color=col,ls='--',
                    label='$n_i$ '+Semiconductor)
line_n = plt.plot(temperature,n/dopant_density,
                  linewidth = 2,
                  color = col,
                  label = 'n '+Semiconductor)
line_p = plt.plot(temperature,
                  p/dopant_density,
                  linewidth = 2,
                  color=col,
                  alpha = 0.5,
                  label = 'p '+Semiconductor)

line_reference = plt.plot()
# end KZ
plt.xlabel('Temperature / K')
plt.ylabel('density/' + dopant_label)
plt.title('Charge Density vs Temperature for ' + Semiconductor + ' ' + dopant_label + ' ' +\
          str(dopant_density) + ' $m^{-3}$')

plt.figtext(0.6, 0.3,
         Semiconductor+' Parameter:\n\n $E_{g}$ = '+str(EC)+' eV\n $m_{n}^{*}$= '+\
         str(m_n_eff)+' $m_{e}$\n $m_{p}^{*}$ = '+str(m_p_eff)+' $m_{e}$',
         bbox={'facecolor':'none', 'alpha':0.3, 'pad':5})

plt.legend()
plt.ylim(0,1.1)

if prompt_fig_show:
    fig2.show() 

#plt.savefig(filename_s2+".svg")
plt.savefig(filename_s2+".pdf")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# plot density of ionized dopants as function of temperature
fig3 = plt.figure(3) # Ionized Dopants

if plot_ylog:
    line_NDi = plt.semilogy(temperature,ND_ionized,lw=2,color=col,
                        label = dopant_label + dopant_ion + ' ' + Semiconductor)
else:
    line_NDi = plt.plot(temperature,ND_ionized,lw=2,color=col,
                        label = dopant_label + dopant_ion + ' ' + Semiconductor)

plt.xlabel('Temperature /K')
#plt.ylabel('N$_{D^+}$/N$_D$ ')
plt.ylabel(dopant_label + dopant_ion +' /' + dopant_label)
plt.title('Number of Ionized Dopants vs Temperature for ' + Semiconductor +\
           ' ' + dopant_label + ' ' + str(dopant_density) + ' $m^{-3}$')
plt.legend()

plt.figtext(0.55, 0.26,
         Semiconductor+' Parameter:\n\n $E_{g}$ = '+str(EC)+' eV\n $m_{n}^{*}$= '+\
         str(m_n_eff)+' $m_{e}$\n $m_{p}^{*}$ = '+str(m_p_eff)+' $m_{e}$',
         bbox={'facecolor':'none', 'alpha':0.3, 'pad':5})

plt.ylim(1E-14,1.2)

if prompt_fig_show:
    fig1.show() 

#plt.savefig(filename_s3+".svg")
plt.savefig(filename_s3+".pdf")
# %%


