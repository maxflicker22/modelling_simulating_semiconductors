# -*- coding: utf-8 -*-
"""
adapted for Lecture Modelling
and Simulation of Semiconductors 2019

by:              Karin Zojer   
original author: Robert Fabbro

"""

# Module with all necessary Functions for Excercise 1 for the Lecture Modelling
# and Simulation of Semiconductors.



import numpy as np
import DOS_Admin
#import scipy as sci
from scipy.integrate import quad

###############################################################################
###############################################################################
def AssignSemiconductor(Semiconductor):
    
    
# Intrinsic Semiconductor Properties
#
#
# input:   label of compound (4 character string)
# output:  E_C     .. conduction band minimum / eV
#          E_V     .. valence band minimum / eV
#          m_n_eff .. effective electron mass / m_e
#          m_h_eff .. effective hole mass / m_e
#
#

# for a collection of material parameters, see
# http://lampx.tugraz.at/~hadley/psd/L3/L3.php
# https://ecee.colorado.edu/~bart/book/book/append/append3.htm
#
#
# within this example we work with the following approximation:
#
# For the valence and conduction band, ONE isotropic parabolic band with
# ONE associated effective mass is assumed
#
# i.e., we disgard
# + that the effective mass depends on the direction of motion, as bands do
#   not necessarily adopt the same curvature along all possible k
#   directions; longitudinal and transversal effective masses are not
#   discriminated
#   ml* = 1.64  ;  mt* = 0.082
#
# + the presence of multiple bands that share the same maximum or minimum
#   energy, e.g., the heavy and light hole bands of silicon
#   Si: mlh* = 0.044;  mhh* = 0.28 
#
# Consequently, the values of effective masses given below describe best a
# square-root-shaped DOS containing all contributing bands and directions


# list of effective masses (in units of electron mass)
#   element          Ge           Si           GaAs
#   -------------+--------+-------------+-------------
#   electron        0.55        1.08           0.067  (Si = 1.18?)
#   hole            0.37        0.81           0.45
#   hole (heavy)

#
# list of band gaps (in eV)
#   element          Ge           Si           GaAs
#   -------------+--------+-------------+-------------
#   band gap        0.66         1.12           1.424

    

    E_V = 0.  # eV
    E_C = 0.
    m_n_eff = 1.
    m_p_eff = 1.
    


    if Semiconductor == 'Si':
        
        m_n_eff = 1.18      # in units of electron mass
        m_p_eff = 0.81      # in units of electron mass
        E_C = E_V + 1.12    # eV
        
    elif Semiconductor == 'Ge':
        
        m_n_eff = 0.55      # in units of electron mass
        m_p_eff = 0.37      # in units of electron mass
        E_C = E_V + 0.66    # eV
        
    elif Semiconductor == 'GaAs':
        
        m_n_eff = 0.067     # in units of electron mass
        m_p_eff = 0.45      # in units of electron mass
        E_C = E_V + 1.424   # eV
     
    else:        # Artificial Semiconductor
        
        m_n_eff = 0.9       # in units of electron mass
        m_p_eff = 0.1       # in units of electron mass
        E_C = E_V + 2.4     # eV

    return E_C, E_V, m_n_eff, m_p_eff


###############################################################################
###############################################################################


def InitializeEnergyAndDOS(E_min, E_max, E_V, E_C, resolution):
    
# provides column vectors for: 
#         energies
#         DOS (returned zero-valued)

# input: 
#        Emin .. lowest considered energy
#        Emax .. highest considered energy
#        E_V  .. energy of valence band maximum
#        E_C  .. energy of conduction band minimum
#        resolution .. number of energy intervals

# equally spaced energy interval
# EnergyInt = linspace(E_min, E_max, resolution).';

#   produce an EnergyInterval being 
#   ... logarithmically divided in energy regions with bands
#   ... linearly divided in the gap region    

    number_of_gap_points = np.round(resolution*(E_C-E_V)/(E_max-E_min)).astype(int)
    number_of_val_points = np.round(resolution*(E_V-E_min)/(E_max-E_min)).astype(int)
    number_of_con_points = np.round(resolution*(E_max-E_C)/(E_max-E_min)).astype(int)
    
    number_of_gap_points = number_of_gap_points.astype(int)
    number_of_val_points = number_of_val_points.astype(int)
    number_of_con_points = number_of_con_points.astype(int)

    EnergyInt_V = E_V - (E_V-E_min)*np.logspace(-1,-5, base = 10, num = number_of_val_points)
    EnergyInt_Gap = np.linspace(E_V, E_C, number_of_gap_points)
    EnergyInt_C = E_C + (E_max-E_C)*np.logspace(-5,-1,number_of_con_points)
    
    # hstack for combining the vectors
    EnergyInt = np.hstack([EnergyInt_V,EnergyInt_Gap,EnergyInt_C])

    DOS = np.zeros_like(EnergyInt)
    occup_vec = np.zeros_like(EnergyInt)
    
    return EnergyInt, DOS, occup_vec


###############################################################################
###############################################################################


def InitializeDOSAdministration(energy_intervals):
    
# Initializes the DOS administration for storing all the necessary density
# if states
# returns a DOS_Admin as data record    
    
    DOS_Admin.Energies = energy_intervals;
    DOS_Admin.Label = 'em'  # em = empty
    DOS_Admin.Type = 'N'

    DOS_Admin.E_ref = 0 
    DOS_Admin.Param = 0
    DOS_Admin.N = 0         # total density or effective density of states

    return DOS_Admin


###############################################################################
###############################################################################


def AddContributionToDOS(DOS_admin, DOS, Label, Type, E_ref,Param,N_DOS):
    
# stores information on new DOS contribution into DOS_admin

# input  DOS_admin .. previous DOS_admin
#        DOS       .. vector containing the same number of bins as energy
#                     interval
#        label     .. 'D', 'A','em' 
#                     will be assigned by AddXXXToDOS 
#                               (XXX <> 'contribution')
#        type      .. 'N', 'P' 
#                     indicator whether states become negatively or  
#                     positively charged upon filling
#        E_ref     .. position of contribution in energy interval (in eV)
#        param     .. optional additional parameter, e.g., Gaussian width
#        N_DOS     .. total density of states (in inverse cubic meters)
#       

# output DOS_admin .. updated DOS_admin
    
    # check for valid entries in label and type
    
    
#     is this the first entry in DOS?
#     yes: overwrite dummy entry
#     no: expand list and set values

    if np.any(DOS_admin.Label == 'em'):
        DOS_admin.Energies = DOS
        DOS_admin.Label = Label
        DOS_admin.Type = Type

        DOS_admin.E_ref = E_ref
        DOS_admin.Param = Param
        DOS_admin.N = N_DOS
        
    # Check if there is only one entry. If yes the Dimensions of the entries 
    # have to be corrected to append with other DOS
    elif np.any(DOS_admin.Label != 'em') and np.size(DOS_admin.Label) < 2:
        
        DOS_admin.Energies = np.append([[DOS_admin.Energies]], [[DOS]], 
                                       axis = 0)
        DOS_admin.Label = np.append([[DOS_admin.Label]],[[Label]], 
                                       axis = 0)
        DOS_admin.Type = np.append([[DOS_admin.Type]], [[Type]], 
                                       axis = 0)

        DOS_admin.E_ref = np.append([[DOS_admin.E_ref]], [[E_ref]], 
                                       axis = 0)
        DOS_admin.Param = np.append([[DOS_admin.Param]], [[Param]], 
                                       axis = 0)
        DOS_admin.N = np.append([[DOS_admin.N]], [[N_DOS]], 
                                       axis = 0)
    
    else:
        DOS_admin.Energies = np.append(DOS_admin.Energies, [[DOS]], 
                                       axis = 0)
        DOS_admin.Label = np.append(DOS_admin.Label, [[Label]], 
                                       axis = 0)
        DOS_admin.Type = np.append(DOS_admin.Type, [[Type]], 
                                       axis = 0)

        DOS_admin.E_ref = np.append(DOS_admin.E_ref, [[E_ref]], 
                                       axis = 0)
        DOS_admin.Param = np.append(DOS_admin.Param, [[Param]], 
                                       axis = 0)
        DOS_admin.N = np.append(DOS_admin.N, [[N_DOS]], 
                                       axis = 0)
    
    return DOS_admin


###############################################################################
###############################################################################



def AddConductionBandToDOS(DOS_admin, E, E_C, eff_mass):
    
    
#     adds DOS of conduction band to a given DOS
# provides vector to store the hole occupation
#
# input: DOS .. previous DOS vector
#        E   .. energy vector
#        E_C .. minimum of conduction band   (eV)
#        eff_mass .. effective electron mass (in units of me)
#
# output: DOS vector with number of states per energy interval 
#         corresponding to requested Gaussian added
#         occupation vector (zero-valued)

# considering the parameters
   me = 9.11e-31 # kg
   h = 6.626E-34 # SI
   q = 1.602176565e-19 # SI needed to convert energies from eV to SI
   
# one arrives at a prefactor

    # bart, pierret, hadley
   prefactor = 8*np.pi*np.sqrt(2)*h**(-3)*me**(3/2)*q**(3/2)

#   when energies given in SI units (J)
#      prefactor = 8*pi*sqrt(2)*h^(-3)*me^(3/2)  ;
#   when energies given in eV
#      prefactor = 8*pi*sqrt(2)*h^(-3)*me^(3/2)*q^(3/2) ;
#
#   this is, because we would like to provide energies in eV and
#   plot the resulting DOS with respect to energy intervals given in eV
#
#   dZ/dE |_eV = dZ/dE |_SI * (dE|_SI)/(dE|_eV) = q dZ/dE |_SI
#   sqrt(E|_SI) = sqrt ( q E|_eV) = sqrt(q) * sqrt(E|_eV)

    
    
   energies_above_E_C = E > E_C
    
   DOS = np.zeros_like(E)

   DOS[energies_above_E_C] = 1*prefactor*eff_mass**(3/2)*np.sqrt(E[energies_above_E_C]-E_C)
   
   DOS_admin = AddContributionToDOS(DOS_admin, DOS, 'CB','N',E_C,0,0)
    
   return DOS_admin



###############################################################################
###############################################################################



def AddValenceBandToDOS(DOS_admin, E, E_V, eff_mass):
    
    
# adds DOS of valence band to a given DOS
# provides vector to store the hole occupation
#
# input: E   .. energy vector
#        E_V .. maximum of valence band   (eV)
#        eff_mass .. effective hole mass (in units of me)
#
# output: entry to DOS with
#         vector with number of states per energy interval 
#         corresponding to requested Gaussian 
#         labeled as VB = valence band
#         type = 'P' for positive charges

# considering the parameters
   me = 9.11e-31 # kg
   h = 6.626E-34 # SI
   q = 1.602176565e-19 # SI needed to convert energies from eV to Si
   
# one arrives at a prefactor

    # bart, pierret, hadley
   prefactor = 8*np.pi*np.sqrt(2)*h**(-3)*me**(3/2)*q**(3/2)


#   when energies given in SI units (J)
#      prefactor = 8*pi*sqrt(2)*h^(-3)*me^(3/2)  ;
#   when energies given in eV
#      prefactor = 8*pi*sqrt(2)*h^(-3)*me^(3/2)*q^(3/2) ;
#
#   this is, because we would like to provide energies in eV and
#   plot the resulting DOS with respect to energy intervals given in eV
#
#   dZ/dE |_eV = dZ/dE |_SI * (dE|_SI)/(dE|_eV) = q dZ/dE |_SI
#   sqrt(E|_SI) = sqrt ( q E|_eV) = sqrt(q) * sqrt(E|_eV)
    
   energies_below_E_V = E < E_V
    
   DOS = np.zeros_like(E)

   DOS[energies_below_E_V] = 1*prefactor*eff_mass**(3/2)*np.sqrt(E_V-E[energies_below_E_V])
   
   DOS_admin = AddContributionToDOS(DOS_admin, DOS, 'VB','P',E_V,0,0)
    
#   occ_vector = np.zeros_like(E)
    
   return DOS_admin




###############################################################################
###############################################################################




def AddLevelToDOS( DOS_admin,E, N_G, E_level, ch_type ):
    
# provides Gaussian-shaped DOS  
# in units of states per energy intervall eV per unit volume in m
# provides vector to store the DOS
#
# input: DOS_admin .. DOS administration data
#        N_G .. density of states
#        E   .. energy vector
#        E_level  .. energy of Gaussian center = mean value (in eV)
#        type .. Donor or acceptor?
#
# output: DOS vector with number of states per energy interval 
#         corresponding to requested Gaussian added
#         occupation vector (zero-valued)
#
# check whether donor or acceptor like


    if ch_type == 'N':      # donor
        label = 'DL'
    elif ch_type == 'P':    # acceptor
        label = 'AL'
    elif ch_type == '0':    # neutral
        label = '0L'


    DOS_temp = np.zeros_like(E)

    # add delta-shaped level
    interval_index = np.max(np.where(E<E_level)[0])
    DOS_temp[interval_index] = N_G

    DOS_admin = AddContributionToDOS(DOS_admin, DOS_temp, label,ch_type,
                                    E_level,0,N_G)

    return DOS_admin



###############################################################################
###############################################################################




def DensityOfBandStates(eff_mass,T):

# calculate the effective density of states for a given effective mass and
# temperature
#
# input:    eff_mass  .. effective mass (units of electron mass me)
#           T         .. temperature (K)
# output:   N         .. effective density of states for given temperature
#                     ..   in cubic meter (!)

#   considering the parameters
#
#   me = 9.11e-31 # kg
#   k = 1.38e-23  # SI
#   h = 6.626E-34 # SI
#
#   one arrives at a prefactor of
#
#   prefactor = 2*(2*pi*k*me)^(3/2)/h^3;

    prefactor = 4.8266e+21;

    return prefactor*(eff_mass*T)**(3/2)



###############################################################################
###############################################################################




def FIntegrationForBands(eta):
    
#    func_handle = @(E)integrandForBands(E, E_offset);

# upper boundary for integration is infinity
# approximate this boundary with a large value
#     y = integration(func_handle,0,1E4)
#
# use here approximate analytical expression taken from 
# https://www.eecis.udel.edu/~kolodzey/courses/ELEG667F06/667F06HW/Pierret_6_FermiDirac.pdf
# Table 4.2 (p. 118)

    if eta >1.3:
        a = (np.power(np.abs(eta-2.13),2.4) + 9.6)
        a = np.power(a,5/12)+eta+2.13
        a = np.power(a,-1.5)*3*np.sqrt(np.pi/2) + np.exp(-eta)
        a = 1/a
    
    else: 
        a = np.exp(-eta) + 2.7 # following wikipedia
        a = 1/a
          
# KZ   
# https://arxiv.org/pdf/0811.0116.pdf
# https://nanohub.org/resources/5475
# References
# [1]D. Bednarczyk and J. Bednarczyk, Phys. Lett. A, 64, 409 (1978)
# [2]J. S. Blakemore, Solid-St. Electron, 25, 1067 (1982)      
    
# Model proposed in [1]
# Expressions from eqs. (22)-(24) of [2]

        mu = np.power(eta,4) + 50 + 33.6 * eta * ( 1 - 0.68 * np.exp( -0.17 * np.power( eta + 1,2 ) ) );
        xi = 3 * np.sqrt( np.pi ) / ( 4 * np.power(mu,3/8) );
        a = np.power( np.exp( - eta ) + xi,-1 );

    return a

###############################################################################
###############################################################################



# Here the charge density within the conduction or valence band is 
# computed with the approximated Fermi integral 


def GetDensityInBand(chemical_potential,E_edge, eff_mass,T):

    if T==0 :
        T = 1E-3
        
    k_eV = 1.3806504E-23 / 1.602176487E-19  # Boltzmann constant in eV/K
    arg = -(np.abs(chemical_potential-E_edge)/k_eV)/T

    charge_density = 2/np.sqrt(np.pi)*DensityOfBandStates(eff_mass,T)*FIntegrationForBands(arg)
    
    #print(charge_density)
    
    return charge_density




###############################################################################
###############################################################################




def chargeNeutralityIntrinsic(Energy, E_C, E_V, m_n_eff,m_p_eff,T):
#   this function emulates the charge neutrality condition
#   desired shape F(E) with F(E) = 0 when condition holds
#   -> one evaluates the expression F(E) = n(E)-p(E)

    n = GetDensityInBand(Energy,E_C, m_n_eff,T)

    p = GetDensityInBand(Energy,E_V, m_p_eff,T)

# rather than checking for the absolute deviation from zero
# look into relative error 
    
    if (n+p) != 0:
        eval = (n-p)/(n+p)
    else:
        eval = 0
        
    return eval



###############################################################################
###############################################################################




def FindRootNestedIntervals( FuncEval,  Interval, initial_Guess, 
                            iter_threshold, max_iter_steps ):
    
# Utilization of method of nested intervals to find the root of a function
# defined on a onedimensional vector
#
# usage : 
#   chemical_potential = 
#   FindRootNestedIntervals(@(E) chargeNeutrality(E,E_C,E_V,m_n_eff,m_p_eff,T),...
#                           energies, (E_C-E_V)/2, 1d-3, 20);
#
#   or
#
#   fh = @(E) chargeNeutralityIntrinsic(E,E_C,E_V,m_n_eff,m_p_eff,T);
#   chemical_potential = FindRootNestedIntervals(fh,energies, ...
#                        (E_C-E_V)/2, 1d-3, 20); 
#
# requires: definition of function FEval
#
#           e.g., FuncEval(x) = chargeNeutralityIntrinsic(energies) 
#           checking for charge neutrality
# input:   
#     Interval        ..  vector containing search interval 
#                         {x} = [xmin,xmax] of target quantity x_target 
#     FuncEval        ..  string quantity f(x) dependent on x
#     initial_guess   ..  initial guess for x_target
#     iter_threshold  ..  threshold how close to f(x)==0 will be iterated 
#                         (until |f(x)| < iter_threshold)
#     max_iter_steps  ..  upper limit of number of interval splitting steps 
#
# output:
#     x_target        ..  one value x_target within interval {x} satisfying 
#                         the condition f(x)==0
#                         x_target is NOT necessarily member of the vector x


    x_max = np.max(Interval)
    x_min = np.min(Interval)

    x_target = 0.

    Iter =  1
    error = 1
    
    while ((np.abs(error) > iter_threshold) and (Iter < max_iter_steps )):
        
        if Iter == 1:
            x_center = initial_Guess
        else:
           x_center = (x_max - x_min)/2 + x_min


        # error is here an absolute error :(
        error = FuncEval(x_center) ;  

        if error > 0:
             x_max = x_center
        else:
             x_min = x_center

        Iter = Iter + 1

    x_target = x_center
    num_iter = Iter
    
    return x_target, num_iter, error
    


###############################################################################
###############################################################################




def FermiDirac(E, CP, T):
    
# Gives occupation probability with FERMI DIRAC distribution
# input: E   .. energy
#        CP .. chemical potential
#        T   .. temperature
# output: prob .. occupation probability

    k = 1.3806504e-23 # in SI
    q = 1.602176565e-19 

# convert eV to SI vy using a temperature equivalent
    T = T/q

    return 1/(np.exp((E-CP)/k/T)+1)



###############################################################################
###############################################################################



def GetDensityInLevel(chemical_potential, DOS_level_E_ref,DOS_level_N, T):
    
# provides electron density in a delta-shaped DOS as a function of the
# chemical potential and temperature
#
# here, there is no discrimination between donor/acceptor/neutral like
# states
#
# input: DOS_level_E_ref    .. energy of level / eV
#        DOS_level_N        .. density of level / eV
#        chemical_potential .. chemical potential / eV        
#        T                  .. temperature / K
# output: determined density

    density = 0

# check if temperature is passed as array or as single number
    if isinstance(T, np.ndarray):
        [ te if (te !=0) else 1E-3 for te in T]     
    else:
        if T==0 :
            T = 1E-3
        
           
#   occupation = f_FD(EL, chem_pot,T) * N_L
    
    density = FermiDirac(DOS_level_E_ref, chemical_potential,T) * DOS_level_N
    # if chemical_potential or temperature were provided as 1D array, 
    # return a 1D aray
    # else a float number
    if not(isinstance(T, np.ndarray)) or (isinstance(chemical_potential, np.ndarray)):
        density = density.item()

    return density



###############################################################################
###############################################################################




def GaussDOS(E, Ecenter, Ewidth):
    
# gives DOS at a specific energy for a GAUSSIAN-shaped DOS
#
# input: E .. energy
#        Ecenter .. shape-defining parameters, here mean value
#        Ewidth  .. shape-defining parameter, here width sigma

    return 1/np.sqrt(2*np.pi)/Ewidth*np.exp(-((E-Ecenter)/Ewidth)**2/2 )



###############################################################################
###############################################################################



def FDIntegrantGauss(E, DOS_Gauss_N,DOS_Gauss_E_ref,DOS_Gauss_Param,chemical_potential,T):
    
# provides the integral kernel of the Fermi Dirac Integral with respect to
# a Gaussian-shaped DOS
#
# note that also 
#         - non-delta (i.e., sharp levels) or 
#         - sqrt-shaped bands 
# can be dealt with this evaluation
#
# input :  E                  .. energy / eV
#          DOS_Gauss_item     .. one data set from DOS admin describing either
#                                a Gauss- or differently shaped DOS (excluding
#                                delta and sqrt-shaped contributions)
#          chemical_potential .. chemical potential / eV
#          T                  .. temperature / K

    int_kernel = 0
    int_kernel = FermiDirac(E,chemical_potential,T)*DOS_Gauss_N*GaussDOS(E,DOS_Gauss_E_ref,DOS_Gauss_Param)
    
    return int_kernel


###############################################################################
###############################################################################





def GetDensityInGauss(chemical_potential,DOS_Gauss_Label,DOS_Gauss_N,DOS_Gauss_E_ref,DOS_Gauss_Param,T):
    
# provides electron density in a Gaussian-shaped DOS as a function of the
# chemical potential and temperature
#
# here, there is no discrimination between donor/acceptor/neutral like
# states
#
# here, the integration is done independent of the energy mesh
#
# input: DOS_Gauss_item     .. data set associated to Gaussian DOS
#        chemical_potential .. chemical potential / eV        
#        T                  .. temperature / K
# output: determined density

    density = 0;

    if T==0 :
        T = 1E-3
# check that neither sqrt- or delta-shaped items are considered

    if ('L' not in DOS_Gauss_Label) and ('B' not in DOS_Gauss_Label):

        func_handle = lambda E: FDIntegrantGauss(E,DOS_Gauss_N,DOS_Gauss_E_ref,DOS_Gauss_Param,chemical_potential,T)
        #density = sci.integrate.quad(func_handle,-np.inf,np.inf)[0]
        density = quad(func_handle,-np.inf,np.inf)[0]
         
    return density



###############################################################################
###############################################################################




#   this function emulates the charge neutrality condition
#   desired shape F(E) with F(E) = 0 when condition holds
#   -> one evaluates the expression F(E) = n(E)-p(E)
#
# input:    Energy    .. reference energy at which condition ought to be
#                        evaluated  / eV
#           DOS_admin .. complete information on DOS
#           m_n_eff   .. effective electron mass
#           m_p_eff   .. effective hole mass
#           T         .. temperature / K
#    
# note:     dT_smear  .. this hidden parameter helps to prevent numerical artefacts
#                        T_smear >> 0K artificially smears out the FD distribution 
#                        at low temperature; for the true FD distribution, the 
#                        integration with respect to the 
#                        delta-shaped DOS gives an ill-defined value          
          

def chargeNeutrality(Energy, DOS_admin, m_n_eff,m_p_eff,T):
    
    
    dT_smear = 350 # in K

    N = 0      # total density associated to negative charges / inverse cubic meters
    P = 0      # total density associated to posiitve charges / inverse cubic meters

    n = 0 
    p = 0 # temporary densities

# here we need to collect all contributions relevant for charge neutrality

    number_of_DOS_entries = np.size(DOS_admin.Label)

# check for type of DOS contribution

    for k in range(number_of_DOS_entries):
        
#             'B'  band state (sqrt-shaped)
        
        if (DOS_admin.Label[k][0] == 'CB'):
            N = GetDensityInBand(Energy,DOS_admin.E_ref[k], m_n_eff,T)
        elif (DOS_admin.Label[k][0] == 'VB'):
            P = GetDensityInBand(Energy,DOS_admin.E_ref[k], m_p_eff,T)
        
        

#             'L'  sharp level (delta-shaped)

        elif ('L' in DOS_admin.Label[k][0]) and (DOS_admin.Type[k][0] == 'P'):
            n = GetDensityInLevel(Energy,DOS_admin.E_ref[k],DOS_admin.N[k],T+dT_smear)
            N = N + n
        elif ('L' in DOS_admin.Label[k][0]) and (DOS_admin.Type[k][0] == 'N'):
            p = DOS_admin.N[k]-GetDensityInLevel(Energy,DOS_admin.E_ref[k],DOS_admin.N[k],T+dT_smear)
            P = P + p

              
#            'G'  Gaussian-shaped

        elif ('G' in DOS_admin.Label[k][0]) and (DOS_admin.Type[k][0] == 'P'):
            n = GetDensityInGauss(Energy,DOS_admin.Label[k][0],DOS_admin.N[k],
                      DOS_admin.E_ref[k],DOS_admin.Param[k],T)
            N = N + n
        elif ('G' in DOS_admin.Label[k][0]) and (DOS_admin.Type[k][0] == 'N'):
            p = DOS_admin.N[k] - GetDensityInGauss(Energy,DOS_admin.Label[k][0],DOS_admin.N[k],
                      DOS_admin.E_ref[k],DOS_admin.Param[k],T)
            P = P + p


# rather than checking for the absolute deviation from zero
# look into relative error 

    if (N+P != 0):
        Eval = (N-P)/(N+P)
    else:
        Eval = 0.0

    return Eval


###############################################################################################################
###############################################################################################################





def GetFullDOS(DOS_admin):
    
# this function provides the full density of states in a vector 
# consistent with the energy interval
# input:    DOS_admin .. data containing DOS information
# output:   DOS       .. vector containing DOS

    DOS = np.zeros_like(DOS_admin.Energies[0])
    max_index = np.size(DOS_admin.Label)

    for k in range(max_index):
        DOS = DOS + DOS_admin.Energies[k]
        
    return DOS



#########################################################################################################
#########################################################################################################



def AddGaussToDOS( DOS_admin,E, N_G, E1, E2, ch_type ):
    
# provides Gaussian-shaped DOS  
# in units of states per energy intervall eV per unit volume in m
# provides vector to store the DOS
#
# input: DOS_admin .. DOS administration data
#        N_G .. density of states
#        E   .. energy vector
#        E1  .. energy of Gaussian center = mean value (in eV)
#        E2  .. width of Gaussian = "sigma"            (in eV)
#        type .. Donor or acceptor?
#
# output: DOS vector with number of states per energy interval 
#         corresponding to requested Gaussian added
#         occupation vector (zero-valued)
#
# check whether donor or acceptor like

    if ch_type == 'N':      # donor
        label = 'DG'
    elif ch_type == 'P':    # acceptor
        label = 'AG'
    elif ch_type == '0':    # neutral
        label = '0G'

    DOS_temp = np.zeros_like(E)
    DOS_temp = N_G*1/np.sqrt(2*np.pi)/E2*np.exp(-((E-E1)/E2)**2/2 )

    DOS_admin = AddContributionToDOS(DOS_admin, DOS_temp, label,ch_type,E1,E2,N_G)
    
    return DOS_admin

#%occ_vector = zeros(size(E),'like',E)

