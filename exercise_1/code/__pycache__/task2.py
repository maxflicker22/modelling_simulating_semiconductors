#%%
import numpy as np
import matplotlib.pyplot as plt
import semiconductor_functions as sf

np.seterr('ignore') # Ignores Over/Underflow warnings/errors

k_eV = 1.3806504E-23 / 1.602176487E-19 # Boltzmann constant in eV/K
T = 300 # Temperature in K
resolution = 1500 # Number of points in the plot

# threshold and maximal number of iterations for root finding
tolerance = 1e-4
max_RF_iter = 30

# doping
dopmin = 1E15 # cm^-3
dopmax = 1E25 # cm^-3
N_steps = 250 # Number of steps in the doping density range

E_doffset = 0.15 # eV
doping_type = 'acceptor' 

doping_densities = np.power(10, np.linspace(np.log10(dopmin),np.log10(dopmax),N_steps))

# Storing vectors
n = np.zeros_like(doping_densities)
n_i = np.zeros_like(doping_densities)
p = np.zeros_like(doping_densities)
trapden = np.zeros_like(doping_densities)

Ndop_ionized = np.zeros_like(doping_densities)
trap_ionized = np.zeros_like(doping_densities)

chemical_potential = np.zeros_like(doping_densities)
chemical_potential_i = np.zeros_like(doping_densities)

EC, EV, m_n_eff, m_p_eff = sf.AssignSemiconductor('GaAs')
E_dopant = EV + E_doffset
dopant_label = 'N$_A$'
dopant_ion = '$^-$'
dopant_marker = 'P'

trap_dens = [0, 1E22, 1E22]
E_toffset = [0, 0.1, -0.1]
sigma = 0.05

E_min = EV - 0.5    # eV
E_max = EC + 0.5    # eV
E_guess = (EC+EV)/2. + 0.2  # Guessing value for the FindRootNestedIntervals routine

for experiment in zip(trap_dens, E_toffset):
    trap_dens = experiment[0]
    E_toffset = experiment[1]
    for k, dopant_density in enumerate(doping_densities):
        
        energies, DOS, occupation = sf.InitializeEnergyAndDOS(E_min, E_max, EV, EC, resolution)
        DOS_Admin = sf.InitializeDOSAdministration(energies)
        DOS_Admin = sf.AddConductionBandToDOS(DOS_Admin, energies, EC, m_n_eff)
        DOS_Admin = sf.AddValenceBandToDOS(DOS_Admin,energies,EV,m_p_eff)
        DOS_Admin = sf.AddLevelToDOS(DOS_Admin, energies, dopant_density, E_dopant, dopant_marker)

        if trap_dens > 0:
            E_trap = (EC-EV)/2 + E_toffset
            DOS_Admin = sf.AddGaussToDOS(DOS_Admin, energies, trap_dens, E_trap, sigma, 'P')

        fh = lambda E: sf.chargeNeutralityIntrinsic(E, EC, EV, m_n_eff, m_p_eff, T)
        chemical_potential_i[k], _, _ = sf.FindRootNestedIntervals(fh, energies, E_guess, tolerance, max_RF_iter)

        n_t = sf.GetDensityInBand(chemical_potential_i[k],EC,m_n_eff, T)
        p_t = sf.GetDensityInBand(chemical_potential_i[k],EV,m_p_eff, T)
        n_i[k] = np.sqrt(n_t*p_t) 

        fh2 = lambda E: sf.chargeNeutrality(E, DOS_Admin, m_n_eff, m_p_eff, T)
        chemical_potential[k], _, _ = sf.FindRootNestedIntervals(fh2, energies, chemical_potential_i[k], tolerance, max_RF_iter)

        n[k] = sf.GetDensityInBand(chemical_potential[k],EC,m_n_eff, T)
        p[k] = sf.GetDensityInBand(chemical_potential[k],EV,m_p_eff, T)

        if trap_dens > 0:
            trapden[k] = sf.GetDensityInGauss(chemical_potential[k], DOS_Admin.Label[3], DOS_Admin.N[3], DOS_Admin.E_ref[3], DOS_Admin.Param[3], T)

        Ndop_ionized[k] = sf.GetDensityInLevel(chemical_potential[k], DOS_Admin.E_ref[2],DOS_Admin.N[2],T)/DOS_Admin.N[2]
        
        if trap_dens> 0:
            trap_ionized[k] = sf.GetDensityInGauss(chemical_potential[k],
                            DOS_Admin.Label[3], 
                            DOS_Admin.N[3],
                            DOS_Admin.E_ref[3],
                            DOS_Admin.Param[3],
                            T)/DOS_Admin.N[3]
    
    col = 'blue'

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,6), sharex=True)

    # Plot 1: Chemical potential
    ax1.semilogx(doping_densities, chemical_potential_i, label='Intrinsic Chemical Potential', linestyle='-.', color='k')
    ax1.semilogx(doping_densities, chemical_potential, label='Chemical Potential', color='k')
    ax1.semilogx(doping_densities, EV*np.ones_like(doping_densities), label='Valence / Conduction Band', color='grey')
    ax1.semilogx(doping_densities, EC*np.ones_like(doping_densities), color='grey')

    ax1.set_xlim(dopmin, dopmax)
    ax1.set_xscale('log')
    ax1.set_xlabel('Doping Density (m$^{-3}$)')
    ax1.set_ylabel('Chemical Potential (eV)')
    ax1.legend(loc='upper right')
    ax1.grid(False)

    # Plot 2: Electron and hole densities
    ax2.plot(doping_densities, n_i, label=r'n$_i$ (Intrinsic)', linestyle='--', color='blue')
    ax2.plot(doping_densities, n, label=r'$n$ (Electrons)', color='blue')
    ax2.plot(doping_densities, p, label=r'$p$ (Holes)', color='red')
    
    ax2.set_xlim(dopmin, dopmax)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Doping Density (m$^{-3}$)')
    ax2.set_ylabel('Carrier Density (m$^{-3}$)')
    ax2.legend(loc='center right')
    ax2.grid(False)

    fig2, ax3 = plt.subplots()
    ax3.plot(doping_densities, Ndop_ionized, label=r'$N_A^-/N_A$ Ionized dopants', color='red')
    if trap_dens > 0:
        ax3.loglog(doping_densities, trap_ionized, label=r'$N_T^-/N_T$ Ionized traps', color='green')
    ax3.set_xlim(dopmin, dopmax)
    ax3.set_xscale('log')
    #ax2.set_yscale('log')
    ax3.set_xlabel('Doping Density (m$^{-3}$)')
    ax3.set_ylabel(r'$N_i^-/N_i$')
    ax3.legend(loc='center left')
    ax3.grid(False)

    plt.tight_layout()
    plt.show()

    #%%