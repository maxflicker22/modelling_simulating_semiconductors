%% 2019: Karin Zojer
%% ADAPTED FROM SECS1D_silicon_material_properties
%% Copyright (C) 2004-2012  Carlo de Falco
%%
%% This file is part of 
%% SECS1D - A 1-D Drift--Diffusion Semiconductor Device Simulator
%%
%% SECS1D is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
%% (at your option) any later version.
%%
%% SECS1D is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with SECS1D; If not, see <http://www.gnu.org/licenses/>.

%% material properties for silicon and silicon dioxide
%%
%% esir       = relative electric permittivity of silicon
%% esio2r     = relative electric permittivity of silicon dioxide
%% esi 	      = electric permittivity of silicon
%% esio2      = electric permittivity of silicon dioxide
%% mn         = effective mass of electrons in silicon
%% mh         = effective mass of holes in silicon
%% 
%% u0n        = low field electron mobility
%% u0p        = low field hole mobility
%% uminn      = parameter for doping-dependent electron mobility
%% betan      = idem
%% Nrefn      = idem
%% uminp      = parameter for doping-dependent hole mobility
%% betap      = idem
%% Nrefp      = idem
%% vsatn      = electron saturation velocity
%% vsatp      = hole saturation velocity
%% tp         = electron lifetime
%% tn         = hole lifetime
%% Cn         = electron Auger coefficient
%% Cp         = hole Auger coefficient
%% an         = impact ionization rate for electrons
%% ap         = impact ionization rate for holes
%% Ecritn     = critical field for impact ionization of electrons
%% Ecritp     = critical field for impact ionization of holes 
%% Nc         = effective density of states in the conduction band
%% Nv         = effective density of states in the valence band
%% Egap       = bandgap in silicon
%% EgapSio2   = bandgap in silicon dioxide
%% 
%% ni         = intrinsic carrier density
%% Phims      = metal to semiconductor potential barrier

function [device_struc] = silicon_material_properties(T0)

    a = 2+1;
    secs1d_physical_constants;
    
    device_struc.esir 	      = 11.7;
    device_struc.esio2r 	  = 3.9;
    device_struc.esi 	      = e0 * device_struc.esir;
    device_struc.esio2 	      = e0 * device_struc.esio2r;
    device_struc.mn           = 0.26*mn0;
    device_struc.mh           = 0.18*mn0;

    device_struc.qsue         = q / device_struc.esi;

    device_struc.u0n          = 1417e-4;
    device_struc.u0p          = 480e-4;
    device_struc.uminn        = device_struc.u0n;            % ref. value: 65e-4;
    device_struc.uminp        = device_struc.u0p;            % ref. value: 47.7e-4;
    device_struc.betan        = 0.72;
    device_struc.betap        = 0.76;
    device_struc.Nrefn        = 8.5e22;
    device_struc.Nrefp        = 6.3e22;
    device_struc.vsatn        = inf;            % ref. value: 1.1e5;
    device_struc.vsatp        = inf;            % ref. value: 9.5e4;

    device_struc.tp           = inf;            % ref. value: 1e-6;
    device_struc.tn           = inf;            % ref. value: 1e-6;

    device_struc.Cn           = 0;              % ref. value: 2.8e-31*1e-12; 
    device_struc.Cp           = 0;              % ref. value: 9.9e-32*1e-12;   
    device_struc.an           = 0;              % ref. value: 7.03e7;
    device_struc.ap           = 0;              % ref. value: 6.71e7;
    device_struc.Ecritn       = 1.231e8; 
    device_struc.Ecritp       = 1.693e8;

    device_struc.mnl          = 0.98*mn0;
    device_struc.mnt          = 0.19*mn0;
    device_struc.mndos        = (device_struc.mnl*device_struc.mnt*...
                                 device_struc.mnt)^(1/3); 

    device_struc.mhh         = 0.49*mn0;
    device_struc.mlh         = 0.16*mn0;
    device_struc.mhdos       = (device_struc.mhh^(3/2)+device_struc.mlh^(3/2))^(2/3);

    device_struc.Nc          = (6/4)*(2*device_struc.mndos*Kb*T0/(hbar^2*pi))^(3/2);   
    device_struc.Nv          = (1/4)*(2*device_struc.mhdos*Kb*T0/(hbar^2*pi))^(3/2);
    device_struc.Eg0         = 1.16964*q;
    device_struc.alfaEg      = 4.73e-4*q;
    device_struc.betaEg      = 6.36e2;
    device_struc.Egap        = device_struc.Eg0-device_struc.alfaEg*((T0^2)/(T0+device_struc.betaEg));
    device_struc.Ei          = device_struc.Egap/2+Kb*T0/2*log(device_struc.Nv/device_struc.Nc);
    device_struc.EgapSio2    = 9*q;
    device_struc.deltaEcSio2 = 3.1*q;
    device_struc.deltaEvSio2 = device_struc.EgapSio2-device_struc.Egap-device_struc.deltaEcSio2;

    device_struc.ni          = sqrt(device_struc.Nc*device_struc.Nv)*exp(-device_struc.Egap/(2*(Kb * T0)));
    device_struc.Phims       = - device_struc.Egap /(2*q);

end