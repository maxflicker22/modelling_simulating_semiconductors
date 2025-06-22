
function [current_den,profile, it_out, res_out] = current4voltage(voltage,...
                                                                  temperature, ...
                                                                  device,...
                                                                  itercontrol)

% obtains current density of a device for a given applied voltage
%
% INPUT: 
%          voltage     .. applied external voltage
%          device      .. data record providing all information on mesh, doping,
%                         and material
%          itercontrol .. data record providing parameters to control
%                         Gummel iterations
%
% OUTPUT:  current_den .. current density
%          profile     .. data record containing 1D profiles of
%          profile.n   .. electron density
%          profile.h   .. hole density
%          profile.Jn  .. electron current density
%          profile.Jp  .. hole current density
%          profile.psi .. electrostatic potential
%          profile.efield .. electric field
%          profile.EFn .. quasi Fermi level for electrons
%          profile.EFp .. quasi Fermi level for holes
%          profile.Ec  .. energy of conduction band minimum incl. band bending
%          profile.Ev  .. energy of valence band maxmum incl. band bending
%          profile.mobilitye .. electron mobility
%          profile.mobilityp .. hole mobility

secs1d_physical_constants;

% initialize output profiles
% necessary ??
profile_elements = device.mesh.Nelements;


profile.n   = zeros(profile_elements,1);
profile.p   = zeros(profile_elements,1);
profile.Jn  = zeros(profile_elements,1);
profile.Jp  = zeros(profile_elements,1);
profile.EFn = zeros(profile_elements,1);
profile.EFp = zeros(profile_elements,1);
profile.Ec  = zeros(profile_elements,1);
profile.Ev  = zeros(profile_elements,1);
profile.psi = zeros(profile_elements,1);
profile.efield = zeros(profile_elements,1);

% initialize profiles for internal handling
% necessary ??
nin = zeros(profile_elements,1);
nout = zeros(profile_elements,1);


%%------ INITIALIZE SIMULATED QUANTITIES AND MESH -------------------------------------------------------

% useful constants
Vth 	 = Kb * temperature / q;

% set geometry
L = device.geometry.length;
L_p_layer = device.geometry.p_layer_length;


% -- applied voltage to be used as boundary condition --
% it is also possible to conceive a symmetric voltage distribution

V_p = voltage;
V_n = 0;

% -- Distributions across device cross section --

% dielectric constant (silicon)
er = device.material.esir * ones (device.mesh.Nelements, 1);

% here it might be useful to reset ni in accord with the given temperature
% by calling
device.material = silicon_material_properties(temperature);
% store intrinsic carrier density
ni = device.material.ni;

% for convenience: remember device geometry
mesh_left = device.mesh.x <= device.mesh.xm;
mesh_right = device.mesh.x > device.mesh.xm;

% set doping profile [m^{-3}]
Na = device.doping.NA * (mesh_left); % set doping densities
Nd = device.doping.ND * (mesh_right);
% avoid zero doping
D = Nd - Na;

% set parameters for iteration control
% set control parameters for simulation flow
% tolerances for convergence checks
toll = itercontrol.tol;
maxit = itercontrol.maxit;
ptoll = itercontrol.ptol;
pmaxit = itercontrol.pmaxit;

% initial guess for n, p, V, phin, phip

% quasi Fermi levels
Fp = V_p * (mesh_left);
Fn = Fp;

% charge carrier densities
p = abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni./abs(D)) .^2)) .* (mesh_left) + ...
ni^2 ./ (abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^2))) .* (mesh_right);

n = abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^ 2)) .* (mesh_right) + ...
ni ^ 2 ./ (abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^2))) .* (mesh_left);

% electrostatic potential
V = Fn + Vth * log (n / ni);

% ----------- SCALE QUANTITIES ---

% provide all necessary scaling factors
xbar = device.geometry.length;      % [m]
nbar = norm(D, 'inf');              % [m^{-3}]
Vbar = Vth;                         % [V]
mubar = max (device.material.u0n, device.material.u0p); % [m^2 V^{-1} s^{-1}]
tbar = xbar^2 / (mubar * Vbar);     % [s]
Rbar = nbar / tbar;                 % [m^{-3} s^{-1}]
Ebar = Vbar / xbar;                 % [V m^{-1}]
Jbar = q * mubar * nbar * Ebar;     % [A m^{-2}]
CAubar = Rbar / nbar^3;             % [m^6 s^{-1}]
abar = 1/xbar;                      % [m^{-1}]


% scaling parameters and quantities

l2 = e0 * Vbar / (q * nbar * xbar^2);
theta = ni / nbar;

xin = device.mesh.x / xbar;
Din = D / nbar;
Nain = Na / nbar;
Ndin = Nd / nbar;
pin = p / nbar;
nin = n / nbar;
Vin = V / Vbar;
Fnin = Vin - log (nin);
Fpin = Vin + log (pin);

tnin = device.material.tn / tbar;
tpin = device.material.tp / tbar;

u0nin = device.material.u0n / mubar;
uminnin = device.material.uminn / mubar;
vsatnin = device.material.vsatn / (mubar * Ebar);

u0pin = device.material.u0p / mubar;
uminpin = device.material.uminp / mubar;
vsatpin = device.material.vsatp / (mubar * Ebar);

Nrefnin = device.material.Nrefn / nbar;
Nrefpin = device.material.Nrefp / nbar;

Cnin = device.material.Cn / CAubar;
Cpin = device.material.Cp / CAubar;

anin = device.material.an / abar;
apin = device.material.ap / abar;
Ecritnin = device.material.Ecritn / Ebar;
Ecritpin = device.material.Ecritp / Ebar;


% solve the problem using the full DD model
[nout, pout, Vout, Fnout, Fpout, Jnout, Jpout, mobnout, mobpout, it, res] = ...
    secs1d_dd_gummel_map (xin, Din, Nain, Ndin, pin, nin, Vin, Fnin, Fpin, ...
                          l2, er, ...
                          u0nin, uminnin, vsatnin, device.material.betan, Nrefnin, ...
                          u0pin, uminpin, vsatpin, device.material.betap, Nrefpin, ...
                          theta, ...
                          tnin, tpin, Cnin, Cpin, anin, apin, ...
                          Ecritnin, Ecritpin, ...
                          toll, maxit, ptoll, pmaxit);

it_out = it;                      
res_out = res;

% Descaling procedure + set ouptput profiles

% charge carrier densities
profile.n   = nout*nbar;
profile.p   = pout*nbar;
% current density profiles
profile.Jn  = Jnout * Jbar;
profile.Jp  = Jpout * Jbar;
% electrostatic potential
profile.psi = Vout*Vbar;
% quasi Fermi levels
profile.EFn = (-1)*(profile.psi - Vth*log(profile.n/device.material.ni));
profile.EFp = (-1)*(profile.psi + Vth*log(profile.p/device.material.ni));
% band structure
profile.Ec  = Vth*log(device.material.Nc./profile.n)+profile.EFn;
profile.Ev  = -Vth*log(device.material.Nv./profile.p)+profile.EFp;
% electrostatic field
dV = diff(profile.psi);
dx = diff(device.mesh.x);
profile.efield = -dV./dx;
% mobilities
profile.mobilitye = mobnout * mubar;
profile.mobilityp = mobpout * mubar;

% total current density
Jtot = profile.Jn + profile.Jp;
current_den = Jtot(1);


end % function
