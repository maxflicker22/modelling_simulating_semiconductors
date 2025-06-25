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

%%
%% Solve the scaled stationary bipolar DD equation system using Gummel algorithm.
%%
%% [n, p, V, Fn, Fp, Jn, Jp, mobn, mobp, it, res] = secs1d_dd_gummel_map (x, D, Na, Nd, 
%%                                                       pin, nin, Vin, Fnin, 
%%                                                       Fpin, l2, er, u0n, 
%%                                                       uminn, vsatn, betan, 
%%                                                       Nrefn, u0p, uminp, vsatp, 
%%                                                       betap, Nrefp, theta, tn, tp, 
%%                                                       Cn, Cp, an, ap, Ecritnin, Ecritpin, 
%%                                                       toll, maxit, ptoll, pmaxit)         
%%
%%     input: 
%%            x                        spatial grid
%%            D, Na, Nd                doping profile
%%            pin                      initial guess for hole concentration
%%            nin                      initial guess for electron concentration
%%            Vin                      initial guess for electrostatic potential
%%            Fnin                     initial guess for electron Fermi potential
%%            Fpin                     initial guess for hole Fermi potential
%%            l2                       scaled Debye length squared
%%            er                       relative electric permittivity
%%            u0n, uminn, vsatn, Nrefn electron mobility model coefficients
%%            u0p, uminp, vsatp, Nrefp hole mobility model coefficients
%%            theta                    intrinsic carrier density
%%            tn, tp, Cn, Cp, 
%%            an, ap, 
%%            Ecritnin, Ecritpin       generation recombination model parameters
%%            toll                     tolerance for Gummel iterarion convergence test
%%            maxit                    maximum number of Gummel iterarions
%%            ptoll                    convergence test tolerance for the non linear
%%                                     Poisson solver
%%            pmaxit                   maximum number of Newton iterarions
%%
%%     output: 
%%             n     electron concentration
%%             p     hole concentration
%%             V     electrostatic potential
%%             Fn    electron Fermi potential
%%             Fp    hole Fermi potential
%%             Jn    electron current density
%%             Jp    hole current density
%%             mobn  electron mobility
%%             mobp  hole mobility
%%             it    number of Gummel iterations performed
%%             res   total potential increment at each step

function [n, p, V, Fn, Fp, Jn, Jp, mobn, mobp, it, res] = secs1d_dd_gummel_map (x, D, Na, Nd, pin, nin, Vin, ...
                                                                    Fnin, Fpin, l2, er, u0n, uminn, ...
                                                                    vsatn,betan,Nrefn, u0p, ...
                                                                    uminp,vsatp,betap,Nrefp, ...
                                                                    theta, tn, tp, Cn, Cp, an, ap, ... 
		                                                    Ecritnin, Ecritpin, toll, maxit, ...
                                                                    ptoll, pmaxit)         

  
 % this function performs Gummel iterations
 % it tracks the progress of the iterations for
 % n,p, electrostatic potential (= V) and the quasi-Fermi levels
 % in the following form:
 % Each of the above quantities will be kepts as a three-column "matrix"
 %  col(1) ... e.g., p(:,1) ... initial guess
 %  col(2) ... e.g., p(:,2) ... updated potential or updated charge density 
 %                              due to updated potential
 %                              e.g., p  = exp (-V(:,2) + Fp);
 %  col(3) ... e.g., p(:,3) ... p due to charge displacement according to
 %                              Gummel discretization and given electrostatic potential
 
% load initial values for the iteration                                                              
  p  = pin;
  n  = nin;
  V  = Vin;
  Fp = Fpin;
  Fn = Fnin;
  % compute the individual separations between adjacent nodes
  % mesh cell = node
  % check whether this implementation allows to incorporate inhomogenous grids
  dx  = diff (x);
  % twice the mean node separation
  dxm = (dx(1:end-1) + dx(2:end));

  Nnodes = numel(x);
  Nelements = Nnodes -1;

  % initialize current density profiles
  Jn = zeros (Nelements, 1);
  Jp = zeros (Nelements, 1);

  for it = 1:maxit
      
    %-----------------
    % UPDATE POTENTIAL
    %------------------  
    % solve Poisson equation for actual charge density
    % update also n and p
    [V(:,2), n(:,2), p(:,2)] = secs1d_nlpoisson_newton (x, [1:Nnodes], V(:,1), n(:, 1), ...
                                                        p(:,1), Fn(:,1), Fp(:,1), D, l2, ...
                                                        er, ptoll, pmaxit); 
    % compute electric field
    dV = diff (V(:, 2));
    E  = - dV ./ dx;
    % compute Bernoulli function for the two arguments
    % Bp = B(psi_i+1-psi_i) and 
    % Bn = B(psi_i-psi_i-1)
    [Bp, Bm] = bimu_bernoulli (dV);
    
    %-----------------
    % UPDATE ELECTRONS
    %-----------------
    % build matrix for electron continuity equation
    
    % compute recombination and generation rates
    % here: Generation counteracts recombination [check]
    [Rn, Rp, Gn, Gp, II] = generation_recombination_model (x, n(:, end), p(:, end),...
	                                                   E, Jn, Jp, tn, tp, theta, ...
                                                           Cn, Cp, an, ap, Ecritnin, Ecritpin); 
    
    % compute field-dependent electron mobility
    mobility = mobility_model (x, Na, Nd, Nrefn, E, u0n, uminn, vsatn, betan);
    
 
    % contributions from divergence of  electron current density
    A = bim1a_advection_diffusion (x, mobility, 1, 1, V(:, 2));
    % contributions due to recombination
    % "reaction" : terms reducing or enhancing n via interactions
    %              with other species (here only p)
    M = bim1a_reaction (x, 1, Rn) + bim1a_reaction (x, II, 1);
    % contributions due to generation
    R = bim1a_rhs (x, 1, Gn);
  
    A = A + M;
  
    % store initial values
    n(:,3) = nin;
    % solve linear system of equations for electron density
    n(2:end-1,3) = A(2:end-1, 2:end-1) \ (R(2:end-1) - A(2:end-1, [1 end]) * nin ([1 end]));
   
    % update (scaled) electron quasi Fermi level
    Fn(:,2) = V(:,2) - log (n(:, 3));
    % update current density
    Jn =  mobility .* (n(2:end, 2) .* Bp - n(1:end-1, 2) .* Bm) ./ dx; 

    % update recombination and generation (current density dependent
    % contributions)
    [Rn, Rp, Gn, Gp, II] = generation_recombination_model (x, n(:, end), p(:, end), ...
	                                                   E, Jn, Jp, tn, tp, theta, ...
                                                           Cn, Cp, an, ap, Ecritnin, Ecritpin);

    % compute field-dependent hole mobility                                                       % compute field-dependent electron mobility
    mobility = mobility_model (x, Na, Nd, Nrefp, E, u0p, uminp, vsatp, betap);

    %-----------------
    % UPDATE HOLES
    %-----------------
    % build matrix for hole continuity equation
    
    % contributions from divergence of current density
    A = bim1a_advection_diffusion (x, mobility, 1, 1, -V(:, 2));
    % contributions due to recombination
    % "reaction" : terms reducing or enhancing p via interactions
    %              with other species (here only n)
    M = bim1a_reaction (x, 1, Rp) + bim1a_reaction (x, II, 1);
    % contributions due to generation
    R = bim1a_rhs (x, 1, Gp);
    
    A = A + M;
  
    % store initial values
    p(:,3) = pin;
    % solve linear system of equations for hole density
    p(2:end-1,3) = A(2:end-1, 2:end-1) \ (R(2:end-1) - A(2:end-1, [1 end]) * pin ([1 end]));
    
    % update (scaled) hole quasi Fermi level
    Fp(:,2) = V(:,2) + log (p(:,3));
    % update hole current density   
    Jp = -mobility .* (p(2:end, 2) .* Bm - p(1:end-1, 2) .* Bp) ./ dx;

    % determine deviation of potentials from previous iteration
    nrfn   = norm (Fn(:,2) - Fn(:,1), inf);
    nrfp   = norm (Fp(:,2) - Fp(:,1), inf);
    nrv    = norm (V(:,2)  - V(:,1),  inf);
    res(it) = max  ([nrfn; nrfp; nrv]);

    if (res(it) < toll)
      break
    end
    
    % set new initial values
    V(:,1)  = V(:,end);
    p(:,1)  = p(:,end) ;
    n(:,1)  = n(:,end);
    Fn(:,1) = Fn(:,end);
    Fp(:,1) = Fp(:,end);  
    
  end
 
  % after convergence is achieved: pass updated values
  n  = n(:,end);
  p  = p(:,end);
  V  = V(:,end);
  Fn = Fn(:,end);
  Fp = Fp(:,end);  

  mobn = mobility_model (x, Na, Nd, Nrefn, E, u0n, uminn, vsatn, betan);
  mobp = mobility_model (x, Na, Nd, Nrefp, E, u0p, uminp, vsatp, betap);
end



%!demo
%! % physical constants and parameters
%! secs1d_physical_constants;
%! secs1d_silicon_material_properties;
%! 
%! % geometry
%! L  = 10e-6;          % [m] 
%! xm = L/2;
%! 
%! Nelements = 1000;
%! x         = linspace (0, L, Nelements+1)';
%! sinodes   = [1:length(x)];
%! 
%! % dielectric constant (silicon)
%! er = esir * ones (Nelements, 1);
%! 
%! % doping profile [m^{-3}]
%! Na = 1e23 * (x <= xm);
%! Nd = 1e23 * (x > xm);
%! 
%! % avoid zero doping
%! D  = Nd - Na;  
%!  
%! % initial guess for n, p, V, phin, phip
%! V_p = -1;
%! V_n =  0;
%! 
%! Fp = V_p * (x <= xm);
%! Fn = Fp;
%! 
%! p = abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni./abs(D)) .^2)) .* (x <= xm) + ...
%!     ni^2 ./ (abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^2))) .* (x > xm);
%! 
%! n = abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^ 2)) .* (x > xm) + ...
%!     ni ^ 2 ./ (abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^2))) .* (x <= xm);
%! 
%! V = Fn + Vth * log (n / ni);
%! 
%! % scaling factors
%! xbar = L;                       % [m]
%! nbar = norm(D, 'inf');          % [m^{-3}]
%! Vbar = Vth;                     % [V]
%! mubar = max (u0n, u0p);         % [m^2 V^{-1} s^{-1}]
%! tbar = xbar^2 / (mubar * Vbar); % [s]
%! Rbar = nbar / tbar;             % [m^{-3} s^{-1}]
%! Ebar = Vbar / xbar;             % [V m^{-1}]
%! Jbar = q * mubar * nbar * Ebar; % [A m^{-2}]
%! CAubar = Rbar / nbar^3;         % [m^6 s^{-1}]
%! abar = 1/xbar;                  % [m^{-1}]
%! 
%! % scaling procedure
%! l2 = e0 * Vbar / (q * nbar * xbar^2);
%! theta = ni / nbar;
%! 
%! xin = x / xbar;
%! Din = D / nbar;
%! Nain = Na / nbar;
%! Ndin = Nd / nbar;
%! pin = p / nbar;
%! nin = n / nbar;
%! Vin = V / Vbar;
%! Fnin = Vin - log (nin);
%! Fpin = Vin + log (pin);
%! 
%! tnin = tn / tbar;
%! tpin = tp / tbar;
%! 
%! u0nin = u0n / mubar;
%! uminnin = uminn / mubar;
%! vsatnin = vsatn / (mubar * Ebar);
%! 
%! u0pin = u0p / mubar;
%! uminpin = uminp / mubar;
%! vsatpin = vsatp / (mubar * Ebar);
%! 
%! Nrefnin = Nrefn / nbar;
%! Nrefpin = Nrefp / nbar;
%! 
%! Cnin     = Cn / CAubar;
%! Cpin     = Cp / CAubar;
%! 
%! anin     = an / abar;
%! apin     = ap / abar;
%! Ecritnin = Ecritn / Ebar;
%! Ecritpin = Ecritp / Ebar;
%! 
%! % tolerances for convergence checks
%! toll  = 1e-3;
%! maxit = 1000;
%! ptoll = 1e-12;
%! pmaxit = 1000;
%! 
%! % solve the problem using the full DD model
%! [nout, pout, Vout, Fnout, Fpout, Jnout, Jpout, it, res] = ...
%!       secs1d_dd_gummel_map (xin, Din, Nain, Ndin, pin, nin, Vin, Fnin, Fpin, ...
%!                             l2, er, u0nin, uminnin, vsatnin, betan, Nrefnin, ...
%! 	                       u0pin, uminpin, vsatpin, betap, Nrefpin, theta, ...
%! 		               tnin, tpin, Cnin, Cpin, anin, apin, ...
%! 		               Ecritnin, Ecritpin, toll, maxit, ptoll, pmaxit); 
%! 
%! % Descaling procedure
%! n    = nout*nbar;
%! p    = pout*nbar;
%! V    = Vout*Vbar;
%! Fn   = V - Vth*log(n/ni);
%! Fp   = V + Vth*log(p/ni);
%! dV   = diff(V);
%! dx   = diff(x);
%! E    = -dV./dx;
%! 
%! % band structure
%! Efn  = -Fn;
%! Efp  = -Fp;
%! Ec   = Vth*log(Nc./n)+Efn;
%! Ev   = -Vth*log(Nv./p)+Efp;
%! 
%! plot (x, Efn, x, Efp, x, Ec, x, Ev)
%! legend ('Efn', 'Efp', 'Ec', 'Ev')
%! axis tight
