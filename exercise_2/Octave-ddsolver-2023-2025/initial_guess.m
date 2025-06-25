%
% usage: 
%     ini_profile = initial_guess(voltage, device, temperature, profile)
%     n_ini = ini_profile.n
%     Fp_ini = ini_profile.EFp ...
%
%  take_prev
%
function ini_profile = initial_guess(take_prev, V_p, device, temperature, profile)
  
  % if profile exists  take initial guess from profile
  % if no profile availabe -> generate one and assign initial values based on
  % thermal equilibrium
  
  secs1d_physical_constants;
  Vth 	 = Kb * temperature / q;
  
  mesh_left = device.mesh.x <= device.mesh.xm;
  mesh_right = device.mesh.x > device.mesh.xm;
    
 
  ni = device.material.ni;
    
  % quasi Fermi levels
  Fp = V_p * (mesh_left);
  Fn = Fp;

  %disp(nargin);
  %disp(take_prev);
  
  
  if ((!take_prev) || (nargin ==4))
    profile_elements = device.mesh.Nelements;

    ini_profile.n   = zeros(profile_elements,1);
    ini_profile.p   = zeros(profile_elements,1);
    ini_profile.Jn  = zeros(profile_elements,1);
    ini_profile.Jp  = zeros(profile_elements,1);
    ini_profile.EFn = zeros(profile_elements,1);
    ini_profile.EFp = zeros(profile_elements,1);
    ini_profile.Ec  = zeros(profile_elements,1);
    ini_profile.Ev  = zeros(profile_elements,1);
    ini_profile.psi = zeros(profile_elements,1);
    ini_profile.efield = zeros(profile_elements,1);
    
    % account for doping profile [m^{-3}]
    D = device.doping.ND * (mesh_right)  - device.doping.NA * (mesh_left);
    
    % charge carrier densities
    p = abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni./abs(D)) .^2)) .* (mesh_left) + ...
    ni^2 ./ (abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^2))) .* (mesh_right);

    n = abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^ 2)) .* (mesh_right) + ...
    ni ^ 2 ./ (abs (D) / 2 .* (1 + sqrt (1 + 4 * (ni ./ abs (D)) .^2))) .* (mesh_left);
    
    % electrostatic potential
    V = Fn + Vth * log(n / ni); 
  
    ini_profile.n = n;
    ini_profile.p = p; 
    
    ini_profile.psi =V;
  else

    ini_profile = profile;
    disp(['@V=',V, 'V, took profile']);
    % inconsistent
    % ini_profile.psi =  Fn + Vth * log(ini_profile.n / ni); 
     

  endif
  
  ini_profile.EFn = Fn;
  ini_profile.EFp = Fp; 
    
endfunction
