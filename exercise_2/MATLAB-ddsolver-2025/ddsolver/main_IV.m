%%
% SIMULATIONS: V = 0 V und V = 0.7 V
V1 = 0.0;
V2 = 0.7;
V3 = -0.2;

[~, profile_0V, ~, ~] = current4voltage(V1, T0, device, itercontrol);
[~, profile_07V, ~, ~] = current4voltage(V2, T0, device, itercontrol);
[~, profile_02V, ~, ~] = current4voltage(V3, T0, device, itercontrol);


x_um = device.mesh.x * 1e6; % x-Achse in µm

% PLOTTING: n(x)
figure(6);
semilogy(x_um, profile_0V.n, 'b-', 'LineWidth', 2); hold on;
semilogy(x_um, profile_07V.n, 'b--', 'LineWidth', 2);
semilogy(x_um, profile_02V.n, 'b-.', 'LineWidth', 2);
xlabel('Position / {\mu}m');
ylabel('Electron Density n / m^{-3}');
legend('n @ 0 V', 'n @ 0.7 V', 'n @ -0.2 V');
title('Electron Density Profile Comparison');
grid on;

% PLOTTING: p(x)
figure(7);
semilogy(x_um, profile_0V.p, 'r-', 'LineWidth', 2); hold on;
semilogy(x_um, profile_07V.p, 'r--', 'LineWidth', 2);
semilogy(x_um, profile_02V.p, 'r-.', 'LineWidth', 2);

xlabel('Position / {\mu}m');
ylabel('Hole Density p / m^{-3}');
legend('p @ 0 V', 'p @ 0.7 V','p @ -0.2 V');
title('Hole Density Profile Comparison');
grid on;

%%
clear all;

% SETUP
T0 = 300;
device.material = silicon_material_properties(T0);
secs1d_physical_constants;

device.doping.NA = 5E23;
device.doping.ND = 1E23;

device.geometry.length = 50e-6;
device.geometry.p_layer_length = 25e-6;

device.mesh.Nelements = 1000;
device.mesh.x = linspace(0, device.geometry.length, device.mesh.Nelements+1)';
device.mesh.sinodes = 1:length(device.mesh.x);
device.mesh.xm = device.geometry.p_layer_length;

itercontrol.tol = 1e-7;
itercontrol.maxit = 5000;
itercontrol.ptol = 1e-15;
itercontrol.pmaxit = 1000;
itercontrol.prev_guess = false;

% SIMULATIONS for V = -0.2, 0.0, 0.7 V
V_all = [-12.2, 0.0, 1.5];
profiles = cell(1, length(V_all));
x_um = device.mesh.x * 1e6;

for k = 1:length(V_all)
    V = V_all(k);
    [~, profiles{k}, ~, ~] = current4voltage(V, T0, device, itercontrol);

    % Extract current profile
    profile = profiles{k};

    figure(k);
    scale = 1;
    grid on
    plot2micron = 1E6; % scale from meter to micrometer
    set(1,'Position', [13 500 435 320]);
    title({'Energy and Potential Profiles @ T= ',num2str(T0), ' V=', num2str(V),' V'}); 
    hold on;
    plot(device.mesh.x*plot2micron, profile.psi, 'LineWidth',2,...
         'Color', [0.5*scale 0.2*scale 0]); 
    plot(device.mesh.x*plot2micron, profile.EFn,  'LineWidth',2,'Color',...
          [0 0 scale]);
    plot(device.mesh.x*plot2micron, profile.EFp, 'LineWidth',2,'Color',...
         [scale 0 0]);
    plot(device.mesh.x*plot2micron, profile.Ec, 'LineWidth',1,'Color',...
          [0 0 scale]);
    plot(device.mesh.x*plot2micron, profile.Ev, 'LineWidth',1,'Color',...
          [scale 0 0]);
    xlabel('Position / {\mu}m');
    ylabel('Potential or Energy'); 
    legend(['{\Psi} /V'],['{E_{F,n}} /eV'],['{E_{F,p}} /eV'],...
           ['{E_{C}} /eV'],['{E_{V}} /eV'])
    my_legend = legend;
    my_legend.Location = 'northwest';    
    axis tight;
    hold off;

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving 1D Poisson + Drift Diffusion semiconductor eqns using
%                    Gummel algorithm
%
%                 
%               
%
%     The code as is will calculate and plot a JV curve
%     as well as carrier densities, current densities, and electric field
%     distributions of a pn junction.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%----- SETUP SIMULATION -----------------------------------------------------

% set temperature
T0       = 300;     % [K]

%secs1d_silicon_material_properties;
device.material = silicon_material_properties(T0);

% physical constants and parameters
secs1d_physical_constants;

device.doping.NA = 5E23; % [m^3]
device.doping.ND = 1E23; % [m^3]

% set device geometry
device.geometry.length         = 50e-6; % [m]
device.geometry.p_layer_length = 25e-6;  % [m]

% set device mesh

% uniform mesh
device.mesh.Nelements = 1000;
device.mesh.x = linspace (0, device.geometry.length, device.mesh.Nelements+1)';
% number of mesh cells
device.mesh.sinodes= [1:length(device.mesh.x)];
% position of junction
device.mesh.xm = device.geometry.p_layer_length; % device.geometry.length /2;

% set control parameters for simulation flow
% tolerances for convergence checks
itercontrol.tol = 1e-4;
itercontrol.maxit =  5000;
itercontrol.ptol = 1e-15;
itercontrol.pmaxit = 1000;

% if TRUE this parameter enable the use of previously calculate quantities as 
% initial guess for self-consistent solution of p,n, and psi
% if FALSE the thermal equilibrium is chosen as initial guess
% we will stay, for the time being with the thermal equilibrium as initial guess
itercontrol.prev_guess = false;
% set external voltage here!
voltage_step = 0.01;      %[V] 
voltage_start = -0.2;     %[V]
voltage_end = 1.2; %[V]

number_voltages = floor((voltage_end-voltage_start)/voltage_step)+1;


% this array stores all voltages for which the currents will be calculated
voltage_ramp = linspace(voltage_start,voltage_end, number_voltages);
find_zero_voltage = find(voltage_ramp == 0);
if (find_zero_voltage > 0)
    voltage_ramp(find(voltage_ramp == 0)) = voltage_step/4;
end;

current_ramp = zeros(1,length(voltage_ramp));

% as the current might be very small for small voltages, it is useful to calculate
% the current for large voltage first and proceed with decreasing voltages
%    you are welcome to test, whether the simulated current changes if one loops 
%    through voltages in increasing order

for bias_voltage_num = length(voltage_ramp):-1:1

    % ask for tighter convergence for small voltages
    if (voltage_ramp(bias_voltage_num) < 0.2)
        itercontrol.tol = 1e-7;
    else
        itercontrol.tol = 1e-4;
    end;

    itercontrol.prev_guess = false;

    disp(['Calculate j at voltage V=',num2str(voltage_ramp(bias_voltage_num)),' V']);

    [current_ramp(bias_voltage_num), profile, it, res] = current4voltage(voltage_ramp(bias_voltage_num),T0,device,itercontrol);
end
% 
% figure(1)
% scale = 1;
% grid on
% plot2micron = 1E6; % scale from meter to micrometer
% set(1,'Position', [13 500 435 320]);
% title({'Energy and Potential Profiles @ T= ',num2str(T0), ' V=', num2str(V),' V'}); 
% hold on;
% plot(device.mesh.x*plot2micron, profile.psi, 'LineWidth',2,...
%      'Color', [0.5*scale 0.2*scale 0]); 
% plot(device.mesh.x*plot2micron, profile.EFn,  'LineWidth',2,'Color',...
%       [0 0 scale]);
% plot(device.mesh.x*plot2micron, profile.EFp, 'LineWidth',2,'Color',...
%      [scale 0 0]);
% plot(device.mesh.x*plot2micron, profile.Ec, 'LineWidth',1,'Color',...
%       [0 0 scale]);
% plot(device.mesh.x*plot2micron, profile.Ev, 'LineWidth',1,'Color',...
%       [scale 0 0]);
% xlabel('Position / {\mu}m');
% ylabel('Potential or Energy'); 
% legend(['{\Psi} /V'],['{E_{F,n}} /eV'],['{E_{F,p}} /eV'],...
%        ['{E_{C}} /eV'],['{E_{V}} /eV'])
% my_legend = legend;
% my_legend.Location = 'northwest';    
% axis tight;
% hold off;
% 
% figure(2)
% scale = 1;
% grid on
% plot2micron = 1E6; % scale from meter to micrometer
% set(2,'Position', [13 100 435 320]);
% title({'Charge Density Profiles @ T= ',num2str(T0),' V=', num2str(V),' V'}); 
% hold on;
% semilogy(device.mesh.x*plot2micron, profile.n, 'LineWidth',2,...
%      'Color', [0 0 1]); 
% xlabel('Position / {\mu}m');
% ylabel('Density / {m^{-3}} '); 
% legend('n - Electron Density');
% my_legend = legend;
% my_legend.Location = 'northwest';    
% axis tight;
% hold off;

% figure(3)
% scale = 1;
% grid on
% hold on
% plot2micron = 1E6; % scale from meter to micrometer
% set(3,'Position', [500 100 435 320]);
% title({'Charge Density Profiles @ T= ',num2str(T0), ' V=', num2str(V),' V'}); 
% semilogy(device.mesh.x*plot2micron, profile.p,  'LineWidth',2,'Color',...
%       [1 0 0]);
% xlabel('Position / {\mu}m');
% ylabel('Density / {m^{-3}} '); 
% legend('p - Hole Density');
% my_legend = legend;
% my_legend.Location = 'northwest';    
% axis tight;
% hold off


figure(4)
set(4,'Position', [1000 500 435 320]);
title({'Current Density vs Voltage @ T= ',num2str(T0),' K'}); 
hold on;
grid on
plot(voltage_ramp, abs(current_ramp), 'LineWidth', 2)
legend ('J_{tot}' );
xlabel('Voltage  / V');
ylabel('Current Density  / A{m^{-2}} '); 
axis tight;
hold off;

figure(5)
set(5,'Position', [500 500 435 320]);
title({'Current Density vs Voltage @ T= ',num2str(T0),' K'}); 
hold on;
grid on
plot(voltage_ramp, abs(current_ramp),'LineWidth',3);
scatter(voltage_ramp, abs(current_ramp));
legend ('J_{tot}' );
xlabel('Voltage  / V');
ylabel('Current Density  / A{m^{-2}} '); 
%this plots the current on a logarithmic scale
set(gca,'yscale','log');
hold off;

