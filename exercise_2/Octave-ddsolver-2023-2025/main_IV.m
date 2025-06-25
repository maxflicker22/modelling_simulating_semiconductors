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

device.doping.NA =  5E23; % [m^3]  5E23
device.doping.ND = 1E23; % [m^3]

% set device geometry
device.geometry.length         = 50e-6; % [m]
device.geometry.p_layer_length = 25e-6;  % [m]

% set device mesh

% uniform mesh
device.mesh.Nelements = 100;
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
voltage_end =    1.2; %[V]

number_voltages = floor((voltage_end-voltage_start)/voltage_step)+1;


% this array stores all voltages for which the currents will be calculated
voltage_ramp = linspace(voltage_start,voltage_end, number_voltages);


disp(['Number of voltage points computed :',num2str(length(voltage_ramp))]);

find_zero_voltage = find(voltage_ramp == 0);
if (find_zero_voltage > 0)
    voltage_ramp(find(voltage_ramp == 0)) = voltage_step/4;
endif;

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
    endif;
    
    disp(['Calculate j at voltage V=',num2str(voltage_ramp(bias_voltage_num)),' V']);
    
    [current_ramp(bias_voltage_num),qprofile, it, res] = ...
    current4voltage(voltage_ramp(bias_voltage_num),T0,device,itercontrol);
      
end


figure(1)
scale = 1;
plot2micron = 1E6; % scale from meter to micrometer
set(1,'Position', [13 500 435 320]);
title({'energy and potential profiles @ T= ',num2str(T0)}); 
hold on;
plot(device.mesh.x*plot2micron, qprofile.psi, 'LineWidth',2,...
     'Color', [0.5*scale 0.2*scale 0]); 
plot(device.mesh.x*plot2micron, qprofile.EFn,  'LineWidth',2,'Color',...
      [0 0 scale]);
plot(device.mesh.x*plot2micron, qprofile.EFp, 'LineWidth',2,'Color',...
      [scale 0 0]);
plot(device.mesh.x*plot2micron, qprofile.Ec, 'LineWidth',1,'Color',...
      [0 0 scale]);
plot(device.mesh.x*plot2micron, qprofile.Ev, 'LineWidth',1,'Color',...
      [scale 0 0]);
xlabel('position / {\mu m}');
ylabel('potential or energy'); 
legend(['{\Psi} /V'],['{E_{F,n}} /eV'],['{E_{F,p}} /eV'],...
       ['{E_{C}} /eV'],['{E_{V}} /eV'])
my_legend = legend;
%my_legend.Location = 'northeastoutside';    
axis tight;
hold off;

figure(2)
scale = 1;
plot2micron = 1E6; % scale from meter to micrometer
set(2,'Position', [13 100 435 320]);
title({'charge density profiles @ T= ',num2str(T0)}); 
semilogy(device.mesh.x*plot2micron, qprofile.n, 'LineWidth',2,...
     'Color', [0 0 1]); 
hold off;
xlabel('position / {\mu m}');
ylabel('density / {m^{-3}} '); 
legend('n');
my_legend = legend;
%my_legend.Location = 'northeastoutside';    
axis tight;
hold on;

figure(3)
scale = 1;
plot2micron = 1E6; % scale from meter to micrometer
set(3,'Position', [500 100 435 320]);
title({'charge density profiles @ T= ',num2str(T0)}); 
semilogy(device.mesh.x*plot2micron, qprofile.p,  'LineWidth',2,'Color',...
      [1 0 0]);
xlabel('position / {\mu m}');
ylabel('density / {m^{-3}} '); 
legend('p');
my_legend = legend;
%my_legend.Location = 'northeastoutside';    
axis tight;


figure(4)
set(4,'Position', [1000 500 435 320]);
title({'current density vs voltage @ T= ',num2str(T0),' K'}); 
hold on;
plot(voltage_ramp, abs(current_ramp),'LineWidth',3);
legend ('Jtot' );
xlabel('voltage  / V');
ylabel('current density  / A{m^{-2}} '); 
axis tight;
hold off;

figure(5)
set(5,'Position', [1000 500 435 320]);
title({'current density vs voltage @ T= ',num2str(T0),' K'}); 
hold on;
plot(voltage_ramp, abs(current_ramp),'LineWidth',3);
scatter(voltage_ramp, abs(current_ramp));
legend ('Jtot' );
xlabel('voltage  / V');
ylabel('current density  / A{m^{-2}} '); 
%this plots the current on a logarithmic scale
set(gca,'yscale','log');
axis tight;
hold off;

 

