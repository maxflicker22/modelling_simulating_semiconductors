%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Solving 1D Poisson + Drift Diffusion semiconductor eqns using
%                    Gummel algorithm
%     +-----------------------------------------------------------------      
%     The code as is will calculate and plot a V-T curve for a given
%     reference current
%     It mimicks the script 
%     The code relies a nested interval method to retrieve the applied
%     voltage whose associated current density matches the reference 
%     current density.
%     The determination of the current density as a function of the applied
%     voltage is performed with the Gummel algorithm.
%     Each current determination is accompanied by the calculation of carrier densities, current densities, and electric field
%     distributions of a generic pn junction. 
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%----- SETUP SIMULATION -----------------------------------------------------

% set temperature
T0       = 300;     % [K]

% provide a reference current






current_ref_value = 2.51E-2; % [A/m2] here 100mA/m^2
current_ref_value = 3.83; % [A/m2] here 100mA/m^
current_ref_value = 38.97; % [A/m2] here 100mA/m^2





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
device.mesh.Nelements = 100; %1000
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

accuracy_nested_interval = 1d-5;

% set external voltage here!
V_applied = 0.1; %[V]

% setup considered voltage range, will serve as search interval 
voltage_step = 0.0001;      %[V] 
voltage_start = 0;     %[V]
voltage_end = 0.6; %[V] 0.4
number_voltages = floor((voltage_end-voltage_start)/voltage_step)+1;

voltage_int = linspace(voltage_start,voltage_end, number_voltages);
% if a voltage entry == 0: we replace 0 by a close-by, small, yet non-zero value
find_zero_voltage = find(voltage_int == 0);
if (find_zero_voltage > 0)
    voltage_int(find(voltage_int == 0)) = voltage_step/4;
end;



% compute the current density for the voltage applied and provide
% the difference difference with respect to the reference current density

current_diff = get_currentdiff(current_ref_value,V_applied,T0,device,itercontrol);

temperatures = linspace(200,400,100);%linspace(250,400,50);
target_voltages = zeros(length(temperatures),1);

for T_count=1:length(temperatures)

    T0 = temperatures(T_count);
    
    % reset intrinsic carrier density by reloading temperature dependent
    % device parameters
    device.material = silicon_material_properties(T0);
    
    % determine voltage at which the reference current is obtained
    voltage = FindRootNestedIntervals(@(V) get_currentdiff(current_ref_value,V,...
                                      T0,device,itercontrol),... 
                                      voltage_int, mean(voltage_int),...
                                      accuracy_nested_interval*current_ref_value, 40);


    % obtain profiles of all quantities associated to this point of operation: 
    % (voltage, current_ref_value)
    [current,profile, it, res] = current4voltage(voltage,T0,device,itercontrol);
    
    target_voltages(T_count) = voltage;
    % plot a few profiles
    plot2micron = 1e6;
    scale  = 1 - T_count / length(temperatures);
    
    if (or(T_count == 1, T_count == length(temperatures)))       
      
    figure(1)
        set(1,'Position', [13 700 435 320]);
        title({'Potential and energy profiles' }); 
        hold on;
        plot(device.mesh.x*plot2micron, profile.psi, 'LineWidth',2,...
             'Color', [0.5*scale 0.2*scale 0],'DisplayName',['{\Psi} /V']); 
        plot(device.mesh.x*plot2micron, profile.EFn,  'LineWidth',2,'Color',...
              [0 0 scale]);% ,'DisplayName',['{E_{F,n}} /eV']); 
        plot(device.mesh.x*plot2micron, profile.EFp, 'LineWidth',2,'Color',...
              [scale 0 0]);% ,'DisplayName',['{E_{F,p}} /eV']); 
        plot(device.mesh.x*plot2micron, profile.Ec, 'LineWidth',1,'Color',...
              [0 0 scale]);% ,'DisplayName',['{E_{C}} /eV']); 
        plot(device.mesh.x*plot2micron, profile.Ev, 'LineWidth',1,'Color',...
              [scale 0 0]);% ,'DisplayName',['{E_{V}} /eV']);
        xlabel('position / {\mu m}');
        ylabel('potential or energy'); 
        legend(['{\Psi} /V'],['{E_{F,n}} /eV'],['{E_{F,p}} /eV'],...
               ['{E_{C}} /eV'],['{E_{V}} /eV'])
        my_legend = legend;
        %my_legend.Location = 'northeastoutside';    
        axis tight;
        jref_str = strrep(sprintf('%.2e', current_ref_value), '.', 'p');
        % Save figure 1: Potential and energy profiles
        filename1 = sprintf('fig_potential_energy_T%03dK_jref_%s.png', round(T0), jref_str);
        saveas(1, filename1);




        % format current for filename: replace "." with "p" for compatibility

        hold off;
        
    figure(3)
        set(3,'Position', [490 700 435 320]);
        title({'Charge carrier density profiles' }); 
        hold on;
        plot(device.mesh.x*plot2micron, profile.n, 'LineWidth',2,'Color',...
              [0 0 scale]); %,'DisplayName',['n / m{^{-3}}']);
        plot(device.mesh.x*plot2micron, profile.p,'LineWidth',2,'Color',...
              [scale 0 0]); %,'DisplayName',['p / m{^{-3}}']);
        set(gca,'yscale','log');
        xlabel('position / {\mu m}');
        ylabel('potential or energy'); 
        legend(['n / m{^{-3}}'],['p / m{^{-3}}']) 
        my_legend = legend;
        % Save figure 3: Charge carrier densities
        filename3 = sprintf('fig_charge_densities_T%03dK_jref_%s.png', round(T0), jref_str);
        saveas(3, filename3);


    
        %my_legend.Location = 'northeastoutside';    
        hold off;
    
    end; % if 
end % for temperatures T

% analyze Voltage temperature behavior
% get slope dV/dT
dV = diff(target_voltages);
dT = diff(temperatures);
slope = dV ./ dT;

% Define fit range
first_index = min(find(temperatures > 285));
index       = min(find(temperatures > 370));

% Extract data
Vdata = target_voltages(first_index:index);
tdata = temperatures(first_index:index);

% Build design matrix
Tdata = [ones(length(tdata), 1), tdata(:)];

% Perform linear fit
theta = (Tdata' * Tdata) \ (Tdata' * Vdata);

% === Figure 5: Voltage vs Temperature with Linear Fit ===
figure(5)
set(5, 'Position', [100 300 800 600]);  % Larger figure size
clf
title({'Voltage vs Temperature in Si at given j_{ref}'}, 'FontSize', 14); 
hold on;

% Simulated V(T) data
plot(temperatures, target_voltages, ...
    'linestyle', 'none', ...
    'marker', 'o', ...
    'markersize', 8, ...
    'Color', [1 0 0], ...
    'MarkerFaceColor', [1 0 0], ...
    'DisplayName', ['Simulated V(T) @ ', num2str(current_ref_value), ' A/m^{2}']);

% Linear Fit with slope and intercept in legend
slope_fit = theta(2);
bias_fit = theta(1);
fit_label = sprintf('Linear Fit (slope = %.2e V/K, bias = %.2e V)', slope_fit, bias_fit);

plot(Tdata(:,2), Tdata * theta, ...
    'LineWidth', 2, ...
    'Color', [0.5 0 0], ...
    'DisplayName', fit_label);

xlim([temperatures(1), temperatures(end)]);
ylim([min(target_voltages), max(target_voltages)]);
xlabel('Temperature [K]', 'FontSize', 12);
ylabel('Voltage [V]', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
grid on;
axis tight;

jref_str = strrep(sprintf('%.2e', current_ref_value), '.', 'p');
filename5 = sprintf('fig_voltage_temperature_curve_jref_%s.png', jref_str);
saveas(5, filename5);

% === Figure 6: dV/dT vs Temperature ===
% Prepare data for slope plot
Tmid = 0.5 * (temperatures(1:end-1) + temperatures(2:end)); % midpoint T

figure(6)
clf                         % Clear the figure window before plotting
set(6, 'Position', [950 300 800 600]);  % Bigger window
plot(Tmid, slope, 'LineWidth', 2, ...
    'DisplayName', 'dV/dT (numerical)');

title('Slope dV/dT as Function of Temperature', 'FontSize', 14);
xlabel('Temperature [K]', 'FontSize', 12);
ylabel('dV/dT [V/K]', 'FontSize', 12);
%legend('Location', 'best', 'FontSize', 10);  % Only one entry
grid on;
axis tight;

filename6 = sprintf('fig_slope_dVdT_jref_%s.png', jref_str);
saveas(6, filename6);

