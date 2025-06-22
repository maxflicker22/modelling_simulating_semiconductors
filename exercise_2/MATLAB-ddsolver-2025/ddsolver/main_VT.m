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
current_ref_value = 1E-5; % [A/m2] here 100mA/m^2

%secs1d_silicon_material_properties;
device.material = silicon_material_properties(T0);

% physical constants and parameters
secs1d_physical_constants;

device.doping.NA = 1E23; % [m^3]
device.doping.ND = 5E23; % [m^3]

% set device geometry
device.geometry.length         = 12e-6; % [m]
device.geometry.p_layer_length = 6e-6;  % [m]

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
accuracy_nested_interval = 1d-5;

% set external voltage here!
V_applied = 0.01; %[V]

% setup considered voltage range, will serve as search interval 
voltage_step = 0.0001;      %[V] 
voltage_start = 0;     %[V]
voltage_end = 0.3; %[V]
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

temperatures = linspace(300,400,100);
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
    
    if (T_count == 1 | T_count == length(temperatures)) 
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
        my_legend.Location = 'northeastoutside';    
        axis tight;
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
        my_legend.Location = 'northeastoutside';    
        hold off;
    end; % if 
end % for temperatures T

% analyze Voltage temperature behavior
% get slope dV/dT
dV = diff(target_voltages);
dT = diff(temperatures);
slope = dV./dT;

figure(5)
set(5,'Position', [13 200 435 320]);
title({'voltage vs temperature in Si at given j_{ref}' }); 
plot(temperatures, target_voltages,'LineWidth',1,'Color',...
         [1 0 0],'DisplayName',...
         ['Voltage @ ',num2str(current_ref_value),' A/m^{2}']);   
hold on;  
%ylim([min(voltage_int) max(voltage_int)]);
xlabel('temperature / K');
ylabel('voltage / V');
axis tight;

figure(6)
set(6,'Position', [490 200 435 320]);
title({'slope dV/dT vs temperature in Si at given j_{ref}' });
plot(temperatures(2:end), slope,'LineWidth',1,'Color',...
         [0.8 0 0]); 
hold on;
legend(['slope @ ',num2str(current_ref_value),' A/m^{2}']);
xlabel('temperature / K');
ylabel('dV/dT / V/K');
axis tight;

