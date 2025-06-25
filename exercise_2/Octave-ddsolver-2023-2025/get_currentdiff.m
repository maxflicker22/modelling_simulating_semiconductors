function [dcurrent_den] =  get_currentdiff(current_ref,voltage,temperature, device,itercontrol)

% function get_currentden
%
% obtain current density of a device for a given applied voltage and computes 
% difference to a reference currentby
% returns a single value and suppresses any other output
%
% INPUT: 
%          current_ref .. reference current
%          voltage     .. applied external voltage
%          device      .. data record providing all information on mesh, doping,
%                         and material
%          itercontrol .. data record providing parameters to control
%                         iterations
%
% OUTPUT:  current_den .. current density
%

 [current_den, temp0,temp1,temp2] = current4voltage(voltage,temperature,device,itercontrol);                                                                  
 dcurrent_den = current_den - current_ref;
 end
