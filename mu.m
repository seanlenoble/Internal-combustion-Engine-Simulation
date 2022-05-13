%% Viscosity of air as a function of temperature (kelvin)
%% if you want viscosity of air at 300k 
  %% mu(300)
  
function out = mu(T)
  out = 3.3e-7 * T^(0.7); %[kg/m.s]