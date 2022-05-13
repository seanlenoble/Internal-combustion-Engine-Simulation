%% gives diffusuion of air at given temperature (kelvin)
%% for diffustion at 300k
  %% k_diffusion(300)

function out = k_diffusion(T)
  a = 1.52e-4;
  b = 4.42e-5;
  c = 8e-9;
  out = a + b*T + c*T^2;  %[W/m.K]