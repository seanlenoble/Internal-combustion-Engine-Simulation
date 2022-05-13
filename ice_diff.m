%% substantial changes to be made here

function out = ice_diff(theta,y,k,B,L,rc,R,N,theta_s, Ach, Ap, m, Rsp, Qin, theta_d, Tinf)
  %Inputs:
  %theta: independent variable, crankshaft angle [radians]
  %y    : vector of dependent variables, y(1) is pressure [kPa], y(2) is
  %work [kJ/kg]
  %B    : Bore   [m]
  %L    : Stroke [m]
  %rc   : compression ratio [unitless]
  %R    : geometric factor  [l/a]
  %N    : RPM
  %theta_s : combustion start angle [radians]
  %Ach  : crown  head area [m^2]
  %Ap   : piston area      [m^2]
  %m    : mass in one cylinder [kg]
  %Rsp  : specific gas constant [kJ/kg.K]
  %Qin  : heat input per kg mix [kJ/kg mixture]
  %theta_d : combustion duration angle [radians]
  %Tinf : cylinder wall temperature/cooling temperature [K]
  
  n = 3;
  a = 5;
  Up_avg = 2*L*N/60;
  P = y(1);   
  W = y(2);   
  out_geometry = geometry_theta(theta, B, L, rc, R, Ach, Ap, N, 0);
  V  = out_geometry(1);      %instantaneous cylinder volume
  dV = out_geometry(2);      %dV/d(theta)
  A  = out_geometry(3);      %instantaneous cylinder area
  dA = out_geometry(4);      %dA/d(theta)
  Up_inst = out_geometry(5); %instantaneous piston velocity
  
  %TASK 2: CALCULATE HEAT LOSS THROUGH MODEL HERE
  Tg   = (P*V)/(m*Rsp);    %gas temperature   %%temp in cylinder
  rho  = P/(Rsp*Tg);       %density           %%density in cylinder
  k_T  = k_diffusion(Tg);  %heat conduction coefficient
  mu_T = mu(Tg);           %viscosity
  Re   = abs((rho*Up_avg*B)/mu_T);           %Reynolds number
                             %(keep the absolute value in case your
                             %velocity is negative)
 
  Tw = 85+273;          %%cylinder wall temp 
  Nu   = 0.49*Re^(0.7); %%0.277    %Nusselt number
  h = k_T/B*Nu;
  
  dQdt = h*A*(Tg-Tinf);          %dQ/dt due to heat transfer [W/m^2.K * m^2*K = W]
  dQ_theta = dQdt*(60/(2*pi*N)); %dQ/dtheta due to heat transfer 
  
  
  %TASK 1: INCLUDE THE FINITE RATE COMBUSTION MODEL HERE
  dQ_comb = 0.0;  %%leave this as 0 for now
  
  if (theta>=theta_s) && (theta<=theta_s+theta_d)

     x_b = 1-exp(-a*((theta-theta_s)/theta_d)^n);
     dQ_comb = n * a *m*Qin*(1-x_b)/theta_d*((theta-theta_s)/theta_d)^(n-1);

  end
  
  
  %TASKS 1 and 2: YOU WILL NEED TO AUGMENT THE MODEL EQUATION FOR PRESSURE 
  %TO INCLUDE THE TERMS FOR HEAT TRANSFER AND COMBUSTION
   
  dP = -k*P*dV/V + ((k-1)/V)*(dQ_comb - (dQ_theta/1000));   
  dW = P*dV;
  out = [dP;dW;dQ_theta/1000];
  
end