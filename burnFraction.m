%% Burn fraction page 414
%  weibe function: 1 - exp[-a *((theta-thetas)/theta_d)^n]
%                = 1 - exp[-a * dum^n]
%                = 1 - exp(temp)
function [result] = burnFraction()
  % computes and plots the cumulative burn fraction and instantenous burnrate
  
  a = 5 ;          % weibe efficiency factor
  n = 3  ;         % weibe form factor
  theta_s = -20;  % start of combustion 
  theta_d = 40;   % duration of combustion
  
  theta = linspace(theta_s,theta_s+theta_d,100);     % crankangle theta vector
  dum   = (theta-theta_s)/theta_d;                   % theta difference vector
  temp  = -a*dum.^n;
  
  xb  = 1.-exp(temp);             %burn fraction
  dxb = n*a*(1-xb).*dum.^(n-1);    %element by element vector multiplication
  
  result(:,1)=theta;
  result(:,2)=xb;
  % plot results
  %plot(theta,xb);
  %set(gca, ’fontsize’, 18,’linewidth’,2);
  %xlabel(’Crank Angle (deg)’,’fontsize’, 18);
  %ylabel(’Cumulative Burn Fraction’,’fontsize’, 18);
  %figure();
  
  %plot(theta,dxb);
  %set(gca, ’fontsize’, 18,’linewidth’,2);
  %xlabel(’Crank Angle (deg)’,’fontsize’, 18);
  %ylabel(’Burn Rate (1/deg)’,’fontsize’, 18);
  
end
  