% Finite Heat release
% gas cycle heat release for 1 engine

function [] = finiteHeatRelease()
  theta_s = -10;     % start of heat release
  theta_d =  40;     % end of heat release
  
  r = 10;      % compression ratio
  k = 1.4;
  q = 34.8;    % dimensionless total heat release Qin/P1V1
  a = 5;       % weibe parameter
  n = 3;       % weibe exponent
  
  step = 1;         % crankangle interval
  NN   = 360/step;  % number of data points
  
  % initialize results
  results.theta    = zeros(NN,1);
  results.vol      = zeros(NN,1);
  results.pressure = zeros(NN,1);
  results.work     = zeros(NN,1);
  
  pinit(1) = 1;     %inital dimensionless pressure P/P1
  theta = -180;           %initial crankangle
  theta_e = theta + step; %final crankangle in step
  fy(1) = pinit(1);       %assign inital pressure to working vector
  fy(2) = 0;              %reset work vector
  
  %pressure and work calculations
  
  for i=1:NN
    [fy,vol] = integrate(theta,theta_e,fy);
    %reset to next interval
    theta = theta_e;
    theta_e = theta+step;
    %copy results to output vector
    results.theta(i) = theta;
    results.vol(i) = vol;
    results.pressure(i) = fy(1);
    results.work(i) = fy(2);
    
  endfor
  
  [pmax1, id_max1] = max(results.pressure(:,1)); %engine max pressure
  thmax1 = results.theta(id_max1);               %engine crank angle
  
  w1 = results.work(NN,1);
  eta1 = w1/q;
  imep1 = eta1*q*(r/(r-1));
  
  % output overall results
##  fprintf(’ Theta_start %5.2f  \n’, thetas(1,1));
##  fprintf(’ Theta_dur %5.2f  \n’, thetad(1,1));
##  fprintf(’ P_max/P_1 %5.2f  \n’, pmax1);
##  fprintf(’ Theta_max %7.1f  \n’,thmax1);
##  fprintf(’ Net Work/P1V1 %7.2f  \n’, w1);
##  fprintf(’ Efficiency %5.3f \n’, eta1);
##  fprintf(’ Imep/P1 %5.2f \n’, imep1);
##  
  
  %plot results
  plot(results.theta,results.pressure(:,1))
  %set(gca, ’fontsize’, 18,’linewidth’,2);
  %legend(’Engine 1’,’Location’,’NorthWest’)
  %xlabel(’Theta (deg)’,’fontsize’, 18)
  %ylabel(’Pressure (bar)’,’fontsize’, 18)
  print -deps2 heatrelpressure
  Computer Programs 417
  figure( );
  plot(results.theta,results.work(:,1))
  %set(gca, ’fontsize’, 18,’linewidth’,2);
  %legend(’Engine 1’,’Location’,’NorthWest’)
  %xlabel(’Theta (deg)’,’fontsize’, 18)
  %ylabel(’Work’,’fontsize’, 18)
    
  
  function[fy,vol] = integrate(theta,theta_e,fy)
  %ode45 integration of the pressure differential equation
  %from theta to thetae with current values of fy as initial conditions
  
  [tt,yy] = ode45(@rates, [theta theta_e], fy);
  fy = yy(length(tt))
  
  %pressure differential equation
  function [yprime] = rates(theta,fy)
    vol=(1.+ (r -1)/2.*(1-cosd(theta)))/r;
    dvol=(r - 1)/2.*sind(theta)/r*pi/180.; %dvol/dthet  
    dx = 0; %set heat release to zero
    
      if (theta > theta_s)
        dum1 = (theta = theta_s)/theta_d;
        x=1.-exp(-(a*dum^n));
        dx = (1-x)*a*n*dum^(n-1).theta_d;  %dx/dtheta
      endif
    term1 = -k*fy(1)*dvol/vol;
    term2 = (k-1)*q*dx/vol;
    yprime(1,1) = term1 + term2;
    yprime(2,1) = fy(1)*dvol;
  endfunction

  end 
  
  
  
 end