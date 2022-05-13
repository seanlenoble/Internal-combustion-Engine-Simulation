%% Dont touch this file
%% You input the degree that you want the outputs at
   %% If you give 0 degrees, expect dV, dA, and Up to be 0

function out = geometry_theta(theta, B, L, r, R, Ach, Ap, N, degrees )
  
  %%%Note: Length units only need to be consistent. If you specify length
  %%%      units of in, then velocity is in in/s, area in in^2, volume in in^3.
  %%%      If you specify meters then velocity is in m/s, area in m^2, volume in m^3.
  
  %Inputs:
  %theta: crank angle [radians];
  %B:     Bore [length];
  %L:     Stroke [length];
  %r:     Compression ratio [unitless];
  %R:     Geometry ratio = l/a [unitless];
  %Ach:   Area of the cylinder head [length^2];
    %% Ach > PI * B^2 /4
    %% can have larger area if surface isnt flat
  %Ap:    Area of the piston crown [length^2];
    %% Ap > PI * B^2 /4
    %% can have larger area if surface isnt flat
  %N:     Rotational crank speed [revs per min (RPM)];
    %% used to find instantaneous piston velocity (Up)
  %degrees: Optional input, default value = 0 (theta given in radians).
  %         Set degrees=1 if theta is specified in degrees.
  
  %Ouputs:
  %V:     instantaneous cylinder volume [length^3];
  %dV:    instantaneous rate of change of cylinder volume, dV/d(theta) [length^3/radian];
  %A:     instantaneous area [length^2];
  %dA:    instantaneous rate of change of cylinder area, dA/d(theta) [length^3/radian];
  %Up:    instantaneous piston velocity [length/second];
  
  T_rad = theta;
  %if needed, convert crank angle to radians
  if (degrees)
    T_rad = theta*pi/180;
  end
  
  %calculate displacement volume
  Vd = pi*B^2*L/4;

  phi = sqrt(R^2-(sin(T_rad)).^2);
  V   = ( Vd / (r-1) ) * ( 1 + 0.5*(r-1) * (R + 1 - cos(T_rad) - sqrt(R^2-(sin(T_rad)).^2) ) );
  dV  = 0.5 * Vd * sin(T_rad) .* (1 + cos(T_rad)./sqrt(R^2-(sin(T_rad)).^2 ) );
  A   = Ach + Ap + (pi*B*L/2)  * (R + 1 - cos(T_rad) - sqrt(R^2-(sin(T_rad)).^2) );
  dA  = (pi*B*L/2) * sin(T_rad) .* (1 + cos(T_rad) ./ sqrt(R^2-(sin(T_rad)).^2 ) );
  Up  = pi * L * (N/60) * sin(T_rad) .* (1 + cos(T_rad) ./ sqrt(R^2 - (sin(T_rad)).^2) );
  
  out = [V, dV, A, dA, Up];
end
  
