function [imep,fmep,bmep,Tb,Pb,qmep,theta,P,V,isfc,bsfc,eta_m] = integration_driver(Pa, P1, T1, Vd, np, B, R, k, rc, N, theta_s, Rsp, Qin, theta_d, L, Vbdc, Ach, Ap, m, Cv, Cp, Tinf, AF, W_p, pmep)
%%integration driver
%inputs:
%Pa : Ambient pressure          [kPa]
%P1 : Initial pressure          [kPa]
%T1 : Initial temperature       [K]
%Vd : Total displacement volume [m^3]
%np : Number of cylinders       [unitless]
%B  : Bore                      [m]
%R  : Geometric ratio, l/a      [unitless]
%k  : specific heat ratio       [unitless]
%rc : compression ratio         [unitless]
%N  : crankshaft rotational speed        [rpm]
%theta_s : combustion start timing angle [radians]
%Rsp : specific gas constant             [kJ/kg.K]
%Qin : Heat input                        [kJ/kg mixture]
%theta_d : combustion duration angle     [radians]
%L   : Stroke               [m]
%Vbdc: BDC volume           [m^3]
%Ach : Cylinder head area   [m^2]
%Ap  : Piston crown area    [m^2]
%m   : mass in the cylinder [kg]
%Cv  : Specific heat at constant volume   [kJ/kg.K]
%Cp  : Specific heat at constant pressure [kJ/kg.K]


%Integrate compression to spark  
ICs   = [ P1, 0, 0];        %P1, W=0, Q=0    %%Q is amount of heat lost so far
tspan = [-pi,theta_s];      %Span of angle, compression
[t,y] = ode45(@(t,y) ice_diff(t,y,k,B,L,rc,R,N,theta_s, Ach, Ap, m, Rsp, Qin, theta_d, Tinf), tspan,ICs);
out_geometry = geometry_theta(t, B, L, rc, R, Ach, Ap, N, 0); %%why call here and inside ice_diff
theta_comp = t;    
P_comp     = y(:,1);
V_comp     = out_geometry(:,1);
V2 = V_comp(end);
P2 = y(end,1);
W2 = y(end,2);
Q2 = y(end,3);
T2 = P2*V2/(m*Rsp); 

%Integrate combustion phase
%Note for task 1: remember to change the limits of integration so that you
%perform the integration in 3 steps:
%(1) -pi to combustion start
%(2) combustion start to combustion end
%(3) combustion end to pi

%% start from P2 not P1
%% Start from W2

%%%%% Everything here must be changed to mirror the integration above
ICs   = [ P2, W2, Q2];  %P1, W=0, Q=0    %%Q is amount of heat lost so far
tspan = [theta_s,theta_s+theta_d];      %Span of angle, compression
[t,y] = ode45(@(t,y) ice_diff(t,y,k,B,L,rc,R,N,theta_s, Ach, Ap, m, Rsp, Qin, theta_d, Tinf), tspan,ICs);
out_geometry = geometry_theta(t, B, L, rc, R, Ach, Ap, N, 0);
theta_comb = t;       %% confusing, where is t declaired? 
P_comb     = y(:,1);
V_comb     = out_geometry(:,1);
V3 = V_comb(end);
P3 = y(end,1);
W3 = y(end,2);    %% wont be true with combustion model
Q3 = y(end,3);    
T3 = P3*V3/(m*Rsp);
%%%%%%


%Integrate expansion phase
ICs   = [ P3, W3, Q3]; %P1, W=W_comp, Q=Q_comp
tspan = [theta_s+theta_d, pi]; %Span of angle, expansion
[t,y] = ode45(@(t,y) ice_diff(t,y,k,B,L,rc,R,N,theta_s, Ach, Ap, m, Rsp, Qin, theta_d, Tinf), tspan,ICs);
out_geometry = geometry_theta(t, B, L, rc, R, Ach, Ap, N, 0);
theta_exp = t;
P_exp     = y(:,1);
V_exp     = out_geometry(:,1);
V4 = V_exp(end);
P4 = y(end,1);
W4 = y(end,2);
Q4 = y(end,3);
T4 = P4*V4/(m*Rsp);

%PRINT THE TEMPERATURE AND PRESSURE AT 4 MAIN STATES
%T1,P1: initial temperature and pressure before compression
%T2,P2: T and P before start of combustion
%T3,P3: T and P after end of combustion
%T4,P4: T and P at the end of the cycle
fprintf('    T1         T2         T3         T4      \n')
fprintf('%10.5f %10.5f %10.5f %10.5f\n',T1,T2,T3,T4)
fprintf('    P1         P2         P3         P4      \n')
fprintf('%10.5f %10.5f %10.5f %10.5f\n',P1,P2,P3,P4)

Wi_total    = np*W4;    %%total indicated work
Qloss_total = np*Q4;
%Compute imep, indicated torque and indicated power
Ti = Wi_total/(4*pi);   %indicated torque [kN.m]
Pi = (N/60)*Wi_total/2; %indicated power  [kW]
imep = Wi_total/Vd;     %imep             [kPa]
qmep = Qloss_total/Vd;  %qmep             [kPa]

%TASK 4: IMPLEMENT FRICTION LOSSES HERE
Up_avg = 2*L*N/60;

%% Models to be found in textbook
%CONSIDER SKIRT, RING AND GAS LOADING

%Skirt Friction 
c_ps = 0.294; %[kPa/s] coefficient suggested by Patton et al. equation (10.27)
fmep_skirt = c_ps*Up_avg/B; % [kPa]

%Ring Friction
c_pr = 0.0406; %[kPa/m2]
fmep_rings = c_pr*(1+1000/N)/B^2; % [kPa] 

%Gas Loading
c_gas = 6.89;   %[unitless] 
K_gas = 0.0238; %[s/m]
fmep_gas = c_gas*((0.088*rc)+0.182*(rc^(1.33-(K_gas*Up_avg))))*P1/Pa; % [kPa]


fmep    = fmep_skirt + fmep_rings + fmep_gas;  % [kPa]
Wdot_f  = fmep*(N/60)*(Vd/2);                  % [kW]
bmep    = imep - pmep - fmep;                  % [kPa]
Tb      = bmep * Vd/(4*pi);                    % [kN.m]
Pb      = (N/60)* bmep * Vd/2;                 % [kW]

%% TASK 5
m_fuel_single = (1000*(m))/(1+AF); %every cycle per cylinder [g]
m_fuel        = m_fuel_single*6;   %every cycle 6 cylinders
mdot_fuel     = m_fuel*((N*60)/2); %total mdot fuel in g/hr

Wdot_i = (imep*(N/60))*(Vd/2); %[kW]
W_in   = Wdot_i*(60/N)*2;      %[kW/cycle]
bsfc   = mdot_fuel/Pb;         %mdot_fuel/Brake Power(Wb)
isfc   = mdot_fuel/Wdot_i;     
Wdot_ig= W_in-W_p;
eta_m=(bmep)/(imep-pmep);

      
%Assemble Output Data
theta = [theta_comp;theta_comb;theta_exp];
P     = [P_comp;P_comb;P_exp];
V     = [V_comp;V_comb;V_exp];


%if you uncomment the following blocks, the function will plot
%the P(theta), V(theta) and P-V plots for each value of N
%You probably want to select a low number of N in driver to not
%be overwhelmed by plots.
%Something like start = 1000, stop = 3000, N = 2


 figure
 plot(theta,P)
 xlabel('Crank Angle [rad]')
 ylabel('Pressure [kPa]')
 
 figure
 plot(theta,V)
 xlabel('Crank Angle [rad]')
 ylabel('Volume [m^3]')

 figure
 plot(V,P)
 title(['Pressure Vs Volume at ', num2str(N), ' rpm'])
 xlabel('Volume [m^3]')
 ylabel('Pressure [kPa]')


