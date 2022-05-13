clear all
clc
close all

%Define gas constants
Rsp = 0.287;     %specific gas constant for air         [kJ/kg.K]
k  = 1.35;       %specific heat ratio                   [unitless]
Cv = Rsp/(k-1);  %specific heat capacity const volume   [kJ/kg.K]
Cp = k*Cv;       %specific heat capacity const pressure [kJ/kg.K]

%Define ambient conditions
Pa = 101.325;         %Ambient pressure          [kPa]
Ta = 333.15;          %Ambient temperature       [K]
rho_a  = Pa/(Rsp*Ta); %Ambient density           [kg/m^3]

%Temperature coming out of the intake pipe
T1 = Ta;       %Initial temperature       [K]

%Engine Geometric Parameters
np = 6;          %Number of cylinders       [unitless]
rc = 10;         %compression ratio         [unitless]

%Specify 2 of the following 3: Vd, B and L
Vd = 0.002;                    %Total displacement volume [m^3]
B  = .5*(4*Vd/(np*pi))^(1/3);     %Bore       [m]
L  =  4*B;%4*Vd/(np*pi*B^2);        %Stroke     [m]
a  = L/2;                      %crankshaft length [m]

%Specify either R (geometric ratio) or l (connecting rod length)
R  = 5.5;                      %Geometric ratio, l/a      [unitless]
l  = R*a;                      %Connecting rod length,    [m]

%Operational parameters
Nstart  = 1000;             %crankshaft rotational speed [rpm]
Nend    = 5000;     
Nnumber = 3;                %don't request fewer than 2 or the code will break
theta_s = -20;              %combustion timing angle [degrees]
theta_s = theta_s * pi/180; %convert angle to radians

%Fuel parameters
Qhv = 43400;                %Fuel heating value  [kJ/kg fuel]
AF  = 22;                   %A/F ratio           [unitless]
Qin = Qhv/(AF+1);           %Heat input,         [kJ/kg mixture]
theta_d = 40;               %combustion duration [degrees]
theta_d = theta_d*pi/180;   %convert angle to radians

%Cooling parameters
Tinf= 85+273;                %Coolant/Wall Temperature [K]
Vbdc = (Vd/np)*(rc/(rc-1));  %BDC volume [m^3]
Vtdc = (Vd/np)*(1/(rc-1));   %TDC volume [m^3]
Ach  = 0.25*pi*B^2;          %Cylinder head surface area [m^2]
Ap   = Ach;                  %Piston crown  surface area [m^2]

%Intake geometry patterns
LD_intake = 1000; %L/D for intake [unitless]
epsD      = 0.01; %epsilon/D (roughness) of intake [unitless]
BD_intake = 1;    %B/D_intake [unitless]

%Valve Parameters
n_iv      = 2;    %Number of intake valves    [unitless]
n_ev      = 2;    %Number of exhaust valves   [unitless]
Te_min = 273+800; %minimum expected exhaust temperature [K]
N_max  = 8000;    %max crankshaft spe2d [RPM]

%TASK 3: CALCULATE THE REQUIRED VALVE DIAMETERS HERE *****LECTURE 10 slide 27
Up_avg_max = 2*L*N_max/60; %max average piston speed based on max crankshaft speed [m/s]
c_in   = sqrt(k*Rsp*1000*Ta);     %speed of sound in the inlet pipe
c_ex   = sqrt(k*Rsp*1000*Te_min); %speed of sound in the exhaust pipe
D_iv   = sqrt(4*1.3*B^2*Up_avg_max/c_in/pi)/n_iv;  %Diameter of intake valves  [m]
D_ev   = sqrt(4*1.3*B^2*Up_avg_max/c_ex/pi)/n_ev;  %Diameter of exhaust valves [m]


%Vectors to store the results as a function of RPM
N_results     = linspace(Nstart,Nend,Nnumber);
imep_results  = 0*N_results;
bmep_results  = 0*N_results;
qmep_results  = 0*N_results;
pmep_results  = 0*N_results;
fmep_results  = 0*N_results;
Tb_results    = 0*N_results;
Pb_results    = 0*N_results;
m_results     = 0*N_results;
etav_results  = 0*N_results;
Up_results    = 0*N_results;
DP_im_results = 0*N_results;
isfc_results  = 0*N_results;
bsfc_results  = 0*N_results;
eta_m_results = 0*N_results;

i = 1;

for N = Nstart: (Nend-Nstart)/(Nnumber-1) : Nend
    disp(N);
    Up_avg = 2*L*N/60;
    
    %TASK 3: CALCULATE INTAKE MANIFOLD PRESSURE DROP FROM TURBULENT WALL FRICTION
    DP_im  = 0.5*rho_a*(Up_avg^2)*f_haaland_turbulent(epsD)*LD_intake/1000; %Pressure drop in intake manifold [kPa] divided by 1000 to go from Pa to kPa
    DP_im_results(i) = DP_im;  %store result to plot later
    P_im   = Pa-DP_im;         % [kPa]    %update pressure: P_im is the pressure just before the intake valves

    
    %TASK 3: CALCULATE INTAKE VALVE PRESURE DROP DP_iv,
    %      EXHAUST VALVE PRESSURE DROP DP_ev,
    %      EXHAUST SYSTEM PRESSURE DROP DP_es
    DP_iv  = (1/1000)*Cv*((P_im/Pa)*(Up_avg/n_iv)*(B^2/D_iv^2))^2;  % [kPa]
    DP_ev  = (1/1000)*Cv*((P_im/Pa)*(Up_avg/n_ev)*(B^2/D_ev^2))^2   % [kPa]
    DP_es  = 0.178*(Up_avg*P_im/Pa)^2;                              % [kPa]
    P1     = P_im - DP_iv;  % [kPa] %update pressure: P1 is the pressure after the intake valve, before compression
    pmep   = DP_im + DP_iv + DP_ev + DP_es;   % [kPa]
    W_p    =(pmep*Vd);
    m      = P1*Vbdc/(Rsp*T1); %  [kg]
    
    
    [imep,fmep,bmep,Tb,Pb,qmep,theta,P,V,isfc,bsfc,eta_m] = integration_driver(Pa, P1, T1, Vd/np, np, B, R, k, rc, N, theta_s, Rsp, Qin, theta_d, L, Vbdc, Ach, Ap, m, Cv, Cp, Tinf, AF, W_p, pmep);
    
    %store the results in the vectors for plotting
    imep_results(i) = imep;
    bmep_results(i) = bmep-pmep;
    Tb_results(i)   = Tb;
    Pb_results(i)   = Pb;
    qmep_results(i) = qmep;
    pmep_results(i) = pmep;
    fmep_results(i) = fmep;
    m_results(i)    = m;
    etav_results(i) = P1/(Rsp*T1*rho_a);
    Up_results(i)   = Up_avg;
    isfc_results(i) = isfc;
    bsfc_results(i) = bsfc;
    eta_m_results(i)= eta_m;
    i=i+1;
    
end    



%plot results after calculations are done for all values of N
figure
plot(N_results, bmep_results, '-b', N_results, imep_results, '-r')
title('MEP')
xlabel('Crankshaft Rotational Speed [RPM]')
ylabel('MEP [kPa]')
legend('BMEP', 'IMEP')

plot(N_results, fmep_results, '-b', N_results, qmep_results, '-r', N_results, pmep_results, '-g')
title('MEP')
xlabel('Crankshaft Rotational Speed [RPM]')
ylabel('MEP [kPa]')
legend('FMEP', 'QMEP', 'PMEP')

figure
plot(N_results, Tb_results, '-')
title('Brake Torque')
xlabel('Crankshaft Rotational Speed [RPM]')
ylabel('Brake Torque [kN.m]')

figure
plot(N_results, Pb_results, '-')
title('Brake Power')
xlabel('Crankshaft Rotational Speed [RPM]')
ylabel('Brake Power [kW]')

figure
plot(N_results, etav_results, '-')
title('Volumetric Efficiency')
xlabel('Crankshaft Rotational Speed [RPM]')
ylabel('Volumetric EFficiency')

figure
plot(N_results, Up_results, '-')
title('Average Piston Velocity')
xlabel('Crankshaft Rotational Speed [RPM]')
ylabel('Velocity')

%%%%%%% Additional Plots %%%%%%%%

figure
plot(N_results, isfc_results, '-b', N_results, bsfc_results, '-r')
title('BSFC and ISFC')
xlabel('Crankshaft Rotational Speed [RPM]')
ylabel(' [g/kW-hr]')
legend('ISFC', 'BSFC')

figure
plot(N_results, eta_m_results, '-b')
title('Mechanical Efficiency')
xlabel('Crankshaft Rotational Speed [RPM]')

