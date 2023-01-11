% CSTR model taken from
%
% Klatt, K.U. and S. Engell, "Ruehrkesselreaktor mit Parallel- und Folgereaktion"  In S. Engell, editor,
%    Nichtlineare Regelung - Methoden, Werkzeuge, Anwendungen.  VDI-Berichte Nr. 1026, pages 101-108,
%    VDI-Verlag, Duesseldorf, 1993.
%
% Description
% This model describes the dynamics of a single CSTR with van der Vusse reaction:
%
% A-->B-->C
% 2A-->D
%
% The inlet is pure cyclopentadiene (substance A) that reacts to form cyclopentenol (B) and two
%    unwanted by-products, cyclopentanediol (C) and dicyclopentadiene (D).  Since C and D are
%    unwanted and do not react further, their concentrations are not computed.

function xdot = cstr(t,x)

global u

% Inputs (2)
% Vdot/V (hr^-1)
FF = u(1,1);
% Cooling Heat Transfer from Reactor (kJ/hr)
QdotK = u(2,1);

% States (4)
% Concentration of A in the reactor (mol/L)
cA = x(1,1);
% Concentration of B in the reactor (mol/L)
cB = x(2,1);
% Temperature of reactor fluid (deg C)
theta = x(3,1);
% Temperature of cooling fluid (deg C)
thetaK = x(4,1);

% Arrehnius law parameters
% Reaction A->B
k10  =  1.287e12; % hr^-1
E1   =  -9758.3; % K
H1   =  4.2; % kJ/mol
% Reaction B->C
k20  =  1.287e12; % hr^-1
E2   =  -9758.3; % K
H2   =  -11.0; % kJ/mol
% Reaction 2A->D
k30  =  9.043e09; % hr^-1
E3   =  -8560; % K
H3   =  -41.85; % kJ/mol

% Density of reactor fluid
rho  =  0.9342; % kg/L
% Heat capacity of reactor fluid
Cp   =  3.01; % kJ/kg-K
% Heat transfer coefficient between jacket and reactor
kw   =  4032; % kJ/hr-m^2-K
% Surface area for cooling
AR   =  0.215; % m^2
% Volume of the CSTR
VR   =  10; % L
% Coolant Mass
mK   =  5.0; % kg
% Coolant Heat Capacity
CPK  =  2.0; % kJ/mol-K
% Feed concentration
cA0  =    5.1; % mol/L
% Feed temperature
theta0 =  104.9 + 273.15; % K

TIMEUNITS_PER_HOUR = 3600.0;

k1=k10*exp(E1/(theta));
k2=k20*exp(E2/(theta));
k3=k30*exp(E3/(theta));

xdot(1,1) = (1/TIMEUNITS_PER_HOUR)* (FF*(cA0-cA) - k1*cA - k3*cA*cA); 

xdot(2,1) = (1/TIMEUNITS_PER_HOUR)* (- FF*cB + k1*cA - k2*cB); 

xdot(3,1) =(1/TIMEUNITS_PER_HOUR) * (FF*(theta0-theta) - ...
   (1/(rho*Cp)) *(k1*cA*H1 + k2*cB*H2 + k3*cA*cA*H3) + ...
   (kw*AR/(rho*Cp*VR))*(thetaK -theta)); 

xdot(4,1) = (1/TIMEUNITS_PER_HOUR) * ((1/(mK*CPK))*(QdotK + kw*AR*(theta-thetaK)));
