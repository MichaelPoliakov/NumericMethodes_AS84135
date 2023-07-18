% save constats for entire project

m0 =        42;         % [kg]
dm0 =       1;          % [kg/sec]
tb =        20;         % [sec]        
mf =        22;         % [kg]
% Basic Constants:

Ae =        pi*(0.0200)^2;  % [m^2] exit area surface
At =        pi*(0.0025)^2;  % [m^2] nezzle area surface
Pa =        101325;         % [Pa] suaraunding pressure
kp =        1.28;           % propellant heat capacity ratio
save pressure_const.mat Ae At Pa kp tb

h =         240;        % [m]
g =         9.81;       % [N/kg]

ue =        1800;       % [m/sec]

save basic_const.mat
fprintf("Create 'basic_const.mat' file in %s", pwd);