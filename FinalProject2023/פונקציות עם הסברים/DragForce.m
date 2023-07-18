function [Dx,Dz] = DragForce(t,V_2rd)
% calculate the drage force acting on the rocket in the vertical and
% horizontal axis
% INPUT: t- time
%        x_2rd - velocity^2
% OUTPUT: Dx- drag force in X axis
%         Dz- drag force in Z axis
load section_e_well_defined_functions.mat CD

M = Mach(V_2rd); %O(1)
Drag_coeff = CD(M);

rho =       1.225;      % [kg/m^3]  air mass dencity
d =         0.132;      % [m]       diameter
A =         pi*(d/2)^2; % [m^2]     aerodinamical slice surface area

D = -rho*A/2 * Drag_coeff * V_2rd;

gamma = GammaAngleWithXAxis(t);

Dx = D.*cosd(gamma);
Dz = D.*sind(gamma);
end