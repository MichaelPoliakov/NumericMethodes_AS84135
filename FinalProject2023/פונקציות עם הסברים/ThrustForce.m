function [Tx,Tz] = ThrustForce(t)
% calculate the thrust force in the X and Z axis
% INPUT: t- time
% OUTPUT: Tx- thrust force in the X axis
%         Tz- thrust force in the X axis

load basic_const.mat tb dm0 ue Ae Pa
load section_e_well_defined_functions.mat Pe

Pe = Pe(t);

T = (t<=tb).*(dm0*ue+Ae.*(Pe-Pa)) ...     % Thrust force.         
  + (t> tb).*0;  

gamma = GammaAngleWithXAxis(t); % [deg]

Tx = T .* cosd(gamma);
Tz = T .* sind(gamma);

end