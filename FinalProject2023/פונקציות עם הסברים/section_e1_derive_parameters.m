function [Dx, Dz, Tx, Tz, m] = section_e1_derive_parameters(t,V_2rd)
% generate the required parameters for section E1
% INPUT: t- time
%        V_2rd- velocity^2
% OUTPUT: Dx- drag force in X axis
%         Dz- drag force in Z axis
%         Tx- thrust force in X axis
%         Tz- thrust force in Z axis
%         m- mass of the rocket
m = MassOfRocket(t);    % O(1)
[Tx,Tz] = ThrustForce(t); %O(1)
[Dx,Dz] = DragForce(t,V_2rd);

end