function [outputArg1] = Mach(V_2rd)
% calculate mach number for a certain velocity
% INPUT: v_2rd- velocity^2
% OUTPUT: mach number
ka =        1.4;        % air heat  capacity ratio
R =         287;        % [J/kg·K]  gas constant for air
Temp =      300;        % [K]       air temprature

outputArg1 = sqrt(V_2rd/(ka*R*Temp));
end