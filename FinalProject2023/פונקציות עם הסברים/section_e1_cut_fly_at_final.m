function [t_proper,y_proper] = section_e1_cut_fly_at_final(t,y)
% cuts the solution when z reaches -h
% INPUT: t- time
% OUTPUT: y- the state vector we derived

load basic_const.mat h
z = y(:,3);

t_proper = t(z>-h);
y_proper = y(z>-h,1:end);
end