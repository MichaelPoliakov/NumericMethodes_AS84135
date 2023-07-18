function m = MassOfRocket(t)
% calculate the mass of the rocket in a certain moment
% INPUT: t- time
% OUTPUT: m- mass
load basic_const.mat m0 mf tb dm0

m = (t<=tb) .* (m0-dm0.*t) ...                % rocket mass.                     
   +(t>tb ) .* mf;

end