function delta_x=section_e7_rocket_harizontal_Delta_x(y)
% find the horizontal distance the rocket reaches
% INPUT: y- the state vector we derived
% OUTPUT: display Delta x.
%        delta_x: returned.

delta_x = abs(y(end,1) - y(1,1));
txt = ['|X2-X1| = ' sprintf('%g',delta_x*1e-3) '[km]'];
disp(txt);
disp(' ');
end