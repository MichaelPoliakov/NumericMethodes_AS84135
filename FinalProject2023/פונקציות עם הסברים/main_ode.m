function dydt = main_ode(t,y)
% calculate the derivative in time of the state equations we derived
% INPUT t- time sampels scalar
%       y- state vector 4x1
% OUT: dydt- d(y)/dt ODE system of 1st order.

load basic_const.mat g

Vsqrd = y(2).^2+y(4).^2;
[Dx, Dz, Tx, Tz, m] = section_e1_derive_parameters(t,Vsqrd);

dydt = zeros(4,1);
    
dydt(1) = y(2)               ;
dydt(2) = Dx ./ m + Tx ./ m    ;
dydt(3) = y(4)               ;
dydt(4) = Dz ./ m + Tz ./ m - g;

end