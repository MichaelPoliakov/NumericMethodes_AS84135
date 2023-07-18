function [Y, R] = RK5solver(F,t,IC,h)
% find the function described by the ODE and the error size for a certian
% moment, initial conditions and step size
% INPUT: F- the derivative of function Y
%        t- time
%        IC- initial conditions
%        h- step size
% OUTPUT: Y- the function
%         R- the error
arguments
    F
    t  (1,:) {mustBeNonnegative, mustBeReal}
    IC (1,:)
    h  double {mustBeScalarOrEmpty} = t(2) - t(1)
end

col = numel(IC);
row = numel(t);
Y = zeros(row,col,'double');
Y(1,:) = IC;

R = 0.*Y;

for i=1:(row-1)                              % calculation loop
    tau = t(i);

    K1 = h .* F( tau         ,       Y(i,:)                            )';
    K2 = h .* F( tau + 1/3*h ,       Y(i,:) + 1/3 .* (K1)              )';
    K3 = h .* F( tau + 1/3*h ,       Y(i,:) + 1/6 .* (K1+K2)           )';
    K4 = h .* F( tau + 1/2*h ,       Y(i,:) + 1/8 .* (K1+3.*K3)        )';
    K5 = h .* F( tau + h     ,       Y(i,:) + 1/2 .* (K1-3.*K3+4.*K4)  )';
   
    dY = 1/6 .* (K1+ 4.*K4 +K5);

    Y(i+1,:) = Y(i,:) + dY ; % main equation
    R(i+1,:) = abs(1/30.*(2.*K1 -9.*K3 +8.*K4 -K5));
end

end
