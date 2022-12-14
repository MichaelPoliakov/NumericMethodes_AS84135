function [Y, R] = RK5solver(F,t, IC,h)
arguments
    F
    t (1,:) {mustBeNonnegative, mustBeReal}
    IC (:,1)
    h double {mustBeScalarOrEmpty} = t(2) - t(1)
end
N = length(IC);
M = length(t);
Y = zeros(N,M);
Y(:,1) = IC;

R = zeros(N,M);

for i=1:(M-1)                              % calculation loop
    K1 = h.* F( t(i)                ,       Y(:,i)                                                   );
    K2 = h.* F( t(i) + 1/3*h ,       Y(:,i) + 1/3 .* (K1)                           );
    K3 = h.* F( t(i) + 1/3*h ,       Y(:,i) + 1/6 .* (K1+K2)                   );
    K4 = h.* F( t(i) + 1/2*h ,       Y(:,i) + 1/8 .* (K1+3.*K3)              );
    K5 = h.* F( t(i) + h         ,       Y(:,i) + 1/2 .* (K1-3.*K3+4.*K4)  );
   
    Y(:,i+1) = Y(:,i) + 1/6.*(K1+4.*K4+K5);  % main equation
    R(:,i+1) = abs(1/30.*(2.*K1-9.*K3+8.*K4-K5));
end



end
