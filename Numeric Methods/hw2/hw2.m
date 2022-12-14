%%
% 
% % 
% $$\displaystyle \int _{a}^{b}f(x)\,dx\approx {\frac {b-a}{n}}\left({f(a) \over 2}+\sum _{k=1}^{n-1}\left(f\left(a+k{\frac {b-a}{n}}\right)\right)+{f(b) \over 2}\right)$$
% 
% Q1.1
clearvars; clc;

E_o     =       2500;       % [V/m]    elctrical field velocity at the regon
L       =       0.15;       % [m]     leangth.
a       =       sqrt(pi/5); % []       const.
b       =       0.1;        % []       const.

x1 = L/3;
x2 = 2*L/3;
E = @(x) Elect_func_of_x(E_o,a,b,L, x);
fnc = @(x) -1 * E(x);
format shortG
I = my_integral(fnc, x1, x2, 1e6)


%%
% Q1.2

h = [2.5e-3 1e-3 5e-4 2.5e-4 1e-4 5e-5 1e-5];
I_hh = [];
for t=h
    I_hh = [I_hh, my_integral2(fnc,x1,x2, t)];
end
semilogx(h, I_hh, 'kO-',LineWidth=1.2 ,MarkerSize=10)
xlabel("h [m]"); ylabel("\Delta\phi [v]");
yline(-1.188015857351312e+02, 'r-.') % evaluated I with N=10^10 as refrence.
box off
title("Convergence of the Numerical Integral as \DeltaX ---> 0")

%%
% Q1.3

phi_o  =        400;        % [V]      elc' potention on the orgen.
h      =        1e-4;

figure(2)
X = linspace(0, L, 200);
plot(X, E(X),'b')
ylabel('E_x [kV/m]','Color','b')
yticklabels(1.8:0.1:2.5)
hold on
yyaxis right
Phi = zeros(1,length(X));
for p = 1:length(X)
    Phi(p) = phi_o + my_integral2(fnc,0,X(p),h);
end
plot(X, Phi);
ylabel '\phi_x [V]'
xlim([0, L])
xlabel 'x [cm]'
xticks(0:0.02:L)
xticklabels(0:2:15)
legend('E(x) [kV/m]','φ(x) [V]')
title('Electric potential & field')
hold off

Ih =@(p) my_integral2(fnc, 0, L, p*1e-6);
phi_f = phi_o + (2/3) * Ih(1) + (1/3) * Ih(2) % fine value for phi calculatet throw integral.

