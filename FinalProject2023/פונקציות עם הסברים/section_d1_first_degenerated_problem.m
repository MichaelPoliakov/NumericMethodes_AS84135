function section_d1_first_degenerated_problem()
% SECTION D.1: take a reduced problem which is easy to solve analyticly and
% compare computer's solver to anlytical.
% generates fig 2
load basic_const.mat g

gamma0 = 45; %[deg]
v0 = 50;     %[m/s]

% initial condition set up
x0 =    0               ;
vx0 =   v0*cosd(gamma0) ;
z0 =    0               ;
vz0 =   v0*sind(gamma0) ;

IC = [x0 vx0 z0 vz0];

% time sampele set up
tf = 100*sind(45)/g;   % known from analytical derivation
step = 1e-4;
t = 0:step:tf;

% ODE's setting:
dydt = @(t,y) [y(2) 0 y(4) -g]';

% numeric
[numeric_y,~] = RK5solver(dydt,t,IC);

% analytic
X = 25*sqrt(2).*t; Z = z0 + 25*sqrt(2).*t-1/2*g.*t.^2;

% ploting results
section_d1_plotting(numeric_y,X,Z)
end