%%
%
% QUESTION E:

% const:
m    = 2;               % [kg]
k    = 50;              % [N/m]
l0   = 2;
rho  = 1.225;
A    = pi*(0.15)^2;     % [m^2]
Cd   = 0.1;
g    = 9.81;            % [m·s^-2]
% some more const.:
az0 = 3;            ax0 = 2;                    % [m·s^-2]
t0  = 10;           t1=7;               t2=11;  % [s]
dzD0 = double(0);   dxD0 = double(0);


h=0.0000001;
N = 100;    % max time [s]
t = 0:h:N;

% drone state set up:
zD = ones(1,length(t)); xD = ones(1,length(t));
zD(1)= 10*c.l0; xD(1)= double(0); % [m]

% consruct D(x(t),y(t)) field
for i=2:(length(t)-1)
    tmp = double(0);
    if t(i) <= t0
        tmp = az0*t(i)^2/2;
    elseif t0<t(i) && t(i)<=2*t0
        tmp = -az0*t0^2 + 2*az0*t0*t(i) - az0*t(i)^2/2;
    else
        tmp = az0*t0^2;
    end
    zD(i) = zD(1) + dzD0*t(i) + tmp;

    tmp = double(0);
    if t0<t(i) && t(i)<=2*t0
        tmp = ax0*(t(i)-t1)^2/2;
    else
        tmp = ax0*(t(i)-t2)^2/2 + ax0*(t2-t1)*(t(i)-t2);
    end
    xD(i) = xD(1) + dxD0*t(i) + tmp;
    clear tmp
end

% package state set up:
x = zeros(1,length(t));                       
z = zeros(1,length(t));
vx= zeros(1, length(t));
vz= zeros(1,length(t));

% initial conditions of package:
x(1) = 0;
z(1) = 9*c.l0 - c.m*c.g/c.k;                       % initial condition
vx(1)= 0;
vz(1)= 0;                                          % initial condition

f = readmatrix("U_wind.txt",'Range','B:B');
w=@(z) u_w(z, f);
[x, z, vx, vz] =  Euler_solver(h,t,x,z,vx,vz,xD,zD,w,m,k,l0,rho,A,Cd,g);

figure()

plot(t,x,t,vx,t,z,t,vz)
legend('X_{pkg} [m]','Vx_{pkg} [m/s]','Z_{pkg} [m]','Vz_{pkg} [m/s]')
xlabel 't [sec]'
ylabel 'fun of t'
title 'fig2.1: State of Package as Function of Time.'




