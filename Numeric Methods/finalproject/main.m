
%%
clc; close all; clear all;
set(groot, 'defaultFigureWindowState', 'maximized');
format shortG
tic
%% const:

        m   =  2;                          % [kg]
        k    = 50;                         % [N/m]
        l0   = 2;                           % [m]
        rho  = 1.225;
        R = 0.15;                       % [m]
        A    = pi*R^2;               % [m^2]
        Cd   = 0.1;
        g    = 9.81;                     % [m·s^-2]


t0 = 0;
tf   = 100;
dt  = 5e-4;
tVec = (t0:dt:tf)';



%% Analytical vs Numerical

x0   = 0        ;                                       % initial condition                                          
z0   = 9*l0  ;                                        % initial condition
vx0 = 0       ;                                        % initial condition
vz0 = 0       ;                                        % initial condition
IC = [x0 z0 vx0 vz0]';                        % initial condition vector

% Y = [Y(1) Y(2) Y(3) Y(4)]' = [x z vx vz]'
Ydot = @(t, y) ODEfcn(t, y,@(t)0,@(t)10*l0,@(z)0,m,k,l0,rho,A,Cd,g);
% [t, Y] = ode113(Ydot, tVec, IC);
[Y, ~] = RK5solver(Ydot, tVec,IC);

om2 = k/m;
Am = g/om2;  
B = 9*l0-Am;
analytical_sol = @(t) Am .* cos(sqrt(om2) * t) + B;
% n= 1+10/dt;
z10 = analytical_sol(10);
disp(['Analytical Solution in t=10 sec is:  z(10) = ', num2str(z10), '[m]'])

%% Ploting analytical vs numerical, anlys which h better

z = Y(2,:);
[fig11, fig12] = plot_fig11(tVec,z,analytical_sol(tVec),t0,tf,B-1.05*Am,B+1.05*Am,z10);

[fig131,fig132,fig133] = plot_fig13(analytical_sol,@RK5solver,Ydot, ...
                                                                IC, t0, tf,B-1.05*Am,B+1.05*Am);
% [fig141, fig142, fig143]=plot_fig14(analytical_sol,@RK5solver, Ydot, IC,t0, tf); %very long run

%%
clearvars x0 vx0 z0 z vz0  IC Ydot Y R
clearvars om2 Am B analytical_sol  z

%% Q-5

% initial conditions of package:
x0  = 0;
z0  = 9*l0-m*g/k;                              
vx0 = 0;
vz0 = 0;
IC = [x0 z0 vx0 vz0]';
disp('IC & Time vectors ready') 

% Wind velocity given data & my interpolation fcn apply.
U_wind = readmatrix("U_wind.txt",'Range','B:B');
w=@(z) u_w(z, U_wind); disp('U_w loaded');

C = {'m=', 'k=', 'l0=', 'ρ=', 'A=', 'C_d=', 'g=','t_0=','t_f=','h=';
          m,     k,     l0,       rho,  A,    Cd,        g,    t0,     tf,      dt  };
disp(['Recall: ', sprintf('  %s%g', C{:})]);

%---------------------------------------------- ODE System Explain ---------------------------------------------
%
% Y = [Y(1) Y(2) Y(3) Y(4)]' = [x z vx vz]
% dY/dt = f(t, Y)
%                
%               ⎡  ⚀  vx 
% dY/dt =⎢   ⚁ vz 
%               ⎢  ⚂  k/m*(x-xD)*(l0/sqrt((z-zD)^2+(x-xD)^2)-1)-rho*A*Cd/(2*m)*(vx-w)^2*sign(vx-w)
%               ⎣  ⚃  k/m*(z-zD)*(l0/sqrt((z-zD)^2+(x-xD)^2)-1)-g
%          
% Ydot = @(t, Y) [Y(3);
%                  Y(4);
%                  k/m*(Y(1)-xD(t))*(l0/sqrt((Y(2)-zD)^2+(Y(1)-xD(t))^2)-1)-rho*A*Cd/(2*m)*(Y(3)-w(Y(2)))^2*sign(Y(3)-w(Y(2)));
%                  k/m*(Y(2)-zD(t))*(l0/sqrt((Y(2)-zD(t))^2+(Y(1)-xD(t))^2)-1)-g];
% -----------------------------------------------------------------------------------------------------------------------
Ydot = @(t, y) ODEfcn(t, y,@xD,@zD,w,m,k,l0,rho,A,Cd,g);
[Y, R] = RK5solver(Ydot, tVec,IC);

%---------------------------------------------Matlab Solver for SelfCheck -------------------------------------
% [~, Y] = ode113(Ydot, tVec, IC); Y=Y';
%---------------------------------------------- Residum Estimation ----------------------------------------------
% sum(R,"all")
% max(R,[],'all')
%-----------------------------------------------------------------------------------------------------------------------
%% Q-5 Ploting

fig41 = plot_fig41(tVec, Y);                                % fig4.1: Package Motion x vx z vz
fig42 = plot_fig42(Y);                                          % fig4.2: Package Trajectory on XZ plane.

xDrone = zeros(1,length(tVec));
zDrone = zeros(1,length(tVec));
for i=1:length(tVec)
    xDrone(i) = xD(tVec(i));
    zDrone(i) = zD(tVec(i));
end

fig43 = plot_fig43(tVec, Y, xDrone, zDrone);  % fig4.3: Length of Rope as Fnc of Time
fig44 = plot_fig44(tVec, Y, xDrone, zDrone);  % fig4.4: Rope Angle as Fnc of time.
% -------------------------------- Add a Drone motion to fig4.1 --------------------------------
% plot_xDrone_yDrone(xDrone, zDrone,tVec,fig41)


%% stability

fig45=plot_fig45(tVec,Y,@RK5solver,Ydot, IC);  % fig4.5: Solution Stability shifted IC

save main_workspase.mat
disp(['time of run:  ', num2str(toc/60), 'min']); 
disp('');

%% asist functions:

function [fig11, fig12] = plot_fig11(t,z,analytical_sol,t0,tf,y0,yf,z10)
disp('☞plotting 2 figueres fig1.1, fig1.2:  Package Numerical vs Analytical Sol.')
fig11 = figure();
plot(t, z,'ob', 'MarkerSize',3, 'LineWidth',1.2,'DisplayName','Numerical RK4 Sol. z(t)')
legend
hold on
p=plot( t, analytical_sol, '--r', 'DisplayName','Analytical Sol.');
% legend('Numerical RK4 Sol.', 'Analytical Sol.')
xlim([0 13])
xlabel 't [sec]'
ylabel 'place magnitude z=z(t) of package [m]'
title 'fig1.1: Position on z axis of Package Numerical vs Analytical Sol. '
ylim([y0 yf])
plot(10, z10, 'kd','LineWidth',4, 'MarkerSize', 15,'DisplayName','z(t=10s)');
datatip(p,10,z10);
hold off
savefig('fig1.1.fig')

fig12 = figure();
subplot(2,1,1)
plot(t, z,'b','DisplayName','Numerical RK4 Sol. z(t)')
title 'Numerical Solution'
ylabel 'z(t) of package [m]'
xlabel 't [sec]'
xlim([t0 tf])
ylim([y0 yf])
box off
hold on
subplot(2,1,2)
plot( t, analytical_sol, 'r','DisplayName','Analytical Sol.')
title 'Analytical Solution'
ylabel 'z(t) of package [m]'
xlabel 't [sec]'
xlim([t0 tf])
ylim([y0 yf])
sgtitle 'fig1.2: Position on z axis of Package'
box off
hold off
savefig('fig1.2.fig')
end

function [f1, f2, f3] = plot_fig13(analytical_sol,solver, Ydot, IC, t0, tf,y0,yf)
disp('☞plotting 3 figueres fig1.3.1, fig1.3.2, fig1.3.3:  How h Influence on Convergence')
f1 = figure(); ax1 = axes ;hold(ax1,'on');   legend;
f2 = figure(); ax2 = axes ; hold(ax2,'on');  legend;
f3 = figure(); ax3 = axes ; hold(ax3,'on');  legend;
dis = @(x, y) abs(x-y);
dt=1e-4; tVec = (t0:dt:tf)';
name = sprintf('Analytical,  |R|=0');
plot(ax1,tVec, analytical_sol(tVec),'r','DisplayName', name,'LineWidth',0.9)
plot(ax2, tVec, analytical_sol(tVec),'r','DisplayName', name,'LineWidth',2.2)
for dt = [1e-3 1e-2 5e-2 1e-1]
    tVec = (t0:dt:tf)';
    [Y, R] = solver(Ydot, tVec,IC);
    z = Y(2,:);
    name = sprintf('h=%0.2g,  |R|≈%0.2g', dt, sum(R,"all"));
    plot(ax1, tVec, z, 'DisplayName', name)
    plot(ax2, tVec, z,'DisplayName', name)
    real_R = dis(z, analytical_sol(tVec)');
    plot(ax3, tVec, real_R, 'DisplayName', sprintf('h=%0.2g', dt))
    clear Y
end
for dt = [2.5e-1 5e-1]
    tVec = (t0:dt:tf)';
    [Y, R] = solver(Ydot, tVec,IC);
    z = Y(2,:);
    name = sprintf('h=%0.2g,  |R|≈%0.2g', dt, sum(R,"all"));
    plot(ax1, tVec, z,'--', 'DisplayName', name)
    plot(ax2, tVec, z,'--', 'DisplayName', name)
    real_R = dis(z, analytical_sol(tVec)');
    plot(ax3, tVec, real_R,'--', 'DisplayName', sprintf('h=%0.2g', dt))
    clear Y
end
tVec = (t0:0.661:tf)';
[Y, R] = solver(Ydot, tVec,IC);
z = Y(2,:);
name = sprintf('h=0.661,  |R|≈%0.2g', sum(R,"all"));
plot(ax1, tVec, z, 'DisplayName', name,'LineWidth',1)
plot(ax2, tVec, z, 'DisplayName', name,'LineWidth',1)
real_R = dis(z, analytical_sol(tVec)');
plot(ax3, tVec, real_R,'.', 'DisplayName', 'h=0.661')
title(ax1, "fig1.3.1: How h Influence on Convergence")
subtitle(ax1, 'Whole-time perspective')
title(ax2, "fig1.3.2: How h Influence on Convergence")
subtitle(ax2, 'Limited-time perspective, beginning of the motion')
title(ax3, "fig1.3.3: Distance form Analytical Sol, Worst h-es")
xlim(ax2,([0 10]))
ylim(ax1, [y0 yf])
ylim(ax2, [y0 yf])
% lg1=legend(ax1); lg2=legend(ax2); lg3=legend(ax3);

 savefig(f1,'fig1.3.1.fig')
 savefig(f2,'fig1.3.2.fig')
 savefig(f3,'fig1.3.3.fig')

end

function [f1, f2, f3] = plot_fig14(analytical_sol,solver, Ydot, IC, t0, tf)
disp('☞plotting 3 figueres fig1.4.1, fig1.4.2, fig1.4.3:   astimated and calculated Residum as fun of h')
f1= figure();
dis = @(x, y) abs(x-y);
h = 1e-2:5e-2:2;
R = zeros(1, length(h));
i = 1;
    for dt = h
        tVec = (t0:dt:tf)';
        [Y, ~] = solver(Ydot, tVec,IC);
        z = Y(2,:);
        R(i) = max(dis(z, analytical_sol(tVec)'),[],"all");
        i = i+1;
        clear Y
    end
semilogy(h, R)
title 'fig1.4.1: Exact residum vailad against analytical solution of z(t)'
subtitle 'h ∈ [10^{-2}, 2]'
xlabel 'h'
ylabel 'max R'
hold off
savefig('fig1.4.1.fig')


% h = 1e-4:5e-3:0.15;
h = logspace(-5,-1);
R = zeros(1, length(h));
i = 1;
    for dt = h
        tVec = (t0:dt:tf)';
        [Y, ~] = solver(Ydot, tVec,IC);
        z = Y(2,:);
        R(i) = max(dis(z, analytical_sol(tVec)'),[],"all");
        i = i+1;
        clear Y
    end

    f2 = NaN;
% f2 = figure();
% semilogy(h, R)
% title 'fig1.4.2: Exact residum vailad against analytical solution'
% subtitle 'h ∈ [10^-5, 10^-1] - semilog print'
% xlabel 'h'
% ylabel 'max R'
% savefig('fig1.4.2.fig')

f3 = figure();
loglog(h,R); grid on
title 'fig1.4.3: Exact residum vailad against analytical solution'
subtitle 'h ∈ [10^{-5}, 10^{-1}]  - loglog print'
xlabel 'h'
ylabel 'max R'
savefig('fig1.4.3.fig')
end

function z = zD(t)
l0   = 2;
% some more const.:
az0 = 3;                               % [m·s^-2]
t0  = 10;                               % [s]
zD0 = 10*l0;
dzD0 = double(0);

if t <= t0
     z = zD0 + dzD0.*t + az0/2.*t.^2;
elseif t<=2*t0
     z = zD0 + dzD0.*t  -  az0*t0^2 + 2*az0*t0.*t - az0/2.*t.^2;
else
     z = zD0 + dzD0*t +az0*t0^2;
end
end

function x = xD(t)
% given:
ax0 = 2;                            % [m·s^-2]
t1=7;               t2=11;      % [s]
xD0 = 0;
dxD0 = 0;

if            t<=t1
    x =  xD0 + dxD0 .* t + 0 ;
elseif    t<=t2
    x =  xD0 + dxD0 .* t + ax0/2.*(t-t1).^2;
 else
    x =  xD0 + dxD0 .* t  + ax0/2*(t2-t1)^2 + ax0*(t2-t1).*(t-t2);
 end
end

function dydt = ODEfcn(t,Y,xD,zD,w,m,k,l0,rho,A,Cd,g)
%--------------------------- Recall the State Vector Structure --------------------------------------------------
%                                   
%                                  Y  =  [Y(1) Y(2) Y(3) Y(4)]'  =  [x z vx vz]'
%
%--------------------------------------------------------------------------------------------------------------------------
x =   Y(1);
z =   Y(2);
vx = Y(3);
vz = Y(4);

X =  x - xD(t); 
Z =  z - zD(t);
V =  vx - w(z);
l = myhypot(X,Z);
r = mysgn(V);
FA = rho*A*Cd/2 .*(V.^2).*r;

dydt = zeros(4,1);
dydt(1) = vx;
dydt(2) = vz;
dydt(3) = k/m .* X .* (l0./l-1) - FA./m ;
dydt(4) = k/m .* Z .* (l0./l-1) - g;

    function out = mysgn(x)
        if x == 0
            out = 0;
        elseif x < 0
            out = -1;
        else
            out = 1;
        end
    end
    function out = myhypot(x, y)
out = sqrt(x.^2 + y.^2);
end

end

function f1 = plot_fig41(t,Y)
disp('☞plotting fig4.1: Package Motion')
x   = Y(1,:);
z   = Y(2,:);                              
vx  = Y(3,:);
vz  = Y(4,:);
f1 = figure();
plot(t,x,"LineWidth",1,DisplayName='$x(t)$')
grid on; hold on;
plot(t,z,'LineWidth',1,DisplayName='$z(t)$')
plot(t,vx,DisplayName='$\dot{x}(t)$')
plot(t,vz,DisplayName='$\dot{z}(t)$')
xlabel 'Time [sec]'
ylabel 'Function of t'
title 'fig4.1: Package Motion'
lgd = legend;
lgd.NumColumns = 2;
lgd.Interpreter = 'latex';
lgd.FontSize = 20;
savefig('fig4.1.fig')
end

function f1 = plot_fig42(Y)
disp( '☞plotting fig4.2: Package Trajectory on XZ plane.')
x = Y(1,:);
z = Y(2,:);
f1 = figure();
plot(x,z, 'LineWidth',2.5)
axis equal
xlabel 'x [m]'
ylabel 'z [m]'
grid on
line([0,0], ylim, 'Color', 'k', 'LineWidth', 1); % Draw line for Y axis.
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 1); % Draw line for X axis.
box off
title 'fig4.2: Package Trajectory on XZ plane.'
savefig('fig4.2.fig')
end

function f1 = plot_fig43(t, Y, xD, zD)
    disp('☞plotting fig4.3: Length of Rope as Func. of Time.')
    x = Y(1,:);
    z = Y(2,:);
    DeltaX = xD - x;
    DeltaZ = zD - z;
    l = 100.*hypot(DeltaX, DeltaZ);
    f1 = figure();
    plot(t, l)
    xlabel 'Time [sec]'
    ylabel 'l=l(t) [cm]'
    title 'fig4.3: Length of Rope as Func. of Time.'
    savefig('fig4.3.fig')
end

function f1 = plot_fig44(t, Y, xD, zD)
    name  = 'fig4.4: Rope Angle as Function of time.';
    disp(['☞ploting ', name])
    x = Y(1,:);
    z = Y(2,:);
    theta = -180/pi .* atan((x-xD)./(zD-z));
    f1 = figure();
    plot(t, theta,'LineWidth',1)
    i = 1;
    while x(i) < 100
        i=i+1; 
    end
    theta_max = max(abs(theta(1:i)));
    disp(['max{Θ(t) | t∈[0,100]} =  ', num2str(theta_max), ' [deg]'])
    xlabel 't [sec]'
    ylabel '\Theta=\Theta (t) [deg]'
    title(name)
    savefig('fig4.4.fig')
end

function fig = plot_xDrone_yDrone(xD, zD, t, fig)
disp('Adding Drone motion to fig4.1')
plot(fig, t, xD,DisplayName = '$x(t)_{Drone}$')
hold on
plot(fig, t, zD,DisplayName = '$z(t)_{Drone}$')
title 'fig5.5: Drone Motion as Function of Time'
xlabel 'Time [sec]'
ylabel 'y=f(t)'
savefig('figEx: Drone Motion.fig')
hold off
end

function f1 = plot_fig45(t,Y,solver,Ydot,IC)
disp('☞plotting fig4.5:   Solution Stability')
f1 = figure();                                                      % fig4.5: Solution Stability shifted IC
ax1 = subplot(2,2,1) ;hold(ax1,'on'); lgd1 = legend;
ax2 = subplot(2,2,2) ;hold(ax2,'on'); lgd2 = legend;
ax3 = subplot(2,2,3) ;hold(ax3,'on'); lgd3 = legend;
ax4 = subplot(2,2,4) ;hold(ax4,'on');  lgd4 = legend;
x   = Y(1,:);
z   = Y(2,:);                              
vx  = Y(3,:);
vz  = Y(4,:);
plot(ax1,t,x,'LineWidth',1.3,DisplayName='$IC\:\, \delta^*=0$');
plot(ax2,t,z,'LineWidth',1.3,DisplayName='$IC\:\, \delta^*=0$' );
plot(ax3,t,vx,'LineWidth',1.3,DisplayName='$IC\:\, \delta^*=0$'); 
plot(ax4,t,vz,'LineWidth',1.3,DisplayName='$IC\:\, \delta^*=0$'); 
xlabel(ax1,'t [s]');  ylabel(ax1, '$x(t) [m]$','Interpreter','latex');
xlabel(ax2,'t [s]');  ylabel(ax2, '$z(t) [m]$','Interpreter','latex');
xlabel(ax3,'t [s]');  ylabel(ax3, '$\dot{x}(t) [m/s]$','Interpreter','latex');
xlabel(ax4,'t [s]');  ylabel(ax4, '$\dot{z}(t) [m/s]$','Interpreter','latex');
lgd1.Interpreter = 'latex';lgd2.Interpreter = 'latex';
lgd3.Interpreter = 'latex';lgd4.Interpreter = 'latex';
sgtitle 'fig4.5: Solution Stability'
a = -0.5:0.2:0.5;
delta = a(a~=0); clear a;
delta = repmat(delta,4,1);
for i = 1:size(delta,2)
    [Y, ~] = solver(Ydot,t,IC+delta(:,i));
    plot_fig(t,Y, ax1,ax2,ax3,ax4,delta(1,i))
end

lgd1.Title.String = 'Shifted IC by:'; lgd1.Location="northwest";
lgd2.Title.String = 'Shifted IC by:'; lgd2.Location="southeast";
lgd3.Title.String = 'Shifted IC by:'; lgd3.Location="southeast";
lgd4.Title.String = 'Shifted IC by:';

savefig('fig4.5.fig')

function plot_fig(t,Y, ax1,ax2,ax3,ax4,d)
    x   = Y(1,:); %TODO: acurate legend names. don"t print all times
    z   = Y(2,:);                              
    vx  = Y(3,:);
    vz  = Y(4,:);
    str_d = num2str(d);
    n=10;  % skip plotting
    plot(ax1,t(1:n:end),x(1:n:end),'DisplayName',['$\delta^*=   \,$', str_d]);
    plot(ax2,t(1:n:end),z(1:n:end),'DisplayName',['$\delta^*=   \,$', str_d]); 
    plot(ax3,t(1:n:end),vx(1:n:end),'DisplayName',['$\delta^*= \,$', str_d]);
    plot(ax4,t(1:n:end),vz(1:n:end),'DisplayName',['$\delta^*= \,$', str_d]); 
end
end


