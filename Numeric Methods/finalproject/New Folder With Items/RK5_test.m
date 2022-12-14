%--------------------------------------------------------------------------
% 
%    Runge-Kutta 5th-order Integration
%
% Last modified:   2015/08/07   M. Mahooti
%
%--------------------------------------------------------------------------
clc
clear
format long g
% constants
GM  = 1;                   % gravitational coefficient
e   = 0.1;                 % eccentricity
Kep = [1, e ,0 ,0 ,0 ,0]'; % (a,e,i,Omega,omega,M)
% initial state
y_0 = State(GM, Kep, 0);
% header
fprintf( '\nRunge-Kutta 5th-order integration\n\n' );
% step size
h = 0.001; % [s]
% Initial values
t_0 = 0;
t_end = 86400; % end time [s]
Steps = t_end/h;
% integration from t=t_0 to t=t_end
for i=1:Steps
    y   = RK5(@deriv, t_0, y_0, h);
    y_0 = y;
    t_0 = t_0+h;
end
y_ref = State(GM, Kep, t_end); % Reference solution
fprintf('Accuracy   Digits\n');
fprintf('%7.3f',norm(y-y_ref));
fprintf('%9.2f\n',-log10(norm(y-y_ref)));