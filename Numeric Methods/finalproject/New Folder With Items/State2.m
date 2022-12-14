%--------------------------------------------------------------------------
%
% State: Computes the satellite state vector from   
%          
%
% Inputs:
%   c         Constants that true throut the project, vector of constants.
%   u_w(z)    Wind v
%   D         Drone elements (x, z) with
%               x      x axis [m]
%               z      z axis [m]
%               
%   P         Package elements (x, z, vx, vz) with
%               x      x axis [m]
%               z      z axis [m]
%               vx     velocity on x axis [m/s]
%               vz     velocity on z axis [m/s]
%   dt        Time since epoch
% 
% Output:
%             State vector [P:(x,z,vx,vz)]
%
%   Reference:
%   O.Montenbruck, E. Gill, "Satellite Orbits - Models, Methods,
%   and Applications", Springer Verlag, Heidelberg, (2005)
%--------------------------------------------------------------------------
function Y = State2(D,P,c,u_w,t)
% Constanst elements at EOMs  
k   = c.k;    C_d  = c.C_d;
l0  = c.l0;   m    = c.m;
rho = c.rho;  g    = c.g;
A   = c.A;

Z = @(t) P.z(t) - D.z(t);  X = @(t) P.z(t) - D.z(t);
alpha = -1 + l0 / sqrt(Z^2 + X^2);
alpha = k * alpha;
P_dVx = @(t) m*(alpha(t)*X(t) +0.5*rho*A*C_d*(P.vx(t)-u_w(t))^3/abs(P.vx(t)-u_w(t)));
P_dVz = @(t) m*(alpha(t)*Z(t) -m*g);

Y = [P.vx;
     P.vz;
     P_dVx;
     P_dVz];
return

% % Mean anomaly  
% if (dt==0)
%   M = M0;
% else
%   n = sqrt(gm/(a*a*a));
%   M = M0 +n*dt;
% end
% % Eccentric anomaly  
% E  = EccAnom(M,e);
% cosE = cos(E);
% sinE = sin(E);
% % Perifocal coordinates
% fac = sqrt((1-e)*(1+e));
% R = a*(1-e*cosE);  % Distance
% V = sqrt(gm*a)/R;    % Velocity
% r = [a*(cosE-e), a*fac*sinE , 0]';
% v = [-V*sinE   , +V*fac*cosE, 0]';
% % Transformation to reference system (Gaussian vectors)
% PQW = R_z(-Omega) * R_x(-i) * R_z(-omega);
% r = PQW*r;
% v = PQW*v;
% Y = [r;v];
end