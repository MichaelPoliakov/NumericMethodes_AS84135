%
%
% const:
c  = struct( ...
    'k'     , 50          , ...          [N/m]
    'l0'    , 2           , ...          [m]
    'rho'   , 1.225       , ...        
    'A'     , pi*(0.15)^2 , ...          [m^2]
    'C_d'   , 0.1         , ...
    'g'     , 9.81)       ; %            [m/s^2]
    


%QUESTION B:
format shortG
wind_u_10030cm = u_w(100.3)


%QUESTION D:






%%
function out = W(~)
m       =       2;          %[kg]       Mass
g       =       9.81;       %[N/kg]     Gravitation acc.
out = m.*g;
end

function out = F_el(l)
l0      =       2;          %[m]
k       =       50;         %[N/m]      String const.
out = (-1).*k.*(l-l0);
end

function out = l(xD, xP, zD, zP)
out = sqrt((xD-xP).^2 + (zD-zP).^2);
end

function out = F_A(z,Vx_p)
rho     =       1.225;      %[kg/m^3]   Air dencity
R       =       0.15;       %[m]        Radious
cD      =       0.1;        %[]         drag coefficient
A = pi*R^2;
V_pw = Vx_p - u_wind(z);
out = 0.5 * rho * A * cD .* V_pw.^2 .* sign(V_pw);
end


function out = u_w(z)
f = readmatrix("U_wind.txt",'Range','B:B');
f = f(2:end);
x = floor(z);
p = z - x;
if x==0 % f(0)=0
    out = p*f(x+1) + 0.5*p*(p-1)*(f(x+2)-2*f(x+1));
    return
end
out = f(x) + p*(f(x+1)-f(x)) + 0.5*p*(p-1)*(f(x+2)-2*f(x+1)+f(x));
end


function out1 = g(t, y)
k = 50;
l0 = 2;
rho = 1.225;
A = pi*(.15)^2;
C_D = .1;
m=2;
g = 9.81;

X = y(2)-y(1);
Z = y(5)-y(4);
l = sqrt((Z)^2 + (X)^2);
p = k * (l0 / l - 1);
out1 = p * Z - m*g;
out1 = m * out1;
end

function out2 = f(t, y)
k = 50;
l0 = 2;
rho = 1.225;
A = pi*(.15)^2;
C_D = .1;
m=2;
g = 9.81;

X = y(2)-y(1);
Z = y(5)-y(4);
l = sqrt((Z)^2 + (X)^2);
p = k * (l0 / l - 1);
V_pw = y(3)- u_w;
out2 = p * X + .5 * rho * A * C_d * V_pw^2 * sgn(V_pw);
out2 = m * out2;
end

