function out = F3(t,x,z,vx,~,xD,zD,w,m,k,l0,rho,A,Cd,~)
a = 0.5*rho*A*Cd;
x = x(t);
z = z(t);
vx = vx(t);
xD = xD(t);
zD = zD(t);
w = w(z); 
out = m^-1 * (k*(x-xD)*(-1+l0/sqrt((z-zD)^2+(x-xD)^2)) + a*(vx-w)^2*sign(vx-w));
end