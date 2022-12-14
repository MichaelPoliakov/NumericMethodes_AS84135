function out = F4(t,x,z,~,~,xD,zD,~,m,k,l0,~,~,~,g)  
x = x(t);
z = z(t);
xD = xD(t);
zD = zD(t);
out = k/m*(z-zD)*(-1+l0/sqrt((z-zD)^2+(x-xD)^2)) - g;
end