function [x, z, vx, vz] =  PredictorCorrectorEuler(t,x,z,vx,vz,xD,zD,w,h,m,k,l0,rho,A,Cd,g)

f1= @(t) F1(t,x,z,vx,vz,xD,zD,w,m,k,l0,rho,A,Cd,g);                                  
f2= @(t) F2(t,x,z,vx,vz,xD,zD,w,m,k,l0,rho,A,Cd,g);
f3= @(t) F3(t,x,z,vx,vz,xD,zD,w,m,k,l0,rho,A,Cd,g);
f4= @(t) F4(t,x,z,vx,vz,xD,zD,w,m,k,l0,rho,A,Cd,g);
% Y  = @(i) [x(i), z(i), vx(i), vz(i)];
% dY = @(i) [f1(i), f2(i), f3(i), f4(i)];

for i = 1:(length(t)-1)
    x(i+1)  = x(i)  + h*f1(i);
    z(i+1)  = z(i)  + h*f2(i);
    vx(i+1) = vx(i) + h*f3(i);
    vz(i+1) = vz(i) + h*f4(i);
% Y(i+1) = Y(i) + h .* dY(i);

% corrector step:
    x(i+1)  = x(i) + h/2*(f1(i)+f1(i+1));
    z(i+1)  = z(i) + h/2*(f2(i)+f2(i+1));
    vx(i+1)  = vx(i) + h/2*(f3(i)+f3(i+1));
    vz(i+1)  = vz(i) + h/2*(f4(i)+f4(i+1));
end

end
