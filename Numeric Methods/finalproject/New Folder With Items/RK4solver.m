function [x,z,vx,vz] = RK4solver(t,x,z,vx,vz,xD,zD,w,h,m,k,l0,rho,A,Cd,g)  
a = 0.5*rho*A*Cd/m;

f1= @(~,~,~,vx,~)          vx;                                  
f2= @(~,~,~,~,vz)          vz;
f3= @(~,x,z,vx,vz,xD,zD,w) k/m*(x-xD)*(-1+l0/sqrt((z-zD)^2+(x-xD)^2)) + a*(vx-w)^2*sign(vx-w);
f4= @(~,x,z,~,vz,xD,zD)    k/m*(z-zD)*(-1+l0/sqrt((z-zD)^2+(x-xD)^2)) - g;


for i=1:(length(x)-1)                              % calculation loop
    a_1 = h*f1( t(i),           x(i),             z(i),             vx(i),             vz(i)                                     );
    b_1 = h*f2( t(i),           x(i),             z(i),             vx(i),             vz(i)                                     );
    c_1 = h*f3( t(i),           x(i),             z(i),             vx(i),             vz(i),            xD(i),  zD(i),  w(z(i)) );
    d_1 = h*f4( t(i),           x(i),             z(i),             vx(i),             vz(i),            xD(i),  zD(i)           );
    
    a_2 = h*f1( t(i)+0.5*h,     x(i)+0.5*a_1,     z(i)+0.5*b_1,     vx(i)+0.5*c_1,     vz(i)+0.5*d_1                              );
    b_2 = h*f2( t(i)+0.5*h,     x(i)+0.5*a_1,     z(i)+0.5*b_1,     vx(i)+0.5*c_1,     vz(i)+0.5*d_1                              );
    c_2 = h*f3( t(i)+0.5*h,     x(i)+0.5*a_1,     z(i)+0.5*b_1,     vx(i)+0.5*c_1,     vz(i)+0.5*d_1,     xD(i),  zD(i),  w(z(i)) );
    d_2 = h*f4( t(i)+0.5*h,     x(i)+0.5*a_1,     z(i)+0.5*b_1,     vx(i)+0.5*c_1,     vz(i)+0.5*d_1,     xD(i),  zD(i)           );

    a_3 = h*f1((t(i)+0.5*h),   (x(i)+0.5*a_2),   (z(i)+0.5*b_2),    (vx(i)+0.5*c_2),  (vz(i)+0.5*d_2)                             );
    b_3 = h*f2((t(i)+0.5*h),   (x(i)+0.5*a_2),   (z(i)+0.5*b_2),    (vx(i)+0.5*c_2),  (vz(i)+0.5*d_2)                             );
    c_3 = h*f3( t(i)+0.5*h,    (x(i)+0.5*a_2),   (z(i)+0.5*b_2),    (vx(i)+0.5*c_2),  (vz(i)+0.5*d_2),    xD(i),  zD(i),  w(z(i)) );
    d_3 = h*f4((t(i)+0.5*h),   (x(i)+0.5*a_2),   (z(i)+0.5*b_2),    (vx(i)+0.5*c_2),  (vz(i)+0.5*d_2),    xD(i),  zD(i)           );
       
    a_4 = h*f1((t(i)+h),       (x(i)+a_3),       (z(i)+b_3),        (vx(i)+c_3),      (vz(i)+d_3)                               );
    b_4 = h*f2((t(i)+h),       (x(i)+a_3),       (z(i)+b_3),        (vx(i)+c_3),      (vz(i)+d_3)                               );
    c_4 = h*f3((t(i)+h),       (x(i)+a_3),       (z(i)+b_3),        (vx(i)+c_3),      (vz(i)+d_3),      xD(i),  zD(i),  w(z(i)) );
    d_4 = h*f4((t(i)+h),       (x(i)+a_3),       (z(i)+b_3),        (vx(i)+c_3),      (vz(i)+d_3),      xD(i),  zD(i)           );
    
    x(i+1)  = x(i)  + (1/6)*(a_1+2*a_2+2*a_3+a_4);  % main equation
    z(i+1)  = z(i)  + (1/6)*(b_1+2*b_2+2*b_3+b_4);  % main equation
    vx(i+1) = vx(i) + (1/6)*(c_1+2*c_2+2*c_3+c_4);  % main equation
    vz(i+1) = vz(i) + (1/6)*(d_1+2*d_2+2*d_3+d_4);  % main equation

end



end
