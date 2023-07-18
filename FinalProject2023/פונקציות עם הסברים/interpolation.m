function yy = interpolation(x, y, xx)
% commit linear interpolation: find an intermediate value in a function
% which contains discrete points
% INPUT: x- X values of the function
%        y- y values of the function
%        xx- X value of the required point
% OUTPUT: yy- Y value of the required point
h = x(2)-x(1); % h is constant

x0 = floor(xx./h).*h;

t0 = int32((x0 / h) + 1);
t1 = t0 + 1;

delta_y = y(t1) - y(t0);
yy = y(t0) + (xx-x(t0))/h * delta_y;
end