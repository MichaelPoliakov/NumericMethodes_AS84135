function out = dFdX(X)
% find the derivative of the squered function
% INPUT: X- the constant ratio Pe/P0
% OUTPUT: out- the derivative of the squered function
load pressure_const.mat kp
out =  - (X^((kp - 1)/kp - 1)*X^(2/kp)*(kp/2 + 1/2)^(2/(kp - 1))*(kp + 1))/kp - ...
(2*X^(2/kp - 1)*(X^((kp - 1)/kp) - 1)*(kp/2 + 1/2)^(2/(kp - 1))*(kp + 1))/(kp*(kp - 1));
end