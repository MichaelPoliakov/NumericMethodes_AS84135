function out = F(X)
% define the squered function
% finding a root of squered function is identical to finding the original
% function's root. We are going to find the root of the squered function to
% avoid negative arguments
% INPUT: X- the constant ratio Pe/P0
% OUTPUT: out- the squered function
load pressure_const.mat Ae At kp
out = -(X^(2/kp)*(X^((kp - 1)/kp) - 1)*(kp/2 + 1/2)^(2/(kp - 1))*(kp + 1))/(kp - 1) - (At/Ae)^2;
end