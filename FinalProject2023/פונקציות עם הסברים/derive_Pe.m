function derive_Pe()
%calculate the fluid pressure exits the engine
%OUTPUT: Pe

load pressure_const.mat tb Pa

Po = @(t)       (t <=tb).*(Pa + 10^7*(atan(t)-atan(t-tb)-1.5)) + ...
                (t > tb).* Pa;   % Cylinder pressure

initial_guess = 0.1;
tolernece = 1e-9;
max_iterations = 1e5;

[X,~] = rootNewtonRaphson(@F,@dFdX,initial_guess,tolernece,max_iterations);

Pe = @(t) X .* Po(t);

save section_e_well_defined_functions.mat Pe -append
disp('Saves Pe(t) to section_e_well_defined_functions.mat')
end