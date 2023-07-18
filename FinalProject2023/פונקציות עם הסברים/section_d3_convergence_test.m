function Y = section_d3_convergence_test(f)
% SECTION D3: Accuracy of V in z direction as ODE solver time step << 1
% INPUT: f- the derivative of unknown function Y
% plot out time step size vs vzEnd

load basic_const.mat tb

step = logspace(-4,1);
vzEndList = zeros(size(step));

t = @(k) 0:step(k):tb;
initial_conditions = zeros(4,1);

for k = 1:numel(step)
    [Y,~] = RK5solver(f,t(k),initial_conditions);
    vzEndList(k) = Y(end,4);
end
save section_d3_data.mat vzEndList

section_d3_plot(step,vzEndList);

end