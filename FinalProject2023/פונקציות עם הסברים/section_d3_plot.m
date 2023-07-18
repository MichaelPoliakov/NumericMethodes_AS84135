function section_d3_plot(step,vz)
% SECTION D3: generate plot time step size vs vzEnd
% INPUT: step- step sizes
%        vz- velocity in vertical direction

load section_d_well_defined_variable.mat vzEndAnalytic 

figure(2)

semilogx(step,vz,'kO-')
title('fig2:   Numerical ODE Solver Accuracy as step-->0')

xlabel('Fixed Step Size [sec]');
ylabel('Calculated  Velocity in \it{z} direction at t_{f}=20s  [m/s]');
yline(vzEndAnalytic, 'r-.') 
set(gca,'XDir','reverse','XAxisLocation','bottom')
grid minor
xlim([1e-4 1])

exportgraphics(gcf,'fig2.pdf','Append',true,'contenttype','vector');
end