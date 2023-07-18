function section_e2_plot_vx_vz_as_function_of_time(t,y)
% generate a graph which shows the velocity in the X and Z axis as a
% function of the time
% INPUT: t- time
%        y- the state vector we derived

figure
vx = y(:,2);
vz = y(:,4);
plot(t,vx,t,vz)
legend('V_x (t)','V_z (t)')
grid
xlabel 'Time [sec]'
title 'fig4: v_x & v_z'
end