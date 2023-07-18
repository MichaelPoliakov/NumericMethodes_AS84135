function section_e1_plot_xz_as_function_of_time(t,y)
% SECTION E.1 plot a graph which shows the rocket's movement in the X and Z
% axis as a function of the time
% INPUT: t- time
%        y- the state vector we derived

load basic_const.mat h
x = y(:,1);
z = y(:,3);

plot(t,x,t,z)

title 'fig3: \it{x}, \it{z} as function of Time'
yline(-h,'--','color',[.5 .5 .5])
legend('x(t) [m/s]','z(t) [m/s]','Location','best')
xlabel 't [sec]'
box off

end