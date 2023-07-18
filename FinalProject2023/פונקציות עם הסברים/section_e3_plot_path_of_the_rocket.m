function p = section_e3_plot_path_of_the_rocket(y)
% generate a graph which presents the path of the rocket in Z axis as a
% function of the path in X axis
% INPUT: y- the state vector we derived

x = y(:,1); z = y(:,3);
figure;
p = plot(x,z);
title 'fig 5: Rocket Path in XZ-Plane'
axis equal
grid minor
xlabel 'x [m]'
ylabel 'z [m]'
end