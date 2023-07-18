function [step,flyingTimeList, deltaXList] = section_e_convergence_tests(f,IC,t_final)
% SECTION E: find the horizontal distance and the flight time of the rocket
% as a function of the step size
% INPUT: f- the function
%           IC- initial conditions
%           t_final

step = logspace(-1,0.5,6);

flyingTimeList = zeros(size(step));
deltaXList     = zeros(size(step));

disp(' ');
disp("SECTION E TESTS:")

disp('step-h V resulted delta x:')
figure
hold on
for k = 1:numel(step)
    t = 0:step(k):t_final;
    [Y,~] = RK5solver(f,t,IC);
    [t,Y] = section_e1_cut_fly_at_final(t,Y);
    h_str = num2str(step(k));
    
    plot(Y(:,1)./10^3,Y(:,3)./10^3,'DisplayName',['\it{h} = ' h_str])
    disp(['h = ' h_str])
    flyingTimeList(k) = t(end);
    deltaXList(k) = section_e7_rocket_harizontal_Delta_x(Y);
end
legend
xlabel 'x [km]'
ylabel 'z [km]'
axis equal
box off
title 'Influence of the step size (h) on horizontal distance \Delta\it{X}'
hold off

figure
semilogx(step,flyingTimeList,'LineStyle','none','Marker','o','MarkerSize',12)
grid minor
xlabel('step size [s]')
ylabel('flying time [s]')
title('Influence of the step size (h) on flying time \Delta t')
box off
set(gca,'XDir','reverse')

disp('step-h V resulted flying time:')
end