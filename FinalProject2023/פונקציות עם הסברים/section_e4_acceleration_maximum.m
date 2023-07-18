function section_e4_acceleration_maximum(t,y)
% SECTION E.4 calculate the max acceleration of the rocket
% INPUT: t- time
%        y- the state vector we derived

dydt = zeros(size(y));
for k = 2:numel(t)
    tau = t(k);
    cnt_dy = main_ode(tau, y(k,:));
    dydt(k,:) = reshape(cnt_dy,[1 4]);
end
ax = dydt(:,2); az = dydt(:,4);
acceleration_magnitude = sqrt(ax.^2 + az.^2);

max = acceleration_magnitude(1);
for k = 1:numel(t)
    if acceleration_magnitude(k)>max
        max = acceleration_magnitude(k);
        tmax = t(k);
    end
end
disp(['Maximus Acceleration Magnitude is: ' num2str(max) '[m·sec^-2]'])
disp(['Happen to be at: ' num2str(tmax) '[sec]'])