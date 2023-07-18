function dydt = section_d2_second_degenerated_problem()
% SECTION D.2: find the final velocity of the rocket for the second
% degenerated problem
% display the final velocity reached by the analytic solution, the numeric
% solution, and the difference between them

load basic_const.mat m0 dm0 ue g mf tb

md2 = @(t) m0-dm0*t;    % rocket mass for this section.                    
Td2 = dm0*ue;           % thrust

IC = [0 0 0 0]';
t = 0:1e-3:tb;

dydt = @(t,y) [ 0
                0 
                y(4) 
                Td2/md2(t)-g ];

[ynumeric,~] = RK5solver(dydt,t,IC);

vzEndAnalytic = ue*log(m0/mf)-g*tb
vzEndNumeric = ynumeric(end,4)

disp(['|V_{Analytical}-V_{Numerical}| = ' num2str(abs(vzEndAnalytic-vzEndNumeric))])

save section_d_well_defined_variable.mat vzEndAnalytic

end