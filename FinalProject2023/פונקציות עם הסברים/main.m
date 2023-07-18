%% Numerical Methods Project
%% 
% * |Ran|
%% 
% * |Michael|

matlab_set_up
%% Constatnts

basic_const_set_up
%% Section B
disp ' '
disp("SECTION B:")
section_b_print_specific_continuous_value_Mach
%% Section D
disp ' '
disp("SECTION D:")

% (D1)
section_d1_first_degenerated_problem

% (D2)
dydt = section_d2_second_degenerated_problem;

% (D3)
section_d3_convergence_test(dydt);
%% 
% 
%% Section E
% Note that $\frac{P_e }{P_0 }$ is constant. We will derive it, by Newton Raphson's 
% method. Once we have this ratio, we esealy may use it for deriving $P_e$ by 
% multiplication with $P_0$.
%
% Newton Raphson's method:
% 
% $$x_{n+1} :=x_n -\frac{f\left(x_n \right)}{f`\left(x_n \right)}$$
%
disp ' '
disp("SECTION E:")
% (E1)
% 
% Full problem Constants:
close all;
derive_Pe

% Initial Condition:
x0 =    0;
vx0 =   0;
z0 =    0;
vz0 =   0;

IC = [x0 vx0 z0 vz0];

% Time sampling
t = 0:1.6:45;

dydt = @main_ode;

[Y,~] = RK5solver(dydt,t,IC);
save
%%
[t,Y] = section_e1_cut_fly_at_final(t,Y);

section_e1_plot_xz_as_function_of_time(t,Y);

% (E.2)
section_e2_plot_vx_vz_as_function_of_time(t,Y);

% (E.3)
plot5 = section_e3_plot_path_of_the_rocket(Y);

% (E.4)
section_e4_acceleration_maximum(t,Y);

% (E.5)
section_e5_find_maximum_hight(Y);

% (E.6)
section_e6_rocket_fly_fime(t);

% (E.7)
section_e7_rocket_harizontal_Delta_x(Y);
%%  Convergence tests for Section E:

[steps, fly_time, delta_x] = section_e_convergence_tests(dydt,IC,45)