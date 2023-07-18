function section_b_print_specific_continuous_value_Mach()
% find the drag coefficient for a given Mach number
Mach = 0.820;
CD = constructCD;
result = CD(Mach);
fprintf('Drag Coefficent at Mach %.3g is: \n CD(Mach=%g) = %g\n',Mach,Mach,result)
save section_e_well_defined_functions.mat CD
disp('Save CD(Mach) to New file: section_e_well_defined_functions.mat')
end