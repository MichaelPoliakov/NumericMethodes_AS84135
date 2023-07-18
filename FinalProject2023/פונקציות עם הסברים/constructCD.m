function CD_continuous = constructCD()
% calculate the drag coefficient for a certain mach
% OUTPUT: drag coefficient
Mach_CD_table = load("CD_vs_M.txt");
Mach = Mach_CD_table(:,1);
CD = Mach_CD_table(:,2);
CD_continuous = @(m) interpolation(Mach,CD,m);
end